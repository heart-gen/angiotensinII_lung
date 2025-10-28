## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.
import os
import numpy as np
import session_info
import pandas as pd
import scanpy as sc
import harmonypy as hm
from pyhere import here
from pathlib import Path
import matplotlib.pyplot as plt

def load_data():
    """Load data for label transfer."""
    # Load AnnData
    input_path = Path('lungmap_dataset.h5ad')
    adata = sc.read_h5ad(Path(input_path))
    # Ensure count layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    return adata


def check_data(adata, outdir="qc_plots"):
    # Preprocess and dimensionality reduction
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Plot & Save
    os.makedirs(outdir, exist_ok=True)

    # Define the filename base
    umap_file = os.path.join(outdir, "umap_overview")

    # Create and save the plot
    sc.pl.umap(adata, color=["donor", "age", "sex", "batch"],
               wspace=0.4, ncols=1, save=None, show=False)

    # Save manually for full control
    plt.tight_layout()
    plt.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{umap_file}.pdf", bbox_inches="tight")
    plt.close()
    return adata


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""
    # Run harmony
    batch_vars = ["donor"]
    meta = adata.obs[batch_vars].copy()

    for col in batch_vars:
        meta[col] = meta[col].astype("category")

    harmony_embed = hm.run_harmony(
        adata.obsm["X_pca"], meta, vars_use=batch_vars,
        max_iter_harmony=max_iter, epsilon_harmony=1e-5,
        random_state=seed
    )
    adata.obsm["X_pca_harmony"] = harmony_embed.Z_corr.T

    return adata


def process_query_data():
    # Preprocessing query data
    adata = load_data()
    adata = check_data(adata, outdir="qc_plots")
    adata = preprocess_data(adata)
    return adata


def load_reference():
    input_path = Path(here("inputs/hlca/_m/hlca_core.h5ad"))
    adata = sc.read_h5ad(input_path)
    # Ensure count layer
    if "counts" not in adata.layers:
        if "soupX" in adata.layers:
            adata.layers["counts"] = adata.layers["soupX"]
        elif adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
        else:
            raise ValueError("No suitable count layer found (expected 'counts', 'soupX', or .raw).")

    # Define categories to add
    required_cats = ["Vascular smooth muscle", "Mesothelium", "Myofibroblasts"]
    if not pd.api.types.is_categorical_dtype(adata.obs["ann_level_4"]):
	adata.obs["ann_level_4"] = adata.obs["ann_level_4"].astype("category")

    for cat in required_cats:
        if cat not in adata.obs["ann_level_4"].cat.categories:
            adata.obs["ann_level_4"] = adata.obs["ann_level_4"].cat.add_categories([cat])

    # Handle annotation issues
    vsm_clusters = {"Smooth muscle", "Smooth muscle FAM83D+",
                    "SM activated stress response"}
    fixes = {
        "Vascular smooth muscle": vsm_clusters,
        "Mesothelium": {"Mesothelium"},
        "Myofibroblasts": {"Myofibroblasts"},
    }

    for label, fine_clusters in fixes.items():
        cond = (
            adata.obs["ann_finest_level"].isin(fine_clusters)
            & (adata.obs["ann_level_4"].isna() | (adata.obs["ann_level_4"].isin(["None", ""])))
        )
        adata.obs.loc[cond, "ann_level_4"] = label

    # Update annotation columns
    obs_map = {
        "subclusters": "ann_finest_level",
        "cell_type": "ann_level_4",
        "clusters": "ann_level_4",
        "compartment": "ann_level_1",
        "patient": "donor_id",
    }
    for new, old in obs_map.items():
        if old in adata.obs:
            adata.obs[new] = adata.obs[old]
        else:
            print(f"Warning: '{old}' not found in obs; skipping {new} mapping.")

    # Filter studies (>= 20 cells)
    if "study" in adata.obs:
        study_counts = adata.obs["study"].value_counts()
        valid_studies = study_counts[study_counts >= 20].index
        adata = adata[adata.obs["study"].isin(valid_studies)].copy()
    else:
        print("Warning: Variable 'study' not found in observation metadata; skipping study filter.")

    return adata


def prepare_data(query_adata, ref_adata):
    # Ensure common genes
    common_genes = ref_adata.var_names.intersection(query_adata.var_names)
    ref_adata = ref_adata[:, common_genes].copy()
    query_adata = query_adata[:, common_genes].copy()

    # Preprocessing
    hvgs = ref_adata.var.highly_variable
    hvgs = hvgs.intersection(query_adata.var_names)
    ref_hvg = ref_adata[:, hvgs].copy()
    query_hvg = query_adata[:, hvgs].copy()
    return ref_hvg, query_hvg


def main():
    # Load data
    query_adata = process_query_data()
    ref_adata = load_reference()
    if "counts" not in ref_adata.layers:
        ref_adata.layers["counts"] = ref_adata.X

    # Prepare data
    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)

    # Save preprocessed objects for GPU step
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    # Reproducibility info
    session_info.show()


if __name__ == "__main__":
    main()

