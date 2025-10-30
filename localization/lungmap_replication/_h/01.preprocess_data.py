## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.
import os
import numpy as np
import session_info
import scanpy as sc
import harmonypy as hm
from pyhere import here
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

def _sanitize_var_for_h5ad(adata):
    if isinstance(adata.var.index.dtype, CategoricalDtype):
        adata.var.index = pd.Index(adata.var.index.astype(str))
    else:
        adata.var.index = pd.Index(adata.var.index.astype(str))

    adata.var.index.name = None
    adata.var_names_make_unique()
    return adata


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
    plt.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{umap_file}.pdf", bbox_inches="tight")
    plt.close()
    return adata


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""
    # Run harmony
    batch_vars = ["donor", "batch"]
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
    if not isinstance(adata.obs["ann_level_4"].dtype, CategoricalDtype):
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

    # Preprocess data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, batch_key='study'
    )
    return adata


def prepare_data(query_adata, ref_adata):
    # Align by gene symbols
    ref_feat = ref_adata.var['feature_name'].astype(str)
    qry_names = query_adata.var_names.astype(str)
    common = pd.Index(ref_feat).intersection(pd.Index(qry_names))

    # Subset both to common genes
    ref = ref_adata[:, ref_feat.isin(common)].copy()
    qry = query_adata[:, query_adata.var_names.isin(common)].copy()

    # Drop duplicate feature_names in reference
    dedup_mask = ~ref.var['feature_name'].duplicated(keep='first')
    ref = ref[:, dedup_mask].copy()

    # Set reference var_names to 'feature_name'
    ref.var_names = ref.var['feature_name'].astype(str).values
    if 'feature_name' in ref.var.columns:
        ref.var = ref.var.rename(columns={'feature_name': 'gene_name'})
        
    # Subset both by the same HVG set (order-preserving)
    hvg_genes = ref.var_names[ref.var['highly_variable'].values]
    hvg_genes = [g for g in hvg_genes if g in qry.var_names]  # safety
    ref_hvg = ref[:, hvg_genes].copy()
    query_hvg = qry[:, hvg_genes].copy()

    # Sanity check same order
    assert (ref_hvg.var_names == query_hvg.var_names).all()

    return ref_hvg, query_hvg


def main():
    # Load data
    query_adata = process_query_data()
    ref_adata = load_reference()
    if "counts" not in ref_adata.layers:
        ref_adata.layers["counts"] = ref_adata.X

    # Prepare data
    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)

    # Sanitize .var for HDF5
    ref_hvg = _sanitize_var_for_h5ad(ref_hvg)
    query_hvg = _sanitize_var_for_h5ad(query_hvg)

    # Write files
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    # Reproducibility info
    session_info.show()


if __name__ == "__main__":
    main()

