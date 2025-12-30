## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.
import numpy as np
import session_info
import scanpy as sc
import harmonypy as hm
from pyhere import here
from pathlib import Path
import matplotlib.pyplot as plt

def set_seed(SEED):
    np.random.seed(13)
    sc.settings.seed = 13


def subset_data(subset_key: str = "cell_type",
                subset_value: str = "Pericyte"):
    """Subset lung single-cell data for label transfer."""
    # Load AnnData
    input_path = Path(here('disease_association/ipf_analysis',
                           '_m/ipf_dataset.h5ad'))
    adata = sc.read_h5ad(input_path)
    
    # Ensure count layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if subset_key not in adata.obs:
        raise KeyError(f"{subset_key} not found in adata.obs")
    mask = adata.obs[subset_key].eq(subset_value)

    return adata[mask].copy()


def check_data(adata, outdir: Path = Path("qc_plots")):
    # Preprocess and dimensionality reduction
    sc.pp.normalize_total(adata, layer="counts")
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
    
    # Plot & Save
    outdir.mkdir(parents=True, exist_ok=True)

    # Create and save the plot
    umap_file = outdir / "umap_overview"
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    fig = sc.pl.umap(adata, color=["Library_Identity", "patient", "disease"],
                     wspace=0.4, ncols=1, show=False, return_fig=True)

    # Save manually for full control
    fig.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{umap_file}.pdf", bbox_inches="tight")
    plt.close(fig)
    return adata


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""    
    # Run harmony
    batch_vars = ["patient"]
    meta = adata.obs[batch_vars].copy()
        
    for col in batch_vars:
        meta[col] = meta[col].astype("category").astype(str)

    harmony_embed = hm.run_harmony(
        adata.obsm["X_pca"], meta, vars_use=batch_vars,
        max_iter_harmony=max_iter, epsilon_harmony=1e-5, 
        random_state=seed
    )
    adata.obsm["X_pca_harmony"] = harmony_embed.Z_corr.T

    return adata


def process_query_data():
    # Preprocessing query data
    adata = subset_data()
    adata = check_data(adata, outdir=Path("qc_plots"))
    adata = preprocess_data(adata)
    return adata


def load_reference():
    input_path = Path(here("localization/pericyte_analysis/_m",
                           "pericyte_with_embeddings.h5ad"))
    adata = sc.read_h5ad(input_path)
    adata.obs["celltype"] = adata.obs["leiden_pericytes"]
    mask = adata.var["feature_name"].notna() & (adata.var["feature_name"] != "")

    adata.var["new_names"] = adata.var_names
    adata.var.loc[mask, "new_names"] = adata.var.loc[mask, "feature_name"]
    adata.var_names = adata.var["new_names"]
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    return adata


def prepare_data(query_adata, ref_adata):
    # Ensure common genes
    common_genes = ref_adata.var_names.intersection(query_adata.var_names)
    ref_adata = ref_adata[:, common_genes].copy()
    query_adata = query_adata[:, common_genes].copy()

    # Preprocessing
    if "highly_variable" not in ref_adata.var:
        raise ValueError("Reference adata lacks highly_variable annotation")

    hvgs = ref_adata.var["highly_variable"].values.astype(bool)
    ref_hvg   = ref_adata[:, hvgs].copy()
    query_hvg = query_adata[:, hvgs].copy()
    return ref_hvg, query_hvg


def main():
    set_seed(13)

    # Load data
    query_adata = process_query_data()
    ref_adata = load_reference()
    if "counts" not in ref_adata.layers:
        ref_adata.layers["counts"] = ref_adata.X.copy()

    # Prepare data
    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)
    
    # Save preprocessed objects for GPU step
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    # Reproducibility info
    session_info.show()


if __name__ == "__main__":
    main()
