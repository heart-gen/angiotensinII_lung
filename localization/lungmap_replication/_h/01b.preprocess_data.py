## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.
import os
import numpy as np
import pandas as pd
import session_info
import scanpy as sc
import harmonypy as hm
from pathlib import Path
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype

def _to_float32(data):
    """Downcast dense or sparse arrays to float32 when possible."""
    if data is None:
        return data

    if hasattr(data, "dtype") and data.dtype != np.float32:
        return data.astype(np.float32, copy=False)

    return data


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
    adata.X = _to_float32(adata.X)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata)
    adata.obsm["X_pca"] = _to_float32(adata.obsm.get("X_pca"))
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
    adata.obsm["X_pca_harmony"] = _to_float32(harmony_embed.Z_corr.T)

    return adata


def process_query_data():
    # Preprocessing query data
    adata = load_data()
    adata = check_data(adata, outdir="qc_plots")
    adata = preprocess_data(adata)
    return adata


def prepare_data(query_adata, ref_adata):
    # Align by gene symbols
    ref_feat = ref_adata.var['feature_name'].astype(str)
    qry_names = pd.Index(query_adata.var_names.astype(str))

    # Build a mask for features present in both datasets
    common_mask_ref = ref_feat.isin(qry_names)
    if not common_mask_ref.any():
        raise ValueError("No overlapping features found between reference and query data.")

    # Drop duplicate feature names while preserving the first occurrence
    common_ref_features = ref_feat[common_mask_ref]
    dedup_mask = ~common_ref_features.duplicated(keep='first')
    final_ref_mask = common_mask_ref.copy()
    final_ref_mask[common_mask_ref] = dedup_mask

    # Subset the reference once using the combined mask
    ref = ref_adata[:, final_ref_mask].copy()

    # Set reference var_names to 'feature_name'
    ref.var_names = ref.var['feature_name'].astype(str).values
    if 'feature_name' in ref.var.columns:
        ref.var = ref.var.rename(columns={'feature_name': 'gene_name'})

    # Determine HVGs that are shared with the query dataset
    hvg_genes = ref.var_names[ref.var['highly_variable'].values]
    hvg_genes = [gene for gene in hvg_genes if gene in qry_names]
    if not hvg_genes:
        raise ValueError("No shared highly variable genes found between reference and query data.")

    # Slice without creating intermediate dense copies
    ref_hvg = ref[:, hvg_genes].copy()
    query_hvg = query_adata[:, hvg_genes].copy()

    ref_hvg.X = _to_float32(ref_hvg.X)
    query_hvg.X = _to_float32(query_hvg.X)

    # Sanity check same order
    assert (ref_hvg.var_names == query_hvg.var_names).all()

    return ref_hvg, query_hvg


def main():
    # Load data
    query_adata = process_query_data()
    ref_adata = sc.read_h5ad(Path("ref_preprocessed.h5ad"))
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

