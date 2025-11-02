"""Reference/query alignment utilities for lungmap label transfer."""

import gc
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
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


def prepare_data(query_adata, ref_adata):
    # Enrsure string
    ref_adata.var_names = ref_adata.var_names.astype(str)
    query_adata.var_names = query_adata.var_names.astype(str)

    # Align on gene symbols
    ref_feat  = ref_adata.var["feature_name"].astype(str)
    qry_names = pd.Index(query_adata.var_names)
    
    # Build final mask
    final_ref_mask = ref_feat.isin(qry_names) & ~ref_feat.duplicated(keep="first")

    # Subset refernece once
    ref = ref_adata[:, final_ref_mask.to_numpy()].copy()

    # Set reference var_names to 'feature_name'
    ref.var_names = ref.var["feature_name"].astype(str).values
    if "feature_name" in ref.var.columns:
        ref.var = ref.var.rename(columns={"feature_name": "gene_name"})

    # Determine HVGs that are shared with the query dataset
    hvg = ref.var_names[ref.var["highly_variable"].astype(bool).values]
    hvg = pd.Index(hvg).intersection(qry_names)

    if len(hvg) == 0:
        raise ValueError("No shared highly variable genes found between reference and query data.")

    # Slice without creating intermediate dense copies
    ref_hvg   = ref[:, hvg].copy()
    query_hvg = query_adata[:, hvg].copy()

    ref_hvg.X = _to_float32(ref_hvg.X)
    query_hvg.X = _to_float32(query_hvg.X)
    assert (ref_hvg.var_names == query_hvg.var_names).all()

    return ref_hvg, query_hvg


def main():
    # Load data
    query_adata = sc.read_h5ad(Path("lungmap_query.h5ad"))
    ref_adata = sc.read_h5ad(Path("ref_preprocessed.h5ad"))

    # Prepare for trainning
    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)

    ref_hvg = _sanitize_var_for_h5ad(ref_hvg)
    query_hvg = _sanitize_var_for_h5ad(query_hvg)

    # Write to file
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("../_m/query_hvg.h5ad", compression="gzip")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
