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


def _load_process_query_module():
    module_path = Path(__file__).with_name("01b.process_query.py")
    spec = importlib.util.spec_from_file_location("process_query", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def prepare_data(query_adata, ref_adata):
    # Align by gene symbols
    ref_feat = ref_adata.var["feature_name"].astype(str)
    qry_names = pd.Index(query_adata.var_names.astype(str))

    # Build a mask for features present in both datasets
    common_mask_ref = ref_feat.isin(qry_names)
    if not common_mask_ref.any():
        raise ValueError("No overlapping features found between reference and query data.")

    # Drop duplicate feature names while preserving the first occurrence
    common_ref_features = ref_feat[common_mask_ref]
    dedup_mask = ~common_ref_features.duplicated(keep="first")
    final_ref_mask = common_mask_ref.copy()
    final_ref_mask[common_mask_ref] = dedup_mask

    # Subset the reference once using the combined mask
    ref = ref_adata[:, final_ref_mask].copy()

    # Set reference var_names to 'feature_name'
    ref.var_names = ref.var["feature_name"].astype(str).values
    if "feature_name" in ref.var.columns:
        ref.var = ref.var.rename(columns={"feature_name": "gene_name"})

    # Determine HVGs that are shared with the query dataset
    ref_hvg_mask = ref.var["highly_variable"].to_numpy()
    if ref_hvg_mask.dtype != bool:
        ref_hvg_mask = ref_hvg_mask.astype(bool, copy=False)

    shared_mask = pd.Index(ref.var_names).isin(qry_names)
    ref_hvg_mask &= shared_mask.to_numpy()
    if not ref_hvg_mask.any():
        raise ValueError("No shared highly variable genes found between reference and query data.")

    # Slice without creating intermediate dense copies
    ref_hvg = ref[:, ref_hvg_mask].copy()
    del ref

    shared_genes = ref_hvg.var_names
    query_hvg = query_adata[:, shared_genes].copy()

    ref_hvg.X = _to_float32(ref_hvg.X)
    query_hvg.X = _to_float32(query_hvg.X)

    # Sanity check same order
    assert (ref_hvg.var_names == query_hvg.var_names).all()

    gc.collect()

    return ref_hvg, query_hvg


def main():
    process_query = _load_process_query_module()
    query_adata = process_query.process_query_data()

    ref_adata = sc.read_h5ad(Path("ref_preprocessed.h5ad"))
    if "counts" not in ref_adata.layers:
        ref_adata.layers["counts"] = ref_adata.X

    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)

    ref_hvg = _sanitize_var_for_h5ad(ref_hvg)
    query_hvg = _sanitize_var_for_h5ad(query_hvg)

    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    session_info.show()


if __name__ == "__main__":
    main()
