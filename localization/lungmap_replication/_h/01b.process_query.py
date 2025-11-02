"""Query-side preprocessing utilities for lungmap label transfer."""

import gc, os
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import harmonypy as hm
from pathlib import Path
import matplotlib.pyplot as plt


def _to_float32(data):
    """Downcast dense or sparse arrays to float32 when possible."""
    if data is None:
        return data

    if hasattr(data, "dtype") and data.dtype != np.float32:
        return data.astype(np.float32, copy=False)

    return data


def load_data():
    """Load the query AnnData object."""
    input_path = Path("lungmap_dataset.h5ad")
    adata = sc.read_h5ad(input_path)

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    return adata


def check_data(adata, outdir: str = "qc_plots"):
    """Run basic QC, dimensionality reduction, and UMAP plotting."""
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.X = _to_float32(adata.X)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata)
    adata.obsm["X_pca"] = _to_float32(adata.obsm.get("X_pca"))
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    os.makedirs(outdir, exist_ok=True)
    umap_file = Path(outdir) / "umap_overview"

    sc.pl.umap(
        adata,
        color=["donor", "age", "sex", "batch"],
        wspace=0.4,
        ncols=1,
        save=None,
        show=False,
    )

    plt.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{umap_file}.pdf", bbox_inches="tight")
    plt.close()

    return adata


def preprocess_data(adata, *, max_iter: int = 30, seed: int = 13):
    """Batch-correct query data with Harmony."""
    batch_vars = ["donor", "batch", "age", "sex"]
    for col in batch_vars:
        if not pd.api.types.is_categorical_dtype(adata.obs[col]):
            adata.obs[col] = adata.obs[col].astype("category")

    meta = adata.obs.loc[:, batch_vars]

    harmony_embed = hm.run_harmony(
        adata.obsm["X_pca"],
        meta,
        vars_use=batch_vars,
        max_iter_harmony=max_iter,
        epsilon_harmony=1e-5,
        random_state=seed,
    )
    adata.obsm["X_pca_harmony"] = _to_float32(harmony_embed.Z_corr.T)
    del harmony_embed
    gc.collect()

    return adata


def process_query_data():
    """Convenience wrapper that executes the full query preprocessing flow."""
    adata = load_data()
    adata = check_data(adata, outdir="qc_plots")
    adata = preprocess_data(adata)
    return adata


def main():
    # Preprocess data
    adata = process_query_data()
    adata.write_h5ad("lungmap_query.h5ad", compression="gzip")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
