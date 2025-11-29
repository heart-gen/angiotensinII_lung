"""
Compute donor mixing metrics (iLISI) + mixing plots for pericytes.
"""
import json
import logging
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from pathlib import Path
from anndata import AnnData
import matplotlib.pyplot as plt
from scib.metrics import ilisi_graph

sns.set_context("talk")
sns.set_style("whitegrid")

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata", required=True, type=Path,
                        help="pericyte_with_embeddings.h5ad")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--use-rep", default="X_pca_harmony")
    parser.add_argument("--seed", type=int, default=13)
    return parser.parse_args()


def save_fig(fig, base: Path):
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_adata(path: Path) -> AnnData:
    return sc.read_h5ad(path)


def make_df(adata: AnnData, emb_key: str, donor_key: str):
    coords = adata.obsm[emb_key]
    df = pd.DataFrame(coords, columns=["dim1", "dim2"], index=adata.obs_names)
    df[donor_key] = adata.obs[donor_key]
    return df


def plot_donor_embedding(df: pd.DataFrame, donor_key: str,
                         title: str, base: Path):
    cat = pd.Categorical(df[donor_key])
    palette = sns.color_palette("tab20", len(cat.categories))

    fig, ax = plt.subplots(figsize=(6, 5))
    for col, cat_val in zip(palette, cat.categories):
        mask = cat == cat_val
        ax.scatter(df.loc[mask, "dim1"], df.loc[mask, "dim2"],
                   s=6, linewidths=0, color=col, label=str(cat_val))
    ax.legend().remove(); ax.set_title(title)
    save_fig(fig, base)


def compute_ilisi(adata: AnnData, donor_key: str, use_rep: str) -> float:
    mask = adata.obs[donor_key].notna()
    ad = adata[mask].copy()
    ad.obs[donor_key] = ad.obs[donor_key].astype(str)
    score = ilisi_graph(
        ad, batch_key=donor_key, use_rep=use_rep, type_="embed", n_cores=1
    )
    return float(score)


def write_json(path: Path, adata: AnnData, score: float,
               donor_key: str, use_rep: str):
    payload = {
        "embedding": {
            "use_rep": use_rep,
            "donor_key": donor_key,
            "n_cells": int(adata.n_obs),
            "n_donors": int(adata.obs[donor_key].nunique()),
            "iLISI_donor": None if not np.isfinite(score) else score,
        }
    }
    path.write_text(json.dumps(payload, indent=2))


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_adata(args.adata)

    # Make mixing plots
    mix_dir = outdir / "donor_mixing"
    mix_dir.mkdir(exist_ok=True)

    for emb in ["X_umap", "X_tsne"]:
        df = make_df(adata, emb, args.donor_key)
        plot_donor_embedding(df, args.donor_key,
                             f"{emb}: donor mixing",
                             mix_dir / f"{emb}_donor")

    # Compute iLISI
    score = compute_ilisi(
        adata, donor_key=args.donor_key, use_rep=args.use_rep,
    )

    # Write JSON
    write_json(outdir / "pericyte_embedding.json",
               adata, score, args.donor_key,
               args.use_rep, args.neighbors, args.seed)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
