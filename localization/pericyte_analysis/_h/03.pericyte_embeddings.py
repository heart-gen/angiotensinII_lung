"""
Compute embeddings + marker overlays for pericytes.
"""
import logging, argparse
from pathlib import Path
from typing import Dict, List, Optional

import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from scipy import sparse
from anndata import AnnData
import matplotlib.pyplot as plt

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
                        help="Subsetted pericyte AnnData")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--markers-yaml", required=True, type=Path)
    parser.add_argument("--cluster-key", default="leiden")
    parser.add_argument("--use-rep", default="X_pca_harmony")
    parser.add_argument("--neighbors", type=int, default=50)
    parser.add_argument("--umap-min-dist", type=float, default=0.9)
    parser.add_argument("--umap-spread", type=float, default=1.5)
    parser.add_argument("--tsne-perplexity", type=float, default=30.0)
    parser.add_argument("--leiden-resolution", type=float, default=0.5,
                        help="Leiden clustering resolution")
    parser.add_argument("--seed", type=int, default=13)
    return parser.parse_args()


def save_figure(fig, base_path: Path):
    fig.savefig(base_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_anndata(path: Path) -> AnnData:
    adata = sc.read_h5ad(path)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    if "feature_name" not in adata.var.columns:
        logging.warning("No feature_name column detected. Using var_names as symbols.")
        adata.var["feature_name"] = adata.var_names
    symbols = adata.var["feature_name"].astype(str)
    adata.var["ensembl_id"] = adata.var_names
    adata.var_names = symbols
    adata.var.index = symbols
    return adata


def compute_embeddings(adata: AnnData, use_rep: str,
                       neighbors=30, min_dist=0.5, spread=1.2,
                       perplexity=30, seed=13):
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=neighbors,
                    random_state=seed) # metric="cosine"
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=seed)
    sc.tl.tsne(adata, use_rep=use_rep, perplexity=perplexity, random_state=seed)


def run_leiden(adata: AnnData, resolution: float = 0.5, seed: int = 13):
    logging.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(
        adata, resolution=resolution, key_added="leiden", random_state=seed,
    )


def add_agtr1_expression(adata: AnnData, gene="AGTR1"):
    if gene not in adata.var_names:
        logging.warning(f"{gene} not in var_names")
        return

    expr = adata[:, gene].layers["logcounts"]
    expr = expr.toarray().ravel() if sparse.issparse(expr) else np.asarray(expr).ravel()

    adata.obs[f"{gene}_expr"] = expr
    adata.obs[f"{gene}_detect"] = (expr > 0).astype(int)


def make_df(adata: AnnData, emb_key: str, cols: List[str]) -> pd.DataFrame:
    coords = adata.obsm[emb_key]
    df = pd.DataFrame(coords, columns=["dim1", "dim2"], index=adata.obs_names)
    for c in cols:
        df[c] = adata.obs[c]
    return df


def plot_seaborn_embedding(df: pd.DataFrame, hue: str,
                           title: str, base: Path, categorical=None):
    values = df[hue]
    if categorical is None:
        categorical = not pd.api.types.is_numeric_dtype(values)

    fig, ax = plt.subplots(figsize=(6, 5))
    if categorical:
        cat = pd.Categorical(values)
        palette = sns.color_palette("tab20", len(cat.categories))
        for col, cat_val in zip(palette, cat.categories):
            mask = cat == cat_val
            ax.scatter(df.loc[mask, "dim1"],
                       df.loc[mask, "dim2"],
                       s=6, linewidths=0, color=col, label=str(cat_val))
        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize="small")
    else:
        sc = ax.scatter(df["dim1"], df["dim2"], c=values, s=6, linewidths=0,
                        cmap="viridis")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    ax.set_title(title)
    save_figure(fig, base)


def load_marker_panels(path: Path) -> Dict[str, List[str]]:
    data = yaml.safe_load(path.read_text())
    if "pericyte" not in data:
        raise ValueError("YAML must contain top-level key 'pericyte'")
    return data["pericyte"]


def plot_marker_panels_umap(adata: AnnData, panels: Dict[str, List[str]],
                            outdir: Path):
    """UMAP feature plots with gene names."""
    for label, mapping in panels.items():
        symbols = panels[label]
        fig = sc.pl.umap(adata, color=symbols, layer="logcounts",
                         show=False, return_fig=True)

        save_figure(fig, outdir / f"umap_markers_{label}")


def extract_dotplot_figure(dp):
    """
    Normalize Scanpy dotplot output to a matplotlib Figure.
    """
    if hasattr(dp, "make_figure"): # DotPlot object
        return dp.make_figure()
    if hasattr(dp, "fig"):
        return dp.fig
    if isinstance(dp, dict) and "fig" in dp: # dict output
        return dp["fig"]
    if isinstance(dp, tuple) and len(dp) == 2: # tuple output: (fig, ax)
        fig, _ = dp
        return fig
    raise RuntimeError(f"Cannot extract figure from object of type {type(dp)}")


def plot_marker_dotplot(adata: AnnData, panels: Dict[str, Dict[str, str]],
                        cluster_key: str, outdir: Path):
    for panel_name, symbols in panels.items():
        dp = sc.pl.dotplot(
            adata, symbols, groupby=cluster_key,
            show=False, return_fig=False, dendrogram=True,
        )
        fig = next(iter(dp.values())).figure
        ax = dp.get("gene_group_ax", list(dp.values())[-1])
        for t in ax.get_xticklabels():
            t.set_ha("center")
            t.set_rotation(90)
            t.set_rotation_mode("anchor")
        fig.suptitle(f"{panel_name} markers", y=1.02)
        save_figure(fig, outdir / f"dotplot_{panel_name}")


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)
    adata.raw = None
    adata.var_names_make_unique()

    # Embedding
    compute_embeddings(
        adata, args.use_rep,
        neighbors=args.neighbors,
        min_dist=args.umap_min_dist,
        spread=args.umap_spread,
        perplexity=args.tsne_perplexity,
        seed=args.seed
    )

    # Leiden clustering
    run_leiden(adata, resolution=args.leiden_resolution, seed=args.seed)

    # AGTR1
    add_agtr1_expression(adata, gene="AGTR1")

    # AGTR1 UMAP / t-SNE
    agtr_dir = outdir / "agtr1"
    agtr_dir.mkdir(exist_ok=True)

    for emb in ["X_umap", "X_tsne"]:
        df = make_df(adata, emb, ["AGTR1_expr", "AGTR1_detect", "leiden"])
        plot_seaborn_embedding(df, "AGTR1_expr", f"{emb}: AGTR1 expression",
                               agtr_dir / f"{emb}_agtr1_expr", categorical=False)
        plot_seaborn_embedding(df, "AGTR1_detect", f"{emb}: AGTR1 detection",
                               agtr_dir / f"{emb}_agtr1_detect", categorical=True)
        plot_seaborn_embedding(df, "leiden", f"{emb}: Leiden clusters",
                               outdir / f"{emb}_leiden", categorical=True)

    # Marker panels
    panels = load_marker_panels(args.markers_yaml)
    for p, genes in panels.items():
        panels[p] = [g for g in genes if g in adata.var_names]

    marker_dir = outdir / "markers"
    marker_dir.mkdir(exist_ok=True)

    plot_marker_panels_umap(adata, panels, marker_dir)
    plot_marker_dotplot(adata, panels, args.cluster_key, marker_dir)

    # Save updated AnnData
    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(
                adata.var.index.to_numpy(),
                adata.var[col].to_numpy()
        ):
            adata.var.index.name = None

    adata.write(outdir / "pericyte_with_embeddings.h5ad")

    session_info.show()


if __name__ == "__main__":
    main()
