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
    parser.add_argument("--neighbors", type=int, default=30)
    parser.add_argument("--umap-min-dist", type=float, default=0.5)
    parser.add_argument("--umap-spread", type=float, default=1.2)
    parser.add_argument("--tsne-perplexity", type=float, default=30.0)
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
    return adata


def compute_embeddings(adata: AnnData, use_rep: str,
                       neighbors=30, min_dist=0.5, spread=1.2,
                       perplexity=30, seed=13):
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=neighbors,
                    metric="cosine", random_state=seed)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=seed)
    sc.tl.tsne(adata, use_rep=use_rep, perplexity=perplexity, random_state=seed)


def symbol_to_varname(adata: AnnData, symbol: str) -> Optional[str]:
    """Return var_name (ENSEMBL) for a gene symbol."""
    if symbol in adata.var["feature_name"].values:
        return adata.var.index[adata.var["feature_name"] == symbol][0]
    elif symbol in adata.var_names:
        return symbol
    logging.warning("Gene symbol '%s' not found.", symbol)
    return None


def resolve_gene_symbol(adata: AnnData, gene: str) -> Optional[str]:
    if "feature_name" in adata.var.columns:
        mask = adata.var["feature_name"] == gene
        if mask.any():
            return adata.var.index[mask][0]
    if gene in adata.var_names:
        return gene
    logging.warning("Marker gene '%s' not found.", gene)
    return None


def add_agtr1_expression(adata: AnnData, gene="AGTR1"):
    var = resolve_gene_symbol(adata, gene)
    if var is None:
        return
    expr = adata[:, var].layers["logcounts"]
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


def convert_panels_to_varnames(adata: AnnData, panels) -> Dict[str, List[str]]:
    out = {}
    for label, symbols in panels.items():
        mapping = {}
        for sym in symbols:
            var = symbol_to_varname(adata, sym)
            if var is not None:
                mapping[sym] = var
        if mapping:
            out[label] = mapping
    return out


def plot_marker_dotplot(adata: AnnData, panels: Dict[str, List[str]],
                        cluster_key: str, outdir: Path):
    groups = list(panels.keys())

    all_symbols = []; all_varnames = []
    group_positions = []; current_idx = 0

    for group in groups:
        mapping = panels[group]
        syms = list(mapping.keys())
        vars_ = [mapping[s] for s in syms]

        start = current_idx
        end = current_idx + len(vars_) - 1
        group_positions.append((start, end))

        all_symbols.extend(syms)
        all_varnames.extend(vars_)
        current_idx += len(vars_)

    dp = sc.pl.dotplot(
        adata, var_names=all_varnames, groupby=cluster_key,
        var_group_labels=groups, var_group_positions=group_positions,
        standard_scale="var", show=False, return_fig=True,
    )

    dp.ax.set_xticklabels(all_symbols, rotation=90)
    save_figure(dp, outdir / "dotplot_markers")


def plot_marker_panels_umap(adata: AnnData, panels: Dict[str, List[str]],
                            outdir: Path):
    """UMAP feature plots with gene SYMBOLS."""
    for label, mapping in panels.items():
        symbols = list(mapping.keys())
        varnames = list(mapping.values())

        fig = sc.pl.umap(adata, color=varnames, layer="logcounts",
                         show=False, return_fig=True)

        for ax, sym in zip(fig.axes, symbols):            
            ax.set_title(sym)

        save_figure(fig, outdir / f"umap_markers_{label}")


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)

    # Embedding
    compute_embeddings(
        adata, args.use_rep,
        neighbors=args.neighbors,
        min_dist=args.umap_min_dist,
        spread=args.umap_spread,
        perplexity=args.tsne_perplexity,
        seed=args.seed
    )

    # AGTR1
    add_agtr1_expression(adata, gene="AGTR1")

    # AGTR1 UMAP / t-SNE
    agtr_dir = outdir / "agtr1"
    agtr_dir.mkdir(exist_ok=True)

    for emb in ["X_umap", "X_tsne"]:
        df = make_df(adata, emb, ["AGTR1_expr", "AGTR1_detect"])
        plot_seaborn_embedding(df, "AGTR1_expr",
                               f"{emb}: AGTR1 expression",
                               agtr_dir / f"{emb}_agtr1_expr", categorical=False)
        plot_seaborn_embedding(df, "AGTR1_detect",
                               f"{emb}: AGTR1 detection",
                               agtr_dir / f"{emb}_agtr1_detect", categorical=True)

    # Marker panels
    raw_panels = load_marker_panels(args.markers_yaml)
    panels = convert_panels_to_varnames(adata, raw_panels)

    marker_dir = outdir / "markers"
    marker_dir.mkdir(exist_ok=True)

    plot_marker_dotplot(adata, panels, args.cluster_key, marker_dir)
    plot_marker_panels_umap(adata, panels, marker_dir)

    # Save updated AnnData
    adata.write(outdir / "pericyte_with_embeddings.h5ad")

    session_info.show()


if __name__ == "__main__":
    main()
