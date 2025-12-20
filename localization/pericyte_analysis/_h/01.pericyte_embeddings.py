"""
Compute embeddings + marker overlays for pericytes.
"""
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from scipy import sparse
from anndata import AnnData
import matplotlib.pyplot as plt
from typing import Dict, List, Optional

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
    parser.add_argument("--cluster-key", default="leiden_pericytes")
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
                    metric="cosine", random_state=seed)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=seed)
    sc.tl.tsne(adata, use_rep=use_rep, perplexity=perplexity, random_state=seed)


def run_leiden(adata: AnnData, resolution: float = 0.5, seed: int = 13):
    logging.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(
        adata, resolution=resolution, key_added="leiden_pericytes", random_state=seed,
    )


def add_gene_expression(adata: AnnData, gene="AGTR1"):
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


def analyze_marker_genes(adata, groupby: str = 'leiden_pericytes',
                         method: str = 'wilcoxon', outdir: Path = Path("./")):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method,
                            layer="logcounts", use_raw=False)
    # Reformat results
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    rank_df = pd.DataFrame({
        group + '_' + key: result[key][group]
        for group in groups
        for key in ['names', 'pvals', 'logfoldchanges']
    })
    rank_df.to_csv(outdir / "rank_genes_groups_results.txt.gz", sep='\t')
    return adata


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


def plot_ranked_marker_dotplot(
        adata: AnnData, cluster_key: str, outdir: Path, n_genes: int = 10,
):
    """Generate a dotplot from marker genes of pericytes subclusters"""
    results = adata.uns["rank_genes_groups"]
    groups = results["names"].dtype.anmes
    top_genes = set()

    for group in groups:
        names = results["names"][group][:n_genes]
        top_genes.update(names)

    top_genes = list(top_genes)
    
    dp = sc.pl.dotplot(
        adata, top_genes, groupby=cluster_key,
        show=False, return_fig=False, dendrogram=True,
    )
    fig = next(iter(dp.values())).figure
    ax = dp.get("gene_group_ax", list(dp.values())[-1])
    for t in ax.get_xticklabels():
        t.set_ha("center")
        t.set_rotation(90)
        t.set_rotation_mode("anchor")
        
    fig.suptitle("To marker genes", y=1.02)
    save_figure(fig, outdir / "dotplot_top_markers")


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

    # Gene UMAP plots
    genes_dir = outdir / "gene_plots"
    genes_dir.mkdir(exist_ok=True)

    for gene in ["AGTR1", "ACTA2"]:
        add_gene_expression(adata, gene=gene)
        df = make_df(adata, "X_umap", [f"{gene}_expr", f"{gene}_detect", "leiden"])
        plot_seaborn_embedding(df, f"{gene}_expr", f"{gene} expression",
                               genes_dir / f"umap_{gene.lower()}_expr", categorical=False)
        plot_seaborn_embedding(df, f"{gene}_detect", f"{gene} detection",
                               genes_dir / f"umap_{gene.lower()}_detect", categorical=True)

    plot_seaborn_embedding(df, args.cluster_key, f"Pericytes Subclusters",
                           outdir / f"umap_leiden", categorical=True)

    # Marker genes
    marker_dir = outdir / "markers"
    marker_dir.mkdir(exist_ok=True)

    adata = analyze_marker_genes(
        adata, args.cluster_key, 'wilcoxon', marker_dir
    )
    plot_ranked_marker_dotplot(adata, args.cluster_key, marker_dir, 2)

    # Save updated AnnData
    df = adata.obs.copy()
    df.to_csv(outdir / "pericytes_metadata.tsv.gz", sep="\t")
    
    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(
                adata.var.index.to_numpy(),
                adata.var[col].to_numpy()
        ):
            adata.var.index.name = None

    adata.write(outdir / "pericyte_with_embeddings.h5ad")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
