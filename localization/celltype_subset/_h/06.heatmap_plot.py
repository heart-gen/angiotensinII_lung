"""
Author: Kynon J Benjamin

Description:
    Loads marker genes and single-cell AnnData,
    selects top N marker genes per cluster,
    computes cluster-averaged expression,
    z-scores them, and plots a heatmap.
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
from typing import List
import matplotlib.pyplot as plt

def load_markers(marker_file: str, n_top: int = 5):
    """Load top N markers per cluster from CSV."""
    markers = pd.read_csv(marker_file)
    top_markers = (
        markers.groupby("cellType.target")
        .apply(lambda x: x.nsmallest(n_top, "MeanRatio.rank"))
        .reset_index(drop=True)
    )
    return top_markers


def filter_genes(adata: sc.AnnData, markers: pd.DataFrame) -> List[str]:
    """Extract marker genes present in AnnData var_names."""
    genes = markers["gene"].unique()
    return [g for g in genes if g in adata.var_names]


def compute_cluster_averages(adata: sc.AnnData, genes: List[str],
                             cluster_key: str = "cluster"):
    """
    Compute cluster-averaged expression of selected genes.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    genes : list
        Gene names to compute averages for.
    cluster_key : str
        Key in adata.obs that contains cluster IDs.
    """
    df = pd.DataFrame(
        adata[:, genes].X.toarray(),
        columns=genes,
        index=adata.obs[cluster_key]
    )
    cluster_avg = df.groupby(level=0).mean()
    return cluster_avg


def scale_expression(cluster_avg: pd.DataFrame):
    """Z-score scale per gene across clusters."""
    return (cluster_avg - cluster_avg.mean()) / cluster_avg.std()


def plot_heatmap(scaled: pd.DataFrame, title: str = "Top Marker Genes per Subcluster",
                 save_pdf: str = None, save_png: str = None):
    """Plot and save heatmap of scaled cluster-averaged gene expression."""
    plt.figure(figsize=(10, 8))
    sns.heatmap(scaled, cmap="vlag", center=0, cbar_kws={'label': 'z-score'})
    plt.title(title)
    plt.xlabel("Genes")
    plt.ylabel("Clusters")
    plt.tight_layout()
    if save_pdf:
        plt.savefig(save_pdf, format='pdf')
    if save_png:
        plt.savefig(save_png, format='png', dpi=300)
    plt.show()


def main(adata_file: str, marker_file: str, cluster_key: str = "cellType.target",
         n_top: int = 5, save_prefix: str = "marker_genes_heatmap"):
    """Main pipeline to run heatmap plotting."""
    # Load inputs
    adata = sc.read_h5ad(adata_file)
    top_markers = load_markers(marker_file, n_top=n_top)
    genes_to_plot = filter_genes(adata, top_markers)

    if not genes_to_plot:
        raise ValueError("No marker genes found in AnnData var_names!")

    # Compute average expression per cluster
    cluster_avg = compute_cluster_averages(adata, genes_to_plot,
                                           cluster_key=cluster_key)
    scaled = scale_expression(cluster_avg)

    # Define file names
    pdf_file = f"{save_prefix}.pdf"
    png_file = f"{save_prefix}.png"

    # Plot results and save
    plot_heatmap(
        scaled,
        title=f"Top {n_top} Marker Genes per Subcluster",
        save_pdf=pdf_file,
        save_png=png_file,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot heatmap of top marker genes per cluster.")
    parser.add_argument("--adata", required=True, help="Path to AnnData (.h5ad) file")
    parser.add_argument("--markers", required=True, help="Path to CSV marker gene file")
    parser.add_argument("--cluster_key", default="cellType.target", help="Column in obs for cluster IDs")
    parser.add_argument("--n_top", type=int, default=5, help="Number of top markers per cluster")
    parser.add_argument("--save_prefix", default="marker_genes_heatmap",
                        help="Prefix for saved heatmap files (PDF and PNG)")
    args = parser.parse_args()

    main(args.adata, args.markers, cluster_key=args.cluster_key, n_top=args.n_top,
         save_prefix=args.save_prefix)
