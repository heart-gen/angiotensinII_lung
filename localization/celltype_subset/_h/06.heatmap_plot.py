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

def load_markers(adata: sc.AnnData, marker_file: str, n_top: int = 5):
    """Load top N markers per cluster from TSV."""
    df = pd.read_csv(marker_file, sep="\t", index_col=0, nrows=100)
    clusters = sorted(set(col.split('_')[0] for col in df.columns))

    # Prepare list to accumulate melted DataFrames per cluster
    dfs = []
    for cluster in clusters:
        # Define column names for this cluster
        gene_col = f"{cluster}_names"
        pval_col = f"{cluster}_pvals"
        lfc_col = f"{cluster}_logfoldchanges"
        # Select and rename columns to unified names
        cluster_df = df[[gene_col, pval_col, lfc_col]].copy()
        cluster_df.columns = ['gene', 'pval', 'logfoldchange']
        # Add cluster identity column (as string)
        cluster_df['cluster'] = cluster
        # Append to list
        dfs.append(cluster_df)

    # Concatenate all clusters long DataFrames
    long_df = pd.concat(dfs, ignore_index=True)
    g_annot = pd.DataFrame(adata.var.loc[:, ["feature_name"]].reset_index())
    markers = g_annot.merge(long_df, left_on="ensembl_id", right_on="gene")

    # Filter out mitochrondria genes
    markers = markers[~markers["feature_name"].str.startswith("MT-")].copy()
    markers = markers[~markers["feature_name"].isin(["MALAT1"])].copy()
    top_markers = (
        markers.groupby("cluster")
        .apply(lambda x: x.nsmallest(n_top, "pval"))
        .reset_index(drop=True)
    )
    return top_markers


def filter_genes(adata: sc.AnnData, markers: pd.DataFrame) -> List[str]:
    """Extract marker genes present in AnnData var_names."""
    genes = markers["feature_name"].unique()
    new_genes = list(genes) + ["AGTR1", "ACTA2"]
    return [g for g in new_genes if g in adata.var.feature_name.values]


def plot_markers(adata: sc.AnnData, genes, groupby="leiden", prefix="markers", cmap="viridis"):
    """
    Generate both heatmap and dotplot for selected genes and save as PNG + PDF.
    """
    sc.pl.heatmap(adata, var_names=genes, groupby=groupby, cmap=cmap, dpi=300,
                  use_raw=False, layers="logcounts", gene_symbols="feature_name",
                  show=False, save=f"_{prefix}_heatmap.png")
    sc.pl.heatmap(adata, var_names=genes, groupby=groupby, cmap=cmap, use_raw=False,
                  layers="logcounts", gene_symbols="feature_name", show=False,
                  save=f"_{prefix}_heatmap.pdf")

    sc.pl.dotplot(adata, var_names=genes, groupby=groupby, use_raw=False, show=False,
                  layers="logcounts", dpi=300, dot_min=0.1, dot_max=1, color_map=cmap,
                  gene_symbols="feature_name", save=f"{prefix}_dotplot.png")
    sc.pl.dotplot(adata, var_names=genes, groupby=groupby, use_raw=False, show=False,
                  layers="logcounts", dot_min=0.1, dot_max=1, color_map=cmap,
                  gene_symbols="feature_name", save=f"{prefix}_dotplot.pdf")


def main(adata_file: str, marker_file: str, cluster_key: str = "leiden",
         n_top: int = 5, save_prefix: str = "markers"):
    """Main pipeline to run heatmap plotting."""
    # Load inputs
    adata = sc.read_h5ad(adata_file)
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X.copy()

    # Select top marker genes
    top_markers = load_markers(adata, marker_file, n_top=n_top)
    genes_to_plot = filter_genes(adata, top_markers)

    if not genes_to_plot:
        raise ValueError("No marker genes found in AnnData var_names!")

    # Plot results and save
    plot_markers(adata, genes_to_plot, prefix=save_prefix)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plot heatmap of top marker genes per cluster.")
    parser.add_argument("--adata", required=True, help="Path to AnnData (.h5ad) file")
    parser.add_argument("--markers", required=True, help="Path to TSV marker gene file")
    parser.add_argument("--cluster_key", default="leiden", help="Column in obs for cluster IDs")
    parser.add_argument("--n_top", type=int, default=5, help="Number of top markers per cluster")
    parser.add_argument("--save_prefix", default="marker_genes_heatmap",
                        help="Prefix for saved heatmap files (PDF and PNG)")
    args = parser.parse_args()

    main(args.adata, args.markers, cluster_key=args.cluster_key, n_top=args.n_top,
         save_prefix=args.save_prefix)
