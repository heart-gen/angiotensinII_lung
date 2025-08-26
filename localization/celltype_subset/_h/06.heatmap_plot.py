"""
Author: Kynon J Benjamin

Description:
    Loads marker genes and single-cell AnnData,
    selects top N marker genes per cluster,
    computes cluster-averaged expression,
    z-scores them, and plots a heatmap.
"""
import session_info
import scanpy as sc
import pandas as pd
from os import path
import seaborn as sns
import matplotlib.pyplot as plt

def load_markers(adata, marker_file, n_top=5):
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
    g_annot = pd.DataFrame(adata.var.loc[:, ["feature_name"]])
    markers = g_annot.merge(long_df, left_index=True, right_on="gene")

    # Filter out mitochrondria genes
    markers = markers[~markers["feature_name"].str.startswith("MT-")].copy()
    markers = markers[~markers["feature_name"].isin(["MALAT1"])].copy()
    top_markers = (
        markers.groupby("cluster")
        .apply(lambda x: x.nsmallest(n_top, "pval"))
        .reset_index(drop=True)
    )
    return top_markers


def filter_genes(adata, markers):
    """Extract marker genes present in AnnData var_names."""
    genes = markers["feature_name"].unique()
    new_genes = list(genes) + ["AGTR1", "ACTA2"]
    return [g for g in new_genes if g in adata.var.feature_name.values]


def plot_markers(adata, genes, groupby="leiden", prefix="markers", model="core",
                 cmap="viridis", outdir="figures", formats=('png', 'pdf')):
    """
    Generate both heatmap and dotplot for selected genes and save as PNG + PDF.
    """
    def save_plot(plot_func, fname, **kwargs):
        plt.figure()
        plot_func(show=False, **kwargs)
        for ext in formats:
            plt.savefig(path.join(outdir, model, f"{fname}.{ext}"), dpi=300,
                        bbox_inches='tight')
        plt.close()

    save_plot(lambda **kwargs: sc.pl.heatmap(
        adata, var_names=genes, groupby=groupby, cmap=cmap, #title="Subclusters",
        use_raw=False, layer="logcounts", gene_symbols="feature_name", **kwargs),
              f'{prefix}.heatmap_subclusters')

    save_plot(lambda **kwargs: sc.pl.dotplot(
        adata, var_names=genes, groupby=groupby, dot_min=0.1, dot_max=1,
        color_map=cmap, use_raw=False, layer="logcounts",
        gene_symbols="feature_name", **kwargs), f'{prefix}.dotplot_subclusters')


def main(adata_file, marker_file, cluster_key="leiden", n_top=5, model="core"):
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
    plot_markers(adata, genes_to_plot, model=model)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot heatmap of top marker genes per cluster.")
    parser.add_argument("--cluster_key", default="leiden",
                        help="Column in obs for cluster IDs")
    parser.add_argument("--n_top", type=int, default=5,
                        help="Number of top markers per cluster")
    parser.add_argument("--model", default="core",
                        help="Model type: 'core' or 'full'. Default: core")
    args = parser.parse_args()

    # Main analysis
    adata_file = f"pericyte.hlca_{args.model}.subclustered.h5ad"
    marker_file = f"figures/{args.model}/rank_genes_groups_results.txt.gz"
    main(adata_file, marker_file, cluster_key=args.cluster_key,
         n_top=args.n_top, model=args.model)

    # Session information
    session_info.show()
