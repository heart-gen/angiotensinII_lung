import session_info
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_data(h5ad_file, cluster_key="cluster"):
    """
    Load AnnData object and extract clustering info.
    """
    adata = sc.read_h5ad(h5ad_file)
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs")
    return adata, adata.obs[cluster_key]


def compute_cluster_avg(adata, clusters):
    """
    Compute mean expression per gene per cluster (using sparse matrix).
    Returns DataFrame (genes x clusters).
    """
    cluster_avg_expr = {}
    for cluster in clusters.unique():
        cells = clusters[clusters == cluster].index
        sub = adata[cells]
        cluster_avg_expr[cluster] = np.asarray(sub.X.mean(axis=0)).ravel()
    cluster_avg_expr_df = pd.DataFrame(cluster_avg_expr, index=adata.var_names)
    return cluster_avg_expr_df


def calc_specificity_ratio(cluster_avg_expr_df):
    """
    Specificity ratio = expression in cluster / sum(expression in other clusters).
    """
    specificity = pd.DataFrame(index=cluster_avg_expr_df.index)
    for cluster in cluster_avg_expr_df.columns:
        cluster_expr = cluster_avg_expr_df[cluster]
        other_expr = cluster_avg_expr_df.drop(columns=cluster).sum(axis=1)
        specificity[cluster] = cluster_expr / (other_expr + 1e-9)
    return specificity


def run_de_analysis(adata, cluster_key="cluster", method="wilcoxon",
                    layer="logcounts"):
    """
    Run Scanpy's DE test for all clusters vs. rest.
    Returns: dict of DataFrames with DE results.
    """
    if 'rank_genes_groups' not in adata.uns:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key,
                                method=method, layer=layer)
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names

    # Convert to tidy DataFrames per cluster
    de_results = {}
    for g in groups:
        df = pd.DataFrame({
            "gene": result["names"][g],
            "logfc": result["logfoldchanges"][g],
            "pval": result["pvals"][g],
            "pval_adj": result["pvals_adj"][g],
            "score": result["scores"][g],
        })
        de_results[g] = df
    return de_results


def select_markers(cluster_avg_expr_df, specificity_ratio, de_results,
                   specificity_threshold=4.0, expression_threshold=1.0,
                   fdr_threshold=0.05):
    """
    Select markers using specificity, expression, and adjusted p-values.
    Returns tidy DataFrame of candidate markers.
    """
    markers = []
    for cluster, de_df in de_results.items():
        de_df = de_df.set_index("gene")

        # Ensure consistent indexing across DataFrames
        common_genes = (
            de_df.index
            .intersection(specificity_ratio.index)
            .intersection(cluster_avg_expr_df.index)
        )

        # Subset everything to common genes
        de_df = de_df.loc[common_genes]
        spec = specificity_ratio.loc[common_genes, cluster]
        avg = cluster_avg_expr_df.loc[common_genes, cluster]

        # Candidate selection
        candidates = common_genes[
            (spec > specificity_threshold) &
            (avg > expression_threshold) &
            (de_df["pval_adj"] < fdr_threshold)
        ]

        for gene in candidates:
            markers.append({
                "cluster": cluster,
                "gene": gene,
                "specificity_ratio": specificity_ratio.at[gene, cluster],
                "mean_expression": cluster_avg_expr_df.at[gene, cluster],
                "pval_adj": de_df.at[gene, "pval_adj"],
                "logfc": de_df.at[gene, "logfc"]
            })

    markers_df = pd.DataFrame(markers).sort_values(
        ["cluster", "specificity_ratio"], ascending=[True, False]
    )
    return markers_df


def plot_marker_expression(adata, gene, cluster_key="cluster",
                           layer="logcounts", outdir="plots"):
    """
    Violin plot of gene expression across clusters from a specified layer.
    Saves PNG and PDF versions of the plot.
    """
    from os import makedirs, path
    makedirs(outdir, exist_ok=True)

    # Get expression values from chosen layer
    X = adata[:, gene].layers[layer] if layer in adata.layers else adata[:, gene].X
    values = X.toarray().ravel() if hasattr(X, "toarray") else X.ravel()

    # Build DataFrame for plotting
    df_plot = pd.DataFrame({
        "expression": values,
        "cluster": adata.obs[cluster_key]
    })

    # Plot
    plt.figure(figsize=(8, 4))
    sns.violinplot(x="cluster", y="expression", data=df_plot, inner="quartile")
    plt.title(f"Expression of {gene} across clusters ({layer})")

    # Save files
    png_path = path.join(outdir, f"{gene.lower()}_expression.png")
    pdf_path = path.join(outdir, f"{gene.lower()}_expression.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()


def main():
    # Set parameters
    cluster_key = "cell_types"
    h5ad_file = "../../_m/clustered_data.with_celltypes.h5ad"
    # Load data
    adata, clusters = load_data(h5ad_file, cluster_key)
    # Calculate specificity scores
    cluster_avg_expr_df = compute_cluster_avg(adata, clusters)
    specificity_ratio = calc_specificity_ratio(cluster_avg_expr_df)
    de_results = run_de_analysis(adata, cluster_key)
    # Identify marker genes
    markers_df = select_markers(cluster_avg_expr_df, specificity_ratio,
                                de_results, specificity_threshold=1.25)

    # Save data
    specificity_ratio.to_csv("specificity_data.tsv", sep="\t")
    markers_df.to_csv("candidate_markers_for_qpcr.tsv", sep="\t", index=False)
    print(f"Saved {len(markers_df)} candidate markers.")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
