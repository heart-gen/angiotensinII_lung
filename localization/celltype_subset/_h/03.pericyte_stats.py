## This script examines differences in expression
import session_info
import scanpy as sc
import pandas as pd
import seaborn as sns
import scikit_posthocs as sp
from os import makedirs, path
from scipy.stats import kruskal
import matplotlib.pyplot as plt

def extract_expression(adata, marker_genes, cluster_key="leiden"):
    # Get corresponding ensembl ids
    ensembl_ids = adata.var[adata.var.feature_name.isin(marker_genes)].index.tolist()
    feature_map = adata.var.loc[ensembl_ids, "feature_name"].to_dict()

    # Subset expression from logcounts layer
    expr = pd.DataFrame(adata[:, ensembl_ids].layers["logcounts"].toarray(),
                        columns=ensembl_ids, index=adata.obs_names)

    # Annotate with metadata
    expr['cluster'] = adata.obs[cluster_key].values
    expr['donor_id'] = adata.obs['donor_id'].values

    # Melt for long format
    df = expr.melt(id_vars=["cluster", "donor_id"],
                   value_vars=ensembl_ids,
                   var_name="gene_id", value_name="expression")\
             .dropna(axis=0)

    # Map gene id to feature name
    df["gene_name"] = df["gene_id"].map(feature_map)
    group_df = df.groupby(["cluster", "donor_id", "gene_name"])["expression"]\
                 .mean().reset_index().dropna(subset=["expression"])
    return group_df


def run_kruskal_test(df):
    grouped = [grp["expression"].values for _, grp in df.groupby("cluster")]
    stat, p = kruskal(*grouped)
    return stat, p


def run_posthoc_dunn(df):
    return sp.posthoc_dunn(df, val_col='expression', group_col='cluster',
                           p_adjust='fdr_bh')


def plot_boxplot_with_jitter(df, gene, stat, p_val, outdir):
    plt.figure(figsize=(7, 6))
    sns.boxplot(data=df, x="cluster", y="expression", whis=1.5,
                showfliers=False, width=0.5, color="gray")
    sns.stripplot(data=df, x="cluster", y="expression", color="black",
                  dodge=True, jitter=True, linewidth=0.5, alpha=0.5)
    plt.title(f"{gene} expression by cluster\nKruskal-Wallis H={stat:.2f}, p={p_val:.1e}",
              fontsize=15)
    plt.xlabel("Pericyte Subclusters", fontsize=20)
    plt.ylabel("Normalized Expression", fontsize=20)
    plt.xticks(fontsize=15); plt.yticks(fontsize=15)
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(path.join(outdir, f"{gene.lower()}_expression_boxplot.{ext}"))
    plt.close()


def calculate_and_plot_stats(adata, marker_genes, outdir=".", cluster_key='leiden'):
    big_df = extract_expression(adata, marker_genes, cluster_key)
    big_df.to_csv(path.join(outdir, "aggregate_data.tsv"), index=False, sep="\t")
    results = []
    for gene in marker_genes:
        df = big_df[(big_df["gene_name"] == gene)].copy()
        stat, p_val = run_kruskal_test(df)
        posthoc_df = run_posthoc_dunn(df)
        df_out = path.join(outdir, f"{gene.lower()}_dunn_posthoc.tsv")
        posthoc_df.to_csv(df_out, sep='\t')

        plot_boxplot_with_jitter(df, gene, stat, p_val, outdir)
        results.append({"gene": gene, "kruskal_H": stat, "p_value": p_val})

    # Export summary stats
    pd.DataFrame(results).to_csv(path.join(outdir, "kruskal_summary.tsv"),
                                 sep='\t', index=False)


def main():
    # General setup
    marker_genes = ['AGTR1', 'ACTA2']
    output_dir = "cluster_marker_stats"
    makedirs(output_dir, exist_ok=True)
    # Load data
    adata = sc.read_h5ad("pericyte.hlca_core.subclustered.h5ad")
    # Calculate and plot statistics
    calculate_and_plot_stats(adata, marker_genes, outdir=output_dir,
                             cluster_key='leiden')
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
