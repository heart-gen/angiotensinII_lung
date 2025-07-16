## This script examines differences in expression
import session_info
import scanpy as sc
import pandas as pd
import seaborn as sns
import scikit_posthocs as sp
from os import makedirs, path
from scipy.stats import kruskal
import matplotlib.pyplot as plt

def extract_expression_by_cluster(adata, gene, cluster_key='leiden'):
    ensembl_id = adata.var[adata.var.feature_name == gene].index.values[0]
    df = adata[:, ensembl_id].to_df()
    df['cluster'] = adata.obs[cluster_key].values
    df['donor_id'] = adata.obs['donor_id'].values
    df.columns = ["expression", "cluster", "donor_id"]
    return df


def run_kruskal_test(df):
    grouped = [grp["expression"].values for _, grp in df.groupby("cluster")]
    stat, p = kruskal(*grouped)
    return stat, p


def run_posthoc_dunn(df):
    return sp.posthoc_dunn(df, val_col='expression', group_col='cluster',
                           p_adjust='fdr_bh')


def plot_boxplot_with_jitter(df, gene, stat, p_val, outdir):
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df, x="cluster", y="expression", whis=1.5, showfliers=False)
    sns.stripplot(data=df, x="cluster", y="expression", hue="donor_id",
                  dodge=True, jitter=True, linewidth=0.5, alpha=0.7)
    plt.title(f"{gene} expression by cluster\nKruskal-Wallis H={stat:.2f}, p={p_val:.1e}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Donor ID')
    plt.tight_layout()
    for ext in ['png', 'pdf']:
        plt.savefig(path.join(outdir, f"{gene.lower()}_expression_boxplot.{ext}"))
    plt.close()


def calculate_and_plot_stats(adata, marker_genes, outdir=".", cluster_key='leiden'):
    results = []
    for gene in marker_genes:
        ensembl_id = adata.var[adata.var.feature_name == gene].index.values[0]
        if ensembl_id not in adata.var_names:
            print(f"Warning: Gene {gene} not found in adata.var_names")
            continue

        df = extract_expression_by_cluster(adata, gene, cluster_key)
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
