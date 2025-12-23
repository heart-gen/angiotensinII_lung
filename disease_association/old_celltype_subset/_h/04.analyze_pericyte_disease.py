#!/usr/bin/env python3
"""
Analyze pericyte subclusters and disease association.
"""
import numpy as np
import session_info
import scanpy as sc
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import scikit_posthocs as sp
from os import makedirs, path
from scipy.stats import kruskal
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

def extract_agtr1_expression(adata, gene_name="AGTR1",
                             cluster_key="leiden", disease_key="Disease_Identity"):
    """Extract AGTR1 expression with cluster and disease annotations."""
    # Extract expression from logcounts layer
    if "logcounts" not in adata.layers:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        adata.layers["logcounts"] = adata.X

    expr = pd.Series(
        adata[:, gene_name].layers["logcounts"].toarray().flatten(),
        index=adata.obs_names, name="expression"
    )

    df = pd.DataFrame({
        "expression": expr,
        "cluster": adata.obs[cluster_key].astype(str).values,
        "disease": adata.obs[disease_key].astype(str).values,
        "donor_id": adata.obs.get("patient", "unknown")
    }).dropna(subset=["expression", "cluster", "disease"])

    # Aggregate per donor × cluster × disease to reduce pseudo-replication
    group_df = (
        df.groupby(["donor_id", "cluster", "disease"], observed=False)["expression"]
          .mean().reset_index()
    )

    # Remove NaNs
    group_df = group_df.dropna(subset=["expression", "cluster", "disease"]).copy()
    group_df = group_df[~(group_df["cluster"] == "4")].copy()
    return group_df


def run_agtr1_stats(df):
    """Run Kruskal–Wallis and ANOVA-like OLS model for AGTR1 expression."""
    # Kruskal–Wallis (cluster)
    grouped = [grp["expression"].values for _, grp in df.groupby("cluster")]
    stat_c, p_c = kruskal(*grouped)

    # Kruskal–Wallis (disease)
    grouped_d = [grp["expression"].values for _, grp in df.groupby("disease")]
    stat_d, p_d = kruskal(*grouped_d)

    # OLS model (ANOVA-style)
    model = ols("expression ~ C(cluster) * C(disease)", data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    return (stat_c, p_c), (stat_d, p_d), anova_table


def plot_agtr1_expression(df, gene, outdir):
    """Plot AGTR1 expression grouped by cluster × disease."""
    plt.figure(figsize=(9, 6))
    sns.boxplot(data=df, x="cluster", y="expression", hue="disease",
                palette="Set2", showfliers=False)
    sns.stripplot(data=df, x="cluster", y="expression", hue="disease",
                  dodge=True, jitter=True, color="black", alpha=0.4)
    plt.title(f"{gene} expression by pericyte subcluster and disease", fontsize=15)
    plt.xlabel("Pericyte Subcluster", fontsize=13)
    plt.ylabel("Normalized Expression", fontsize=13)
    plt.legend(title="Disease", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    for ext in ["png", "pdf"]:
        plt.savefig(path.join(outdir, f"{gene.lower()}_expression_by_disease.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()


def analyze_agtr1_expression(adata, gene_name="AGTR1",
                             cluster_key="leiden", outdir="."):
    """Wrapper to extract, analyze, and plot AGTR1 expression statistics."""
    df = extract_agtr1_expression(adata, gene_name, cluster_key)
    df.to_csv(path.join(outdir, f"{gene_name.lower()}_aggregate_data.tsv"),
              sep="\t", index=False)

    (stat_c, p_c), (stat_d, p_d), anova_table = run_agtr1_stats(df)
    print(f"\nKruskal–Wallis (cluster): H={stat_c:.2f}, p={p_c:.2e}")
    print(f"Kruskal–Wallis (disease): H={stat_d:.2f}, p={p_d:.2e}")
    print("\nANOVA-like model:\n", anova_table)

    anova_table.to_csv(path.join(outdir, f"{gene_name.lower()}_anova.tsv"), sep="\t")
    plot_agtr1_expression(df, gene_name, outdir)

    # Posthoc Dunn tests
    posthoc_cluster = sp.posthoc_dunn(df, val_col='expression', group_col='cluster')
    posthoc_cluster.to_csv(path.join(outdir, f"{gene_name.lower()}_dunn_cluster.tsv"), sep="\t")

    posthoc_disease = sp.posthoc_dunn(df, val_col='expression', group_col='disease')
    posthoc_disease.to_csv(path.join(outdir, f"{gene_name.lower()}_dunn_disease.tsv"), sep="\t")


def calculate_cluster_proportions(adata, cluster_key="leiden",
                                  disease_key="Disease_Identity"):
    """Compute per-donor per-cluster proportions across disease states."""
    df = adata.obs[[cluster_key, disease_key, "patient"]].copy()
    df["count"] = 1
    df = df.dropna(subset=[cluster_key, disease_key, "patient"]).copy()

    donor_cluster = (
        df.groupby(["patient", disease_key, cluster_key], observed=False)["count"]
          .sum().reset_index()
    )

    donor_totals = donor_cluster.groupby(["patient", disease_key], observed=False)["count"].transform("sum")
    donor_cluster["proportion"] = donor_cluster["count"] / donor_totals

    donor_cluster.rename(columns={cluster_key: "cluster"}, inplace=True)
    donor_cluster = donor_cluster.dropna(subset=["proportion"]).copy()
    donor_cluster = donor_cluster[~(donor_cluster["cluster"] == "4")].copy()
    donor_cluster["proportion_transformed"] = np.arcsin(np.sqrt(donor_cluster["proportion"]))
    return donor_cluster


def test_cluster_proportion_differences(prop_df, outdir="."):
    """ANOVA proportions across diseases."""
    model = ols("proportion_transformed ~ C(cluster) * C(disease)", data=prop_df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_table.to_csv(path.join(outdir, "cluster_proportion_anova.tsv"), sep="\t")
    return anova_table


def plot_cluster_proportions(prop_df, outdir="."):
    """Visualize pericyte subcluster proportions by disease."""
    plt.figure(figsize=(9, 6))
    sns.boxplot(data=prop_df, x="cluster", y="proportion_transformed", hue="disease",
                palette="Set2", showfliers=False)
    sns.stripplot(data=prop_df, x="cluster", y="proportion", hue="disease",
                  dodge=True, jitter=True, color="black", alpha=0.5)
    plt.title("Pericyte subcluster proportions by disease", fontsize=15)
    plt.xlabel("Pericyte Subcluster", fontsize=13)
    plt.ylabel("Transformed Proportion (per donor)", fontsize=13)
    plt.legend(title="Disease", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(path.join(outdir, "cluster_proportions_by_disease.png"),
                dpi=300, bbox_inches="tight")
    plt.savefig(path.join(outdir, "cluster_proportions_by_disease.pdf"),
                bbox_inches="tight")
    plt.close()


def main():
    outdir = "pericyte_disease_analysis"
    makedirs(outdir, exist_ok=True)

    print("\n[INFO] Loading data:")
    adata = sc.read_h5ad("clustered_data.h5ad")

    print("\n[STEP 1] Analyzing AGTR1 expression...")
    analyze_agtr1_expression(adata, gene_name="AGTR1",
                             cluster_key="predicted_labels", outdir=outdir)

    print("\n[STEP 2] Analyzing pericyte subcluster proportions...")
    prop_df = calculate_cluster_proportions(adata, cluster_key="predicted_labels",
                                            disease_key="disease")
    test_cluster_proportion_differences(prop_df, outdir=outdir)
    plot_cluster_proportions(prop_df, outdir=outdir)

    print("\n[INFO] Analysis complete. Results saved to:", outdir)
    session_info.show()


if __name__ == "__main__":
    main()
