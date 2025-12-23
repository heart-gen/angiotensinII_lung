"""
Compute donor-level age correlation with AGTR1 expression globally and per cell type.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from pyhere import here
from pathlib import Path
import logging, argparse
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from pandas.api.types import CategoricalDtype
from statsmodels.stats.multitest import fdrcorrection

sns.set_context("talk")
sns.set_style("whitegrid")

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--outdir", default=Path("./"), type=Path,
                        help="Output directory for results")
    parser.add_argument("--gene", default="AGTR1", type=str, help="Gene to correlate with age")
    parser.add_argument("--age-key", default="age_or_mean_of_age_range", 
                        type=str, help="Column in .obs with donor age")
    parser.add_argument("--donor-key", default="donor_id", type=str,
                        help="Column in .obs for donor ID")
    parser.add_argument("--celltype-key", default="cell_type", 
                        type=str, help="Column in .obs for cell type")
    return parser.parse_args()


def load_adata():
    adata = sc.read_h5ad(here("inputs/hlca/_m/hlca_core.h5ad"))
    symbols = adata.var["feature_name"].astype(str)
    adata.var["ensembl_id"] = adata.var_names
    adata.var_names = symbols
    adata.var.index = symbols

    # Fix annotation
    required_cats = [
        "Vascular smooth muscle", "Mesothelium", "Myofibroblasts",
        "AT1", "AT2", "Hematopoietic stem cells", "Lymphatic EC",
        'Mast cells', "Monocyte-derived Mph",
    ]

    if not isinstance(adata.obs["ann_level_4"].dtype, CategoricalDtype):
        adata.obs["ann_level_4"] = adata.obs["ann_level_4"].astype("category")

    for cat in required_cats:
        if cat not in adata.obs["ann_level_4"].cat.categories:
            adata.obs["ann_level_4"] = adata.obs["ann_level_4"].cat.add_categories([cat])

    # Handle annotation issues
    vsm_clusters = {"Smooth muscle", "Smooth muscle FAM83D+",
                    "SM activated stress response"}
    lec_clusters = {"Lymphatic EC differentiating",
                    "Lymphatic EC mature",
                    "Lymphatic EC proliferating"}
    fixes = {
        "Vascular smooth muscle": vsm_clusters,
        "Mesothelium": {"Mesothelium"},
        "Myofibroblasts": {"Myofibroblasts"},
        "AT1": {"AT1"}, "AT2": {"AT2"},
        "Hematopoietic stem cells": {"Hematopoietic stem cells"},
        "Lymphatic EC": lec_clusters,
        "Mast cells": {"Mast cells"},
        "Monocyte-derived Mph": {"Monocyte-derived Mph"},
    }

    for label, fine_clusters in fixes.items():
        cond = (
            adata.obs["ann_finest_level"].isin(fine_clusters)
            & (adata.obs["ann_level_4"].isna() |
               (adata.obs["ann_level_4"].isin(["None", ""])))
        )
        adata.obs.loc[cond, "ann_level_4"] = label

    # Update annotation columns
    obs_map = {
        "subclusters": "ann_finest_level",
        "cell_type": "ann_level_4",
        "clusters": "ann_level_4",
        "compartment": "ann_level_1",
        "patient": "donor_id",
    }
    for new, old in obs_map.items():
        if old in adata.obs:
            adata.obs[new] = adata.obs[old]
        else:
            print(f"Warning: '{old}' not found in obs; skipping {new} mapping.")

    # Filter studies (>= 20 cells)
    study_counts = adata.obs["study"].value_counts()
    valid_studies = study_counts[study_counts >= 20].index
    adata = adata[adata.obs["study"].isin(valid_studies)].copy()
    return adata


def ensure_logcounts(adata, threshold=20):
    """Check if adata.X is log-transformed. If not, apply log1p transformation."""
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else X
    max_val = np.percentile(X, 99)
    logging.info(f"99th percentile of expression matrix: {max_val:.2f}")

    if max_val >= threshold:
        logging.warning("Expression values appear to be raw counts. Applying log1p transformation.")
        adata.X = np.log1p(X)
    else:
        logging.info("Expression matrix appears to be log-transformed (logcounts).")
    return adata


def compute_donor_averages(adata, gene, age_key, donor_key):
    df = adata.obs.copy()
    df[gene] = adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else adata[:, gene].X.flatten()
    donor_df = df.groupby(donor_key, observed=False)\
                 .agg({gene: "mean", age_key: "first"}).dropna()
    return donor_df


def plot_correlation(df, gene, age_key, title, output_path):
    rho, pval = spearmanr(df[age_key], df[gene])
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.regplot(data=df, x=age_key, y=gene, scatter_kws={"s": 40}, 
                ax=ax, ci=None, line_kws={"color": "red"})
    ax.set_title(f"{title}\nSpearman r={rho:.2f}, p={pval:.2e}")
    sns.despine()
    fig.tight_layout()
    fig.savefig(output_path.with_suffix(".png"), dpi=300)
    fig.savefig(output_path.with_suffix(".pdf"))
    plt.close(fig)


def main():
    args = parse_args()
    configure_logging()

    args.outdir.mkdir(parents=True, exist_ok=True)
    adata = load_adata()
    adata = ensure_logcounts(adata)

    logging.info(f"Analyzing gene: {args.gene}")

    # Global correlation
    logging.info("Computing global donor-level correlation...")
    global_df = compute_donor_averages(adata, args.gene, args.age_key, args.donor_key)
    plot_correlation(global_df, args.gene, args.age_key, f"Global: Age vs {args.gene}", 
                     args.outdir / "global_correlation")

    # Per cell-type correlation
    logging.info("Computing per cell-type donor-level correlation...")
    results = []
    for cell_type in adata.obs[args.celltype_key].unique():
        subset = adata[adata.obs[args.celltype_key] == cell_type]
        df = compute_donor_averages(subset, args.gene, args.age_key, args.donor_key)
        if len(df) < 3:
            logging.warning(f"Skipping {cell_type} (too few donors)")
            continue
        rho, pval = spearmanr(df[args.age_key], df[args.gene])
        results.append({
            "cell_type": cell_type, "n_donors": len(df),
            "spearman_rho": rho, "p_value": pval
        })
        if pval < 0.05:
            plot_correlation(df, args.gene, args.age_key, f"{cell_type}",
                             args.outdir / f"{cell_type}_correlation")

    # Save summary table
    result_df = pd.DataFrame(results).sort_values("p_value")
    fdr_corrected, _ = fdrcorrection(result_df["p_value"])
    result_df["fdr"] = fdr_corrected
    result_df.to_csv(args.outdir / "age_correlation_summary.tsv", sep="\t", index=False)

    # Session info
    session_info.show()


if __name__ == "__main__":
    main()
