import scanpy as sc
import session_info
import pandas as pd
import logging, argparse
from pathlib import Path
from anndata import AnnData
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata-file", type=str,
                        default="airspace_pericyte_states.h5ad",
                        help="Filename of AnnData with airspace score")
    parser.add_argument("--outdir", type=Path, default=Path("results"),
                        help="Output directory for QC and updated AnnData")
    parser.add_argument("--cluster-key", type=str, default="subclusters",
                        help="obs column with fine cell types (default: subclusters)")
    parser.add_argument("--pericyte-label", type=str, default="Pericytes",
                        help="Value in cluster-key used for pericytes")
    return parser.parse_args()


def load_anndata(path: Path) -> AnnData:
    logging.info("Loading AnnData from %s", path)
    adata = sc.read_h5ad(path)
    return adata


def analyze_pericytes_airspace(
    df: pd.DataFrame, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", leiden_key: str = "leiden_pericytes",
    score_key: str = "airspace_score", imputed_key: str = "AGTR1_scvi",
    donor_key: str = "donor_id", min_cells_per_donor: int = 20,
) -> None:

    cols_needed = [cluster_key, leiden_key, score_key, imputed_key,
                   donor_key, "sex", "disease", "self_reported_ethnicity",
                   "age_or_mean_of_age_range"]
    cols_existing = [col for col in cols_needed if col in df.columns]
    peri = df[df[cluster_key].astype(str) == pericyte_label].copy()
    # Per donor summaries
    per_donor_counts = peri.groupby(donor_key, observed=False).size()
    keep_donors = per_donor_counts[per_donor_counts >= min_cells_per_donor].index
    dperi = peri.loc[peri[donor_key].isin(keep_donors)].copy()
    # Airspace vs imputed AGTR1
    new_df = dperi.dropna(subset=[score_key, leiden_key, imputed_key]).copy()
    formula = f"{score_key} ~ {imputed_key} + C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range"
    fit = smf.ols(formula, data=new_df).fit(cov_type="HC3") # unequal variances
    anova_table = anova_lm(fit, typ=2)
    anova_table.to_csv(outdir / f"anova_ols.{score_key}.{imputed_key}.tsv", sep="\t")
        
    # Airspace vs leident clusters
    valid = new_df.groupby(donor_key, observed=False)[leiden_key].nunique() >= 2
    df2 = new_df[new_df[donor_key].isin(valid[valid].index)].copy()
    if not df2.empty and df2[leiden_key].nunique() > 1 and df2[donor_key].nunique() > 1:
        formula = f"{score_key} ~ C({leiden_key}) + C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range"
        fit = smf.ols(formula, data=df2).fit(cov_type="HC3") # unequal variances
        anova_table = anova_lm(fit, typ=2)
        anova_table.to_csv(outdir / f"anova_{score_key}_ols.tsv", sep="\t")
        with open(outdir / f"donor_aware_ols_{score_key}.txt", "w") as fh:
            fh.write(fit.summary().as_text())
            fh.write("\n\nDonor-adjusted effect of Leiden:\n")
        pairwise = fit.t_test_pairwise(f"C({leiden_key})")
        summary  = pairwise.result_frame
        summary.to_csv(outdir / f"pairwise_leiden_{score_key}_ols.csv", index=False)


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading data.")
    adata_path = outdir / args.adata_file
    adata = load_anndata(adata_path)
    
    # Imputed denoising per donor
    impute_dir = outdir / "airspace"
    imputed_path = impute_dir / "pericytes_airspace_denoising.tsv"
    imputed_df = pd.read_csv(imputed_path, sep="\t", index_col=0)
    imputed_df = imputed_df[imputed_df["Model"] == "Pericyte-only-trained"].copy()

    obs = adata.obs.copy()
    df  = obs.merge(imputed_df, left_index=True, right_index=True, how="inner")
    analyze_pericytes_airspace(df, impute_dir)

    logging.info("Pericyte analysis complete.")

    # Session information
    session_info.show()
    

if __name__ == "__main__":
    main()
