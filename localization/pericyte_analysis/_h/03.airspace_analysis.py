"""
Airspace proximity scoring for pericytes:
- Compute class centroids in latent space
- Airspace score = mean cosine similarity to {AT1, AT2, EC aerocyte, EC general capillary}
- Compare AGTR1+ vs AGTR1- with donor random intercept (LMM)
Outputs:
  adata.obs['airspace_score']
  lmm_results.csv
"""
import warnings
import logging, argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from scipy import sparse
from anndata import AnnData
import statsmodels.api as sm
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from pandas.api.types import CategoricalDtype
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.regression.mixed_linear_model import MixedLMResults

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata", type=Path, required=True,
                        help="HLCA pericyte/airspace subset AnnData")
    parser.add_argument("--use-rep", default="X_pca_harmony",
                        help="Latent representation for centroids")
    parser.add_argument("--cluster-key", default="subclusters")
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()



def add_agtr1_expression(adata: AnnData, gene="AGTR1"):
    if gene not in adata.var_names:
        logging.warning(f"{gene} not in var_names")
        return
    expr = adata[:, gene].layers["logcounts"]
    expr = expr.toarray().ravel() if sparse.issparse(expr) else np.asarray(expr).ravel()
    adata.obs[f"{gene}_expr"] = expr
    adata.obs[f"{gene}_detect"] = (expr > 0).astype(int)


def load_adata(path: Path) -> AnnData:
    adata = sc.read_h5ad(path)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X

    # Change var_names to gene names
    symbols = adata.var["feature_name"].astype(str)
    adata.var["ensembl_id"] = adata.var_names
    adata.var_names = adata.var.index = symbols

    adata.raw = None
    adata.var_names_make_unique()

    if "AGTR1_detect" not in adata.obs:
        add_agtr1_expression(adata, gene="AGTR1")

    return adata


def compute_centroids(adata: AnnData, rep: str, key: str,
                      classes=("AT1", "AT2", "EC aerocyte capillary",
                               "EC general capillary")):
    """Compute mean latent vector per class."""
    X = adata.obsm[rep]
    centroids = {}
    for cls in classes:
        mask = adata.obs[key] == cls
        if mask.sum() == 0:
            logging.warning(f"No cells for {cls}")
            continue
        centroids[cls] = X[mask].mean(axis=0)
    return centroids


def compute_airspace_scores(adata: AnnData, rep: str, centroids: dict,
                            pericyte_label="Pericytes", key="subclusters"):
    """Mean cosine similarity to AT1/AT2/EC centroids."""
    X = adata.obsm[rep]
    pericyte_mask = adata.obs[key] == pericyte_label

    scores = np.zeros(X.shape[0]) * np.nan
    cx = np.stack(list(centroids.values()))
    cx /= np.linalg.norm(cx, axis=1, keepdims=True)

    for idx in np.where(pericyte_mask)[0]:
        v = X[idx]
        v = v / np.linalg.norm(v)
        scores[idx] = (cx @ v).mean()

    adata.obs["airspace_score"] = scores
    return adata


def fit_lmm(adata: AnnData, outdir: Path, pericyte_label="Pericytes", key="subclusters"):
    """
    Donor-level analysis for airspace score analysis in pericytes.
    """
    # Filter to pericytes
    df = adata.obs.copy()
    df = df[df[key] == pericyte_label].copy()

    # Drop incomplete rows
    df = df.dropna(subset=[
        "airspace_score", "AGTR1_detect", "donor_id", 
        "sex", "age_or_mean_of_age_range",
    ]).copy()

    df["AGTR1_detect"] = df["AGTR1_detect"].astype(int)
    df["sex"] = df["sex"].astype("category")
    df["age"] = df["age_or_mean_of_age_range"].astype(int)

    # Keep donors with both AGTR1+ AND AGTR1-
    donor_levels = df.groupby("donor_id", observed=False)["AGTR1_detect"].nunique()
    keep_donors = donor_levels[donor_levels >= 2].index
    df = df[df["donor_id"].isin(keep_donors)].copy()

    # Aggregate to donor level
    donor_df = df.groupby(["donor_id"], observed=False)\
                 .agg(
                     mean_airspace_score=("airspace_score", "mean"),
                     frac_AGTR1_pos=("AGTR1_detect", "mean"),
                     n_cells=("AGTR1_detect", "size"),
                     age=("age", "first"),
                     sex=("sex", "first"),
                 ).reset_index()
    donor_df.to_csv(outdir / "airspace_donor_summary.csv", index=False)

    # Fit OLS at donor level
    formula = "mean_airspace_score ~ frac_AGTR1_pos + age + C(sex)"
    ols_res = smf.ols(formula, data=donor_df).fit()

    # Save coefficient table
    effect_row = {
        "term": "frac_AGTR1_pos",
        "estimate": ols_res.params.get("frac_AGTR1_pos", float("nan")),
        "se": ols_res.bse.get("frac_AGTR1_pos", float("nan")),
        "pval": ols_res.pvalues.get("frac_AGTR1_pos", float("nan")),
    }
    pd.DataFrame([effect_row]).to_csv(
        outdir / "airspace_effect_AGTR1.csv", index=False
    )

    # Diagnostic plot
    plt.figure(figsize=(5, 4))
    sns.scatterplot(data=donor_df, x="frac_AGTR1_pos", y="mean_airspace_score")
    sns.regplot(
        data=donor_df, x="frac_AGTR1_pos", y="mean_airspace_score",
        scatter=False
    )
    plt.xlabel("Fraction AGTR1+ Pericytes\n(per donor)")
    plt.ylabel("Mean Airspace Proximity Score\n(per donor)")
    plt.tight_layout()
    plt.savefig(outdir / "airspace_score_scatter.png", dpi=300)
    plt.savefig(outdir / "airspace_score_scatter.pdf")
    plt.close()

    # Save text summary
    with open(outdir / "airspace_model_summary.txt", "w") as fh:
        fh.write(ols_res.summary().as_text())

    return ols_res


def fit_lmm_per_cell(
    adata: AnnData, outdir: Path, pericyte_label="Pericytes", key="subclusters",
    include_dataset_vc=True, min_cells_per_class_per_donor=1, extra_covariates=None,
):
    outdir.mkdir(exist_ok=True, parents=True)

    df = adata.obs.copy()
    df = df[df[key] == pericyte_label].copy()

    # standardize
    if "age_or_mean_of_age_range" in df.columns and "age" not in df.columns:
        df["age"] = df["age_or_mean_of_age_range"]

    required = ["airspace_score", "AGTR1_detect", "donor_id", "sex", "age"]
    df = df.dropna(subset=required).copy()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(int)
    df["sex"] = df["sex"].astype("category")
    df["age"] = df["age"].astype(float)

    if "dataset" in df.columns:
        df["dataset"] = df["dataset"].astype("category")

    # QC covariates (optional)
    qc_covs = [c for c in ["n_counts", "n_genes", "pct_mito"] if c in df.columns]
    if "pct_mito" in qc_covs and df["pct_mito"].max() > 1:
        df["pct_mito"] /= 100.0

    # filter donors
    if min_cells_per_class_per_donor > 0:
        counts = (
            df.groupby(["donor_id", "AGTR1_detect"]).size().unstack(fill_value=0)
        )
        keep = counts[
            (counts.get(0, 0) >= min_cells_per_class_per_donor) &
            (counts.get(1, 0) >= min_cells_per_class_per_donor)
        ].index
        df = df[df["donor_id"].isin(keep)].copy()

    df.groupby(["donor_id", "AGTR1_detect"]).size().rename("n_cells").reset_index() \
      .to_csv(outdir / "per_cell_donor_counts.csv", index=False)

    # extra covariates
    extra_covariates = extra_covariates or []
    extra_covs = []
    for c in extra_covariates:
        if c in df.columns:
            if df[c].dtype.name in ("object", "string"):
                df[c] = df[c].astype("category")
            extra_covs.append(c)

    # formula
    terms = ["AGTR1_detect", "age", "C(sex)"]
    terms.extend(qc_covs)
    for c in extra_covs:
        if isinstance(df[c].dtype, CategoricalDtype):
            terms.append(f"C({c})")
        else:
            terms.append(c)
    formula = "airspace_score ~ " + " + ".join(terms)

    # MixedLM
    groups = df["donor_id"].astype("category")
    vc_formula = None
    if include_dataset_vc and "dataset" in df.columns and df["dataset"].nunique() > 1:
        vc_formula = {"dataset_vc": "0 + C(dataset)"}

    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            md = sm.MixedLM.from_formula(
                formula=formula, data=df, groups=groups,
                re_formula="1", vc_formula=vc_formula,
            )
            model = md.fit(method="lbfgs", maxiter=500, disp=False)
            fallback = False
    except Exception:
        model = sm.OLS.from_formula(formula, data=df).fit(
            cov_type="cluster", cov_kwds={"groups": df["donor_id"]}
        )
        fallback = True

    # Save results
    params = model.params
    bse = model.bse
    pvals = model.pvalues

    try:
        ci = model.conf_int().rename(columns={0: "ci_lower", 1: "ci_upper"})
        ci_agtr1 = ci.loc["AGTR1_detect"].to_dict()
    except Exception:
        ci_agtr1 = {"ci_lower": np.nan, "ci_upper": np.nan}

    effect = {
        "term": "AGTR1_detect",
        "estimate": params.get("AGTR1_detect", np.nan),
        "se": bse.get("AGTR1_detect", np.nan),
        "pval": pvals.get("AGTR1_detect", np.nan),
        **ci_agtr1,
        "n_cells": df.shape[0],
        "n_donors": df["donor_id"].nunique(),
        "fallback": fallback,
        "formula": formula,
    }

    pd.DataFrame([effect]).to_csv(outdir / "airspace_effect_AGTR1_per_cell.csv", index=False)

    if hasattr(model, "random_effects"):
        pd.DataFrame.from_dict(model.random_effects, orient="index") \
            .to_csv(outdir / "airspace_lmm_per_cell_random_effects.csv")

    with open(outdir / "airspace_lmm_per_cell_summary.txt", "w") as fh:
        fh.write(model.summary().as_text())

    return model


def plot_ridge(adata: AnnData, outdir: Path):
    df = adata.obs.loc[adata.obs["subclusters"] == "Pericytes"].copy()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(str)

    plt.figure(figsize=(6, 5))
    sns.kdeplot(
        data=df, x="airspace_score", hue="AGTR1_detect",
        fill=True, common_norm=False, alpha=0.5
    )
    plt.xlabel("Airspace Proximity Score")
    plt.tight_layout()
    plt.savefig(outdir / "airspace_ridge.png", dpi=300)
    plt.savefig(outdir / "airspace_ridge.pdf")
    plt.close()


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    adata = load_adata(args.adata)

    # Compute centroids
    centroids = compute_centroids(
        adata, rep=args.use_rep, key=args.cluster_key,
    )

    # Airspace score
    adata = compute_airspace_scores(
        adata, rep=args.use_rep, centroids=centroids,
        pericyte_label="Pericytes",
        key=args.cluster_key,
    )

    # LMM analysis
    airspace_dir = outdir / "airspace"
    airspace_dir.mkdir(exist_ok=True)

    # Donor-level
    result = fit_lmm(adata, airspace_dir)
    print(result.summary())

    # Per-cell
    cell_res = fit_lmm_per_cell(
        adata, outdir=airspace_dir, pericyte_label="Pericytes",
        key=args.cluster_key, extra_covariates=["pericyte_state"],
        include_dataset_vc=True,
    )
    print(cell_res.summary())

    # Plots
    plot_ridge(adata, airspace_dir)

    # Save AnnData with airspace scores
    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(
                adata.var.index.to_numpy(),
                adata.var[col].to_numpy()
        ):
            adata.var.index.name = None

    adata.write(outdir / "pericytes_with_airspace_score.h5ad")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
