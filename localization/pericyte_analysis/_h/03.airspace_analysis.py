"""
Airspace proximity scoring for pericytes:
- Compute class centroids in latent space
- Airspace score = mean cosine similarity to {AT1, AT2, EC aerocyte, EC general capillary}
- Compare AGTR1+ vs AGTR1- with donor random intercept (LMM)
Outputs:
  adata.obs['airspace_score']
  lmm_results.csv
"""

import logging, argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from scipy import sparse
from anndata import AnnData
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from sklearn.metrics.pairwise import cosine_similarity

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
    # normalize centroid matrix for speed
    cx_norm = cx / np.linalg.norm(cx, axis=1, keepdims=True)

    for idx in np.where(pericyte_mask)[0]:
        v = X[idx]
        v_norm = v / np.linalg.norm(v)
        sims = cx_norm @ v_norm
        scores[idx] = sims.mean()

    adata.obs["airspace_score"] = scores
    return adata


def fit_lmm(adata: AnnData, outdir: Path, pericyte_label="Pericytes", key="subclusters"):
    """
    LMM: airspace_score ~ AGTR1_detect + age + sex + disease + (1|donor_id).
    Restrict to donors with both AGTR1+ and AGTR1- pericytes.
    """
    df = adata.obs.copy()
    df = df[df[key] == pericyte_label].copy()

    # Drop missing values
    df = df.dropna(subset=["airspace_score", "AGTR1_detect", "donor_id", 
                           "sex", "age_or_mean_of_age_range", "disease"]).copy()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(int)
    df["sex"] = df["sex"].astype("category")
    df["disease"] = df["disease"].astype("category")
    df["age"] = df["age_or_mean_of_age_range"].astype(int)

    # Keep donors with both AGTR1+ AND AGTR1- pericytes
    n_levels = df.groupby("donor_id", observed=False)["AGTR1_detect"].nunique()
    keep_donors = n_levels[n_levels >= 2].index
    df = df[df["donor_id"].isin(keep_donors)].copy()
    n_donors = df["donor_id"].nunique()

    # Save summary table
    donor_summary = (
        df.groupby(["donor_id", "AGTR1_detect"], observed=False)
          .size()
          .unstack(fill_value=0)
          .rename(columns={0: "AGTR1_neg_cells", 1: "AGTR1_pos_cells"})
    )
    donor_summary.to_csv(outdir / "airspace_lmm_donor_summary.csv")

    results = {}
    # Donor fixed-effect model
    formula = "airspace_score ~ AGTR1_detect + age + C(sex) + C(disease) + C(donor_id)"
    ols_res = smf.ols(formula, data=df).fit()
    results["model_type"] = "ols_donor_fixed_effects"
    results["summary"] = ols_res.summary().as_text()
        
    coef = ols_res.params.get("AGTR1_detect", np.nan)
    se = ols_res.bse.get("AGTR1_detect", np.nan)
    pval = ols_res.pvalues.get("AGTR1_detect", np.nan)

    est_df = pd.DataFrame(
        {"term": ["AGTR1_detect"], "estimate": [coef], "se": [se], "pval": [pval]}
    )
    est_df.to_csv(outdir / "airspace_ols_effect_AGTR1.csv", index=False)

    # Quick diagnostic plot
    plt.figure(figsize=(5, 4))
    sns.violinplot(
        data=df, x="AGTR1_detect", y="airspace_score", inner="quartile", cut=0,
    )
    plt.xticks([0, 1], ["AGTR1−", "AGTR1+"])
    plt.ylabel("Airspace score")
    plt.tight_layout()
    plt.savefig(outdir / "airspace_score_by_AGTR1_detect.png", dpi=300)
    plt.savefig(outdir / "airspace_score_by_AGTR1_detect.pdf")
    plt.close()

    # Save text summary
    with open(outdir / "airspace_model_summary.txt", "w") as fh:
        fh.write(f"Model type: {results['model_type']}\n\n")
        fh.write(results["summary"])
    return result


def plot_violin(adata: AnnData, outdir: Path):
    df = adata.obs.loc[adata.obs["subclusters"] == "Pericytes"].copy()
    plt.figure(figsize=(5, 5))
    sns.violinplot(data=df, x="AGTR1_detect", y="airspace_score")
    sns.stripplot(data=df, x="AGTR1_detect", y="airspace_score",
                  color="black", size=2, alpha=0.5)
    plt.xticks([0, 1], ["AGTR1−", "AGTR1+"])
    plt.tight_layout()
    plt.savefig(outdir / "airspace_violin.png", dpi=300)
    plt.savefig(outdir / "airspace_violin.pdf")
    plt.close()


def plot_ridge(adata: AnnData, outdir: Path):
    df = adata.obs.loc[adata.obs["subclusters"] == "Pericytes"].copy()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(str)

    plt.figure(figsize=(6, 5))
    sns.kdeplot(
        data=df, x="airspace_score", hue="AGTR1_detect",
        fill=True, common_norm=False, alpha=0.5
    )
    plt.xlabel("Airspace proximity score")
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

    # LMM
    airspace_dir = outdir / "airspace"
    airspace_dir.mkdir(exist_ok=True)
    
    result = fit_lmm(adata, airspace_dir)
    print(result.summary())

    # Plots
    plot_violin(adata, airspace_dir)
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
