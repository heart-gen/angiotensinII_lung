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
from anndata import AnnData
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf

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
                      classes=("AT1", "AT2", "EC aerocyte capillary", "EC general capillary")):
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


def fit_lmm(adata: AnnData, outdir: Path):
    """LMM: airspace_score ~ AGTR1_detect + age + sex + disease + (1|donor_id)."""
    df = adata.obs.dropna(subset=["airspace_score", "AGTR1_detect", "donor_id", 
                                  "sex", "age", "disease"]).copy()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(int)
    df["sex"] = df["sex"].astype("category")
    df["disease"] = df["disease"].astype("category")

    # Mixed model: donor_id as random effect
    model = smf.mixedlm("airspace_score ~ AGTR1_detect + age + sex + disease",
                        data=df, groups=df["donor_id"])
    result = model.fit(reml=True)
    res_df = pd.DataFrame({
        "param": result.params.index,
        "coef": result.params.values,
        "pval": result.pvalues.values
    })
    res_df.to_csv(outdir / "airspace_lmm_results.csv", index=False)
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
