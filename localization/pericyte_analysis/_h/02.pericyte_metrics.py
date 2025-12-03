"""
Compute donor mixing metrics (iLISI) + mixing plots for pericytes.
"""
import json
import logging
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from pathlib import Path
from anndata import AnnData
import matplotlib.pyplot as plt
from scib.metrics import ilisi_graph

sns.set_context("talk")
sns.set_style("whitegrid")

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata", required=True, type=Path,
                        help="pericyte_with_embeddings.h5ad")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--donor-key", default="donor_id")
    parser.add_argument("--use-rep", default="X_pca_harmony")
    return parser.parse_args()


def save_fig(fig, base: Path):
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_adata(path: Path) -> AnnData:
    return sc.read_h5ad(path)


def poisson_binom_pmf_fft(p: np.ndarray) -> np.ndarray:
    """Exact Poisson-binomial PMF via FFT."""
    p = np.asarray(p, dtype=float)
    n = p.size
    m = n + 1
    theta = 2 * np.pi * np.arange(m) / m
    exp_i_theta = np.exp(1j * theta)
    G = ((1 - p)[:, None] + p[:, None] * exp_i_theta[None, :]).prod(axis=0)
    pmf = np.real(np.fft.fft(G)) / m
    pmf = np.clip(pmf, 0, None)
    pmf /= pmf.sum()
    return pmf


def poisson_binom_pvalue(p: np.ndarray, k_obs: int, side="greater"):
    """Tail p-value under H0: K ~ Poisson-binomial(p)."""
    p = np.asarray(p, dtype=float)
    pmf = poisson_binom_pmf_fft(p)
    n = pmf.size - 1
    k_obs = int(np.clip(k_obs, 0, n))
    cdf = pmf.cumsum()
    if side == "greater":
        pval = pmf[k_obs:].sum()
    elif side == "less":
        pval = cdf[k_obs]
    else:
        p_less = cdf[k_obs]
        p_greater = pmf[k_obs:].sum()
        pval = min(1.0, 2.0 * min(p_less, p_greater))
    return {"pval": float(pval), "n": int(n)}


def matched_gene_zero_probs(X, gene_index: int, mean_tol=0.25,
                            cv_tol=0.25, top_n=150):
    # Per-gene mean & CV in pericytes
    gx = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    mu = gx.mean(0) + 1e-12
    sd = gx.std(0)
    cv = sd / mu
    mu_g, cv_g = float(mu[gene_index]), float(cv[gene_index])

    # Candidate matched genes
    keep = (np.abs(np.log1p(mu) - np.log1p(mu_g)) <= mean_tol) & (np.abs(cv - cv_g) <= cv_tol)
    keep[gene_index] = False
    matched = np.where(keep)[0]

    # per-cell expected zero probability
    p_i = (gx[:, matched] == 0).mean(1)
    return p_i, matched


def dropout_expectation_pvalue_matched(
        adata: AnnData, outdir: Path, gene="AGTR1", side="greater"
):
    Xp  = adata.X
    genes = adata.var_names
    gix = np.where(genes == gene)[0][0]

    # observed zero count
    GX = Xp.toarray() if hasattr(Xp, "toarray") else np.asarray(Xp)
    k_obs = int((GX[:, gix] == 0).sum())

    # expected per-cell dropout probabilities
    p_i, matched = matched_gene_zero_probs(Xp, gix)

    # p-value
    pval = poisson_binom_pvalue(p_i, k_obs, side=side)["pval"]

    # Save data
    fp = outdir / "dropout_expectation_results.tsv"
    rec = dict({
        "gene": gene, "n_cells": int(GX.shape[0]),
        "obs_zeros": k_obs, "exp_zeros": float(p_i.sum()),
        "ratio_obs_to_exp": float(k_obs / (p_i.sum() + 1e-9)),
        "pval": float(pval), "side": side, "method": "FFT",
        "n_matched": int(len(matched)),
    })
    pd.DataFrame([rec]).to_csv(fp, sep="\t", index=False)


def make_df(adata: AnnData, emb_key: str, donor_key: str):
    coords = adata.obsm[emb_key]
    df = pd.DataFrame(coords, columns=["dim1", "dim2"], index=adata.obs_names)
    df[donor_key] = adata.obs[donor_key]
    return df


def plot_donor_embedding(df: pd.DataFrame, donor_key: str,
                         title: str, base: Path):
    cat = pd.Categorical(df[donor_key])
    palette = sns.color_palette("tab20", len(cat.categories))

    fig, ax = plt.subplots(figsize=(6, 5))
    for col, cat_val in zip(palette, cat.categories):
        mask = cat == cat_val
        ax.scatter(df.loc[mask, "dim1"], df.loc[mask, "dim2"],
                   s=6, linewidths=0, color=col, label=str(cat_val))
    ax.legend().remove(); ax.set_title(title)
    save_fig(fig, base)


def compute_ilisi(adata: AnnData, donor_key: str, use_rep: str) -> float:
    mask = adata.obs[donor_key].notna()
    ad = adata[mask].copy()
    ad.obs[donor_key] = ad.obs[donor_key].astype(str)
    score = ilisi_graph(
        ad, batch_key=donor_key, use_rep=use_rep, type_="embed", n_cores=1
    )
    return float(score)


def write_json(path: Path, adata: AnnData, score: float,
               donor_key: str, use_rep: str):
    payload = {
        "embedding": {
            "use_rep": use_rep,
            "donor_key": donor_key,
            "n_cells": int(adata.n_obs),
            "n_donors": int(adata.obs[donor_key].nunique()),
            "iLISI_donor": None if not np.isfinite(score) else score,
        }
    }
    path.write_text(json.dumps(payload, indent=2))


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_adata(args.adata)

    # Check dropout
    dropout_expectation_pvalue_matched(
        adata, outdir, gene="AGTR1", side="greater"
    )

    # Make mixing plots
    mix_dir = outdir / "donor_mixing"
    mix_dir.mkdir(exist_ok=True)

    for emb in ["X_umap", "X_tsne"]:
        df = make_df(adata, emb, args.donor_key)
        plot_donor_embedding(df, args.donor_key, f"{emb}: donor mixing",
                             mix_dir / f"{emb}_donor")

    # Compute iLISI
    score = compute_ilisi(
        adata, donor_key=args.donor_key, use_rep=args.use_rep,
    )

    # Write JSON
    write_json(outdir / "pericyte_embedding.json",
               adata, score, args.donor_key, args.use_rep)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
