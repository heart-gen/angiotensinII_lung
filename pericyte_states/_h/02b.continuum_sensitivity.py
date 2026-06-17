"""
Robustness of the pericyte expression continuum (DPT) to analyst choices.

02.continuum_dpt.py reports the headline result: along a diffusion-pseudotime
axis rooted at the vascular-stabilizing pole, the injury / activation / ECM and
AGTR1 scores rise monotonically. A fair reviewer asks whether that trend is an
artifact of the specific root cell, neighborhood size, number of diffusion
components, or the particular set of cells. This script re-runs DPT across a grid
of those choices and reports the distribution of the pseudotime-vs-score Spearman
correlations.

Two questions, two blocks of runs:
  (1) ROOT-pole invariance of the AXIS. Rooting at the opposite (injury) pole just
      reverses pseudotime, so the sign of rho flips but |rho| -- the strength of
      each score's association with the continuum -- should be preserved. We root
      at each state in turn (+ a latent PC1 extreme) at default params.
  (2) PARAMETER robustness at the canonical (vascular-stabilizing) root: sweep
      n_neighbors x n_dcs x subsample fraction (with seeds). Here both sign and
      magnitude of rho should be stable.

Output:
  continuum_sensitivity_runs.tsv      one row per (setting x feature): rho, p, n
  continuum_sensitivity_summary.tsv   per feature, over the canonical-root runs:
                                      mean/sd/min/max rho + sign-consistency
  figures/continuum_sensitivity.{pdf,png}
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

TREND_COLS = [
    "vascular_stabilizing_score", "inflammatory_score",
    "synthetic_contractile_score", "activated_migratory_score",
    "fibroblast_like_score", "AGTR1_expr",
]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad (from 00.state_discovery.py)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--use-rep", default="X_pca_harmony")
    p.add_argument("--canonical-root", default="vascular_stabilizing")
    return p.parse_args()


def ensure_agtr1(adata):
    if "AGTR1_expr" not in adata.obs and "AGTR1" in adata.var_names:
        expr = adata[:, "AGTR1"].layers["logcounts"]
        adata.obs["AGTR1_expr"] = (expr.toarray().ravel()
                                   if sparse.issparse(expr) else np.asarray(expr).ravel())


def pick_root(obs_states, X, root_state):
    """Cell of root_state nearest that state's latent centroid; PC1 extremes for
    the two pseudo-roots pc1min / pc1max."""
    if root_state == "pc1min":
        return int(np.argmin(X[:, 0]))
    if root_state == "pc1max":
        return int(np.argmax(X[:, 0]))
    mask = (obs_states == root_state).to_numpy()
    if mask.sum() == 0:
        return int(np.argmin(X[:, 0]))
    centroid = X[mask].mean(axis=0)
    d = np.linalg.norm(X - centroid, axis=1)
    d[~mask] = np.inf
    return int(np.argmin(d))


def one_run(adata, rep, root_state, neighbors, n_dcs, frac, seed):
    """Subsample (if frac<1), run neighbors->diffmap->dpt, return per-feature rho."""
    if frac < 1.0:
        rng = np.random.default_rng(seed)
        idx = rng.choice(adata.n_obs, size=int(round(frac * adata.n_obs)), replace=False)
        sub = adata[np.sort(idx)].copy()
    else:
        sub = adata.copy()
    root_idx = pick_root(sub.obs["pericyte_state"], sub.obsm[rep], root_state)
    sc.pp.neighbors(sub, use_rep=rep, n_neighbors=neighbors, random_state=seed)
    sc.tl.diffmap(sub, n_comps=max(n_dcs, 15))
    sub.uns["iroot"] = root_idx
    sc.tl.dpt(sub, n_dcs=n_dcs)
    pt = sub.obs["dpt_pseudotime"].to_numpy()
    rows = []
    for col in TREND_COLS:
        if col not in sub.obs:
            continue
        rho, pval = spearmanr(pt, sub.obs[col].to_numpy(), nan_policy="omit")
        rows.append(dict(feature=col, spearman_rho=rho, p_value=pval, n=sub.n_obs))
    return rows


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.outdir / "figures"; fig_dir.mkdir(exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    ensure_agtr1(adata)
    rep = args.use_rep
    states = list(pd.unique(adata.obs["pericyte_state"]))

    runs = []

    # ---- block (2): parameter robustness at the canonical root -------------
    for neighbors in (15, 30, 50):
        for n_dcs in (10, 15):
            for frac, seed in [(1.0, 13), (0.8, 13), (0.8, 7)]:
                logging.info(f"[param] root={args.canonical_root} nbr={neighbors} "
                             f"n_dcs={n_dcs} frac={frac} seed={seed}")
                for r in one_run(adata, rep, args.canonical_root, neighbors, n_dcs, frac, seed):
                    r.update(block="param", root=args.canonical_root,
                             neighbors=neighbors, n_dcs=n_dcs, frac=frac, seed=seed)
                    runs.append(r)

    # ---- block (1): root-pole invariance (default params) ------------------
    alt_roots = [s for s in states if s != args.canonical_root] + ["pc1min", "pc1max"]
    for root_state in alt_roots:
        logging.info(f"[root] root={root_state} (default params)")
        for r in one_run(adata, rep, root_state, 30, 10, 1.0, 13):
            r.update(block="root", root=str(root_state),
                     neighbors=30, n_dcs=10, frac=1.0, seed=13)
            runs.append(r)

    df = pd.DataFrame(runs)
    df.to_csv(args.outdir / "continuum_sensitivity_runs.tsv", sep="\t", index=False)

    # ---- summary over canonical-root parameter sweep -----------------------
    can = df[df["block"] == "param"]
    summ = (can.groupby("feature")["spearman_rho"]
            .agg(n_settings="count", mean_rho="mean", sd_rho="std",
                 min_rho="min", max_rho="max")
            .reset_index())
    # fraction of canonical-root settings with the same sign as the mean
    sign = (can.assign(sgn=np.sign(can["spearman_rho"]))
            .groupby("feature")["sgn"]
            .apply(lambda s: (s == np.sign(s.mean())).mean())
            .rename("sign_consistency").reset_index())
    summ = summ.merge(sign, on="feature")
    summ.to_csv(args.outdir / "continuum_sensitivity_summary.tsv", sep="\t", index=False)
    print("\n== Canonical-root parameter robustness (Spearman rho vs DPT) ==")
    print(summ.to_string(index=False))

    # ---- figure: rho across all settings, faceted by block -----------------
    order = [c for c in TREND_COLS if c in set(df["feature"])]
    g = sns.catplot(data=df, x="spearman_rho", y="feature", hue="block",
                    order=order, kind="strip", dodge=True, height=4.5, aspect=1.5,
                    s=7, alpha=0.7)
    g.refline(x=0, color="grey", lw=1)
    g.set_axis_labels("Spearman rho (pseudotime vs score)", "")
    g.fig.suptitle("DPT continuum: robustness to root / neighbors / n_dcs / subsample",
                   y=1.02, fontsize=11)
    g.savefig(fig_dir / "continuum_sensitivity.png", dpi=300, bbox_inches="tight")
    g.savefig(fig_dir / "continuum_sensitivity.pdf", bbox_inches="tight")
    plt.close(g.fig)

    session_info.show()


if __name__ == "__main__":
    main()
