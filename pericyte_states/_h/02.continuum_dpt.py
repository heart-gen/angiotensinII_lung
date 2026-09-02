"""
Expression-only state continuum for pericytes (no RNA velocity / no BAM).

Tests whether lung pericytes occupy a continuum from vascular-supporting to
activated/fibrotic-like states using diffusion pseudotime (DPT) + PAGA on the
harmonized representation. We DO NOT claim a temporal direction: the ordering
reflects transcriptional similarity, rooted at the vascular-stabilizing state
purely as a reference pole.

Outputs:
  - pericyte_continuum.h5ad         (obs gains dpt_pseudotime)
  - continuum_metadata.tsv.gz       (per-cell pseudotime + scores for donor stats)
  - correlation of pseudotime vs state/AGTR1 scores (cell-level + donor-level)
  - UMAP / diffmap colored by pseudotime; PAGA graph; trend scatters
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from scipy import sparse, stats
from scipy.stats import spearmanr, rankdata
from anndata import AnnData
import matplotlib.pyplot as plt

sns.set_context("talk")
sns.set_style("whitegrid")

# Curated-program module scores (annotation scores from 00.state_discovery.py).
# Study is handled by the integrated embedding the DPT runs on, so scores are used
# at face value; the total-counts partial correlation guards against the trend
# being a sequencing-depth artifact.
TREND_COLS = [
    "vascular_stabilizing_score", "inflammatory_score",
    "synthetic_contractile_score", "activated_migratory_score",
    "fibroblast_like_score", "basement_membrane_score",
    "AGTR1_expr", "AGTR1_scvi",
]

# scVI-denoised AGTR1, merged from the same table and model that
# pericyte_states/_h/03.agtr1_lenses.R uses as its robust lens. Raw AGTR1 and the
# program scores both scale with capture depth, so a raw-only continuum trend
# cannot distinguish a receptor gradient from a depth gradient -- which is the
# exact confound 03.agtr1_lenses.R showed reverses AGTR1's apparent program bias.
# The denoised lens is carried alongside so figure_pericyte_layer panel F shows
# both and the reader can see the raw trend for what it is.
DENOISE_DEFAULT = ("../../localization/airspace_analysis/_m/airspace/"
                   "pericytes_airspace_denoising.tsv")


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad (from 00.state_discovery.py)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--use-rep", default="X_pca_harmony")
    p.add_argument("--neighbors", type=int, default=30)
    p.add_argument("--n-dcs", type=int, default=10)
    p.add_argument("--root-state", default="vascular_stabilizing")
    p.add_argument("--denoise", type=Path, default=Path(DENOISE_DEFAULT),
                   help="pericytes_airspace_denoising.tsv (scVI-denoised AGTR1)")
    p.add_argument("--den-model", default="Pericyte-only-trained",
                   help="which scVI model in --denoise to use as the denoised lens")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def save_figure(fig, base: Path):
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def pick_root(adata: AnnData, root_state: str, rep: str) -> int:
    """Root = cell of root_state closest to that state's centroid in latent space.

    KNOWN AND UNRESOLVED (recorded 2026-09-01, not fixed here): `--root-state`
    defaults to the PROGRAM name `vascular_stabilizing`, but `pericyte_state` holds
    the numeric Leiden cluster id, so the match is empty and the root silently falls
    back to the global PC1 minimum. Every sign in the continuum results -- BM
    falling, AGTR1 rising, the switch index -- is therefore root-dependent and
    currently set by PC1, not by the intended biology. The fallback is left in
    place deliberately so this run reproduces the published pseudotime; it is now
    logged at ERROR level and written to root_selection.tsv so it cannot be
    mistaken for the intended behaviour. Fixing it changes every published rho and
    is a separate, deliberate decision.
    """
    mask = (adata.obs["pericyte_state"].astype(str) == str(root_state)).to_numpy()
    if mask.sum() == 0:
        logging.error(
            "ROOT FALLBACK IN EFFECT: no cells match --root-state=%r in "
            "obs['pericyte_state'] (values: %s). Using the global PC1 minimum "
            "instead. Every continuum sign below is PC1-rooted, NOT rooted on the "
            "requested state.", root_state,
            ", ".join(map(str, pd.unique(adata.obs["pericyte_state"])[:10])))
        return int(np.argmin(adata.obsm[rep][:, 0]))
    X = adata.obsm[rep]
    centroid = X[mask].mean(axis=0)
    d = np.linalg.norm(X - centroid, axis=1)
    d[~mask] = np.inf
    return int(np.argmin(d))


def run_dpt(adata, rep, neighbors, n_dcs, root_idx, seed):
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=neighbors, random_state=seed)
    sc.tl.diffmap(adata, n_comps=max(n_dcs, 15))
    adata.uns["iroot"] = root_idx
    sc.tl.dpt(adata, n_dcs=n_dcs)
    sc.tl.paga(adata, groups="pericyte_state")


def ensure_total_counts(adata):
    """Per-cell depth used as the partial-correlation control."""
    if "total_counts" in adata.obs.columns:
        return
    for key in ("n_counts",):
        if key in adata.obs.columns:
            adata.obs["total_counts"] = np.asarray(adata.obs[key], dtype=float)
            return
    counts = adata.layers["counts"] if "counts" in adata.layers else adata.X
    tot = counts.sum(axis=1)
    adata.obs["total_counts"] = np.asarray(tot).ravel().astype(float)


def partial_spearman(x, y, z):
    """Spearman partial correlation of x,y controlling for z: rank all three,
    residualize ranks of x and y on ranks of z, correlate residuals.
    Returns (rho, p) with df = n - 3. All scores fall together along DPT when
    pseudotime tracks library depth; partialling z out tests for a depth-free
    association."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 5:
        return np.nan, np.nan
    rx, ry, rz = rankdata(x[m]), rankdata(y[m]), rankdata(z[m])
    Z = np.column_stack([np.ones(n), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    ex -= ex.mean(); ey -= ey.mean()
    denom = np.sqrt((ex @ ex) * (ey @ ey))
    if denom == 0:
        return np.nan, np.nan
    r = float((ex @ ey) / denom)
    df = n - 3
    if df <= 0 or abs(r) >= 1:
        return r, np.nan
    t = r * np.sqrt(df / (1 - r ** 2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p


def correlate_trends(adata, outdir: Path):
    ensure_total_counts(adata)
    pt = adata.obs["dpt_pseudotime"].to_numpy()
    tc = np.log10(adata.obs["total_counts"].to_numpy() + 1.0)
    rows = []
    for col in TREND_COLS:
        if col not in adata.obs:
            continue
        y = adata.obs[col].to_numpy()
        rho, pval = spearmanr(pt, y, nan_policy="omit")
        prho, ppval = partial_spearman(pt, y, tc)
        rows.append({"feature": col, "level": "cell", "spearman_rho": rho,
                     "p_value": pval, "partial_rho": prho, "partial_p": ppval,
                     "n": np.isfinite(pt).sum()})
    # donor-level (mean per donor): more conservative, donor-aware
    df = adata.obs.copy()
    agg = {**{"dpt_pseudotime": "mean", "total_counts": "mean"},
           **{c: "mean" for c in TREND_COLS if c in df.columns}}
    donor = df.groupby("donor_id", observed=True).agg(agg)
    donor = donor.dropna(subset=["dpt_pseudotime"])
    tcd = np.log10(donor["total_counts"].to_numpy() + 1.0)
    for col in TREND_COLS:
        if col not in donor.columns:
            continue
        sub = donor[["dpt_pseudotime", col, "total_counts"]].dropna()
        if sub.shape[0] < 5:
            continue
        rho, pval = spearmanr(sub["dpt_pseudotime"], sub[col])
        prho, ppval = partial_spearman(
            sub["dpt_pseudotime"].to_numpy(), sub[col].to_numpy(),
            np.log10(sub["total_counts"].to_numpy() + 1.0))
        rows.append({"feature": col, "level": "donor", "spearman_rho": rho,
                     "p_value": pval, "partial_rho": prho, "partial_p": ppval,
                     "n": sub.shape[0]})
    out = pd.DataFrame(rows)
    out.to_csv(outdir / "pseudotime_trend_correlations.tsv", sep="\t", index=False)
    return out


def plot_trends(adata, outdir: Path):
    coords = adata.obsm["X_umap"]
    pt = adata.obs["dpt_pseudotime"].to_numpy()
    fig, ax = plt.subplots(figsize=(6, 5))
    sca = ax.scatter(coords[:, 0], coords[:, 1], c=pt, s=6, linewidths=0, cmap="magma")
    fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title("DPT pseudotime (continuum, not time)")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    save_figure(fig, outdir / "umap_pseudotime")

    # Trend scatters with lowess
    for col in TREND_COLS:
        if col not in adata.obs:
            continue
        fig, ax = plt.subplots(figsize=(5, 4))
        sns.regplot(x=pt, y=adata.obs[col].to_numpy(), lowess=True, ax=ax,
                    scatter_kws=dict(s=4, alpha=0.2), line_kws=dict(color="red"))
        ax.set_xlabel("DPT pseudotime"); ax.set_ylabel(col)
        save_figure(fig, outdir / f"trend_{col}")

    # PAGA graph
    try:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.paga(adata, color="pericyte_state", ax=ax, show=False)
        save_figure(fig, outdir / "paga_states")
    except Exception as e:
        logging.warning(f"PAGA plot failed: {e}")


def add_denoised_agtr1(adata, path: Path, model: str):
    """Attach the scVI-denoised AGTR1 lens as obs['AGTR1_scvi'].

    Missing is tolerated and warned about, not fatal: the continuum result itself
    does not depend on it. But figure_pericyte_layer panel F asserts the row is
    present, so a silent absence there is caught downstream rather than here.
    """
    if not Path(path).exists():
        logging.warning(f"denoising table not found at {path}; the denoised AGTR1 "
                        "trend will be absent from pseudotime_trend_correlations.tsv")
        return
    den = pd.read_csv(path, sep="\t", usecols=["index", "Model", "AGTR1_scvi"])
    den = den.loc[den["Model"] == model]
    if den["index"].duplicated().any():
        raise ValueError(f"denoising table has duplicate cells for model {model}")
    den = den.set_index("index")["AGTR1_scvi"]
    adata.obs["AGTR1_scvi"] = den.reindex(adata.obs_names).to_numpy()
    n = int(np.isfinite(adata.obs["AGTR1_scvi"]).sum())
    logging.info(f"Denoised AGTR1 ({model}): {n}/{adata.n_obs} cells "
                 f"({100 * n / adata.n_obs:.1f}%)")


def main():
    args = parse_args()
    configure_logging()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if "AGTR1_expr" not in adata.obs and "AGTR1" in adata.var_names:
        expr = adata[:, "AGTR1"].layers["logcounts"]
        adata.obs["AGTR1_expr"] = (expr.toarray().ravel()
                                   if sparse.issparse(expr) else np.asarray(expr).ravel())
    add_denoised_agtr1(adata, args.denoise, args.den_model)

    root_idx = pick_root(adata, args.root_state, args.use_rep)
    root_matched = bool((adata.obs["pericyte_state"].astype(str)
                         == str(args.root_state)).any())
    logging.info(f"Root cell index={root_idx} (requested state={args.root_state}, "
                 f"matched={root_matched})")
    pd.DataFrame([{"requested_root_state": args.root_state,
                   "root_state_matched": root_matched,
                   "root_cell_index": root_idx,
                   "root_barcode": adata.obs_names[root_idx],
                   "method": "state centroid" if root_matched
                             else "GLOBAL PC1 MINIMUM (fallback)"}]).to_csv(
        outdir / "root_selection.tsv", sep="\t", index=False)
    run_dpt(adata, args.use_rep, args.neighbors, args.n_dcs, root_idx, args.seed)

    fig_dir = outdir / "figures"
    fig_dir.mkdir(exist_ok=True)
    trends = correlate_trends(adata, outdir)
    print(trends.to_string(index=False))
    plot_trends(adata, fig_dir)

    keep = [c for c in ["donor_id", "disease", "lung_condition", "pericyte_state",
                        "dpt_pseudotime", "total_counts"] + TREND_COLS
            if c in adata.obs.columns]
    adata.obs[keep].to_csv(outdir / "continuum_metadata.tsv.gz", sep="\t")

    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(adata.var.index.to_numpy(), adata.var[col].to_numpy()):
            adata.var.index.name = None
    adata.write(outdir / "pericyte_continuum.h5ad")

    session_info.show()


if __name__ == "__main__":
    main()
