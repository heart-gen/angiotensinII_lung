"""
Robustness of the pericyte expression continuum (DPT) to analyst choices.

02.continuum_dpt.py reports the headline result: along a diffusion-pseudotime
axis rooted at the vascular-stabilizing pole, the state / ECM and AGTR1 scores
order monotonically. A fair reviewer asks whether that ordering is an artifact of
the specific root cell, neighborhood size, number of diffusion components, or the
particular set of cells. This script re-runs DPT across a grid of those choices
and reports the distribution of the pseudotime-vs-score Spearman correlations.

Two questions, two blocks of runs:
  (1) ROOT-pole invariance of the AXIS. Rooting at the opposite pole just reverses
      pseudotime, so the sign of rho flips but |rho| -- the strength of each
      score's association with the continuum -- should be preserved. We root at
      each state PROGRAM in turn, at each Leiden CLUSTER in turn, and at the two
      latent PC1 extremes, all at default params.
  (2) PARAMETER robustness at the canonical (vascular-stabilizing) root: sweep
      n_neighbors x n_dcs x subsample fraction (with seeds). Here both sign and
      magnitude of rho should be stable.

FIXED 2026-09-07 (P1-3, second instance). This script carried its OWN copy of
`pick_root`, and that copy matched `--canonical-root` against `pericyte_state`
alone. `pericyte_state` holds Leiden ids ('0'..'5'); the canonical root is the
PROGRAM name `vascular_stabilizing`, which lives in `state_program`. The match was
therefore always empty and every run in block (2) silently fell back to the global
PC1 minimum -- a `basement_membrane` cell, the OPPOSITE pole. The runs table and
the log both then labelled those runs `root=vascular_stabilizing`, so the output
did not merely lose the root, it misreported it, and
`continuum_sensitivity_summary.tsv` -- computed over block (2) alone -- described
the wrong-rooted sweep. The old file is reproducible as the `root=pc1min` rows,
which were bit-identical to the "canonical" run at nbr=30/n_dcs=10/frac=1.0.

The structural fix is that root selection is no longer duplicated here: this
script imports `pick_root` from 02.continuum_dpt.py, which hard-fails on an
unmatched root. Only the two PC1 extremes are resolved locally, because there they
are deliberately requested pseudo-roots rather than a fallback. Every run also
records the column its root matched on and the `state_program` of the root cell
it actually used, so a mislabelled root cannot survive unread again.

Every run is scored twice: at the CELL level, and DONOR-aware (Spearman of the
per-donor mean pseudotime against the per-donor mean score). The donor level is the
conservative read -- it is the unit a reviewer will ask for -- and it is aggregated
exactly as the headline script does (`02.continuum_dpt.py::correlate_trends`): plain
per-donor means, no minimum-cell filter, so the canonical setting reproduces the
donor rows of `pseudotime_trend_correlations.tsv`.

Output:
  continuum_sensitivity_runs.tsv      one row per (setting x feature x level):
                                      rho, p, n, plus the root's provenance;
                                      `level` is "cell" or "donor"
  continuum_sensitivity_summary.tsv   per (feature, level), over the canonical-root
                                      runs: mean/sd/min/max rho + sign-consistency
  continuum_sensitivity_root_audit.tsv one row per distinct root: what was
                                      requested, what column matched, and which
                                      program the root cell belongs to
  figures/continuum_sensitivity.{pdf,png}
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import importlib.util
import logging, argparse, sys
from pathlib import Path
from scipy import sparse
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

HEADLINE = Path(__file__).resolve().parent / "02.continuum_dpt.py"

# The two latent extremes are the only roots resolved in this file. They are
# pseudo-roots on purpose -- "what if the analyst had picked no biological pole at
# all" -- not a fallback for a root that failed to match.
PC1_PSEUDOROOTS = {"pc1min": np.argmin, "pc1max": np.argmax}


def load_headline(path: Path):
    """Import 02.continuum_dpt.py by path (the leading digits make it an invalid
    module name, so a plain import will not work). Same idiom as
    00b.annotation_support.py::load_discovery.

    Importing rather than copying is the point: `pick_root` and the trend feature
    list must be the ones the headline actually used, or this script measures the
    robustness of a different analysis than the one being defended.

    The sys.modules registration is not optional (it is what tables/_h/02.genesets.py
    already does): `session_info.show()` walks the caller's globals and looks each
    referenced module up by name in sys.modules, so a module-level handle to an
    unregistered module raises KeyError on the script's last line -- after every
    output has been written. That is a job that reports FAILED while its results are
    complete, which is the worst of both.
    """
    name = "continuum_dpt"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


HEAD = load_headline(HEADLINE)
TREND_COLS = list(HEAD.TREND_COLS)


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad (from 00.state_discovery.py)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--use-rep", default="X_pca_harmony")
    p.add_argument("--canonical-root", default="vascular_stabilizing")
    p.add_argument("--denoise", type=Path, default=Path(HEAD.DENOISE_DEFAULT),
                   help="pericytes_airspace_denoising.tsv (scVI-denoised AGTR1); "
                        "carried so the sweep covers the same features as the "
                        "headline run")
    p.add_argument("--den-model", default="Pericyte-only-trained")
    return p.parse_args()


def ensure_agtr1(adata):
    if "AGTR1_expr" not in adata.obs and "AGTR1" in adata.var_names:
        expr = adata[:, "AGTR1"].layers["logcounts"]
        adata.obs["AGTR1_expr"] = (expr.toarray().ravel()
                                   if sparse.issparse(expr) else np.asarray(expr).ravel())


def resolve_root(adata, root_state, rep):
    """(root_idx, matched_on, root_cell_program) for one requested root.

    Delegates to the headline `pick_root` with the PC1 fallback OFF, so a root
    that matches nothing raises instead of quietly becoming the PC1 minimum.
    """
    if root_state in PC1_PSEUDOROOTS:
        idx = int(PC1_PSEUDOROOTS[root_state](adata.obsm[rep][:, 0]))
        matched_on = f"{root_state.upper()}_PSEUDOROOT"
    else:
        idx, matched_on = HEAD.pick_root(adata, root_state, rep, allow_fallback=False)
    program = (str(adata.obs["state_program"].iloc[idx])
               if "state_program" in adata.obs.columns else "NA")
    return idx, matched_on, program


def donor_rows(obs, feats):
    """Donor-aware rho: Spearman of per-donor mean pseudotime vs per-donor mean score.

    Aggregation mirrors 02.continuum_dpt.py::correlate_trends exactly (plain group
    means, drop donors with no pseudotime, require >= 5 donors) so that the canonical
    setting here reproduces the donor rows of pseudotime_trend_correlations.tsv.
    """
    if "donor_id" not in obs:
        return []
    agg = {"dpt_pseudotime": "mean", **{c: "mean" for c in feats}}
    donor = obs.groupby("donor_id", observed=True).agg(agg)
    donor = donor.dropna(subset=["dpt_pseudotime"])
    rows = []
    for col in feats:
        d = donor[["dpt_pseudotime", col]].dropna()
        if d.shape[0] < 5:
            continue
        rho, pval = spearmanr(d["dpt_pseudotime"], d[col])
        rows.append(dict(feature=col, level="donor", spearman_rho=rho,
                         p_value=pval, n=d.shape[0]))
    return rows


def one_run(adata, rep, root_state, neighbors, n_dcs, frac, seed):
    """Subsample (if frac<1), run neighbors->diffmap->dpt, return per-feature rho
    at both the cell and the donor level, each row tagged with the root actually
    used."""
    if frac < 1.0:
        rng = np.random.default_rng(seed)
        idx = rng.choice(adata.n_obs, size=int(round(frac * adata.n_obs)), replace=False)
        sub = adata[np.sort(idx)].copy()
    else:
        sub = adata.copy()
    root_idx, matched_on, root_program = resolve_root(sub, root_state, rep)
    logging.info("    root %r -> cell %d (matched_on=%s, root cell program=%s)",
                 root_state, root_idx, matched_on, root_program)
    sc.pp.neighbors(sub, use_rep=rep, n_neighbors=neighbors, random_state=seed)
    sc.tl.diffmap(sub, n_comps=max(n_dcs, 15))
    sub.uns["iroot"] = root_idx
    sc.tl.dpt(sub, n_dcs=n_dcs)
    pt = sub.obs["dpt_pseudotime"].to_numpy()
    feats = [c for c in TREND_COLS if c in sub.obs]
    rows = []
    for col in feats:
        rho, pval = spearmanr(pt, sub.obs[col].to_numpy(), nan_policy="omit")
        rows.append(dict(feature=col, level="cell", spearman_rho=rho,
                         p_value=pval, n=sub.n_obs))
    rows += donor_rows(sub.obs, feats)
    prov = dict(root_matched_on=matched_on, root_cell_program=root_program,
                root_cell_index=root_idx)
    for r in rows:
        r.update(prov)
    return rows


def alt_root_list(adata, canonical):
    """Roots for block (1), in the order they are reported.

    PROGRAMS first: rooting at the opposite biological pole is the invariance test
    the block exists for, and before the fix this script could not perform it at
    all -- its root list held only Leiden ids. Leiden CLUSTERS are kept because the
    published sweep used them. The PC1 extremes come last, and `pc1min` is
    deliberately retained: it is exactly the root the broken canonical silently
    used, so it stays in the table as the reference for what the old numbers were.
    """
    roots = []
    for col in ("state_program", "pericyte_state"):
        if col in adata.obs.columns:
            roots += [s for s in map(str, pd.unique(adata.obs[col]))
                      if s != canonical and s not in roots]
    return roots + list(PC1_PSEUDOROOTS)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.outdir / "figures"; fig_dir.mkdir(exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    ensure_agtr1(adata)
    HEAD.add_denoised_agtr1(adata, args.denoise, args.den_model)
    rep = args.use_rep

    # Gate: the canonical root must resolve to a cell of the canonical program.
    # This is the assertion whose absence let the sweep run at the wrong pole for
    # months while its own log said otherwise.
    can_idx, can_col, can_prog = resolve_root(adata, args.canonical_root, rep)
    if can_prog != args.canonical_root:
        raise ValueError(
            f"canonical root {args.canonical_root!r} resolved to cell {can_idx}, "
            f"whose state_program is {can_prog!r}. Refusing to sweep parameters "
            "around a root at the wrong pole -- that is the exact defect (P1-3) "
            "this gate exists to catch.")
    logging.info("Canonical root OK: %r matched on obs[%r], root cell %d, program %s",
                 args.canonical_root, can_col, can_idx, can_prog)

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
    for root_state in alt_root_list(adata, args.canonical_root):
        logging.info(f"[root] root={root_state} (default params)")
        for r in one_run(adata, rep, root_state, 30, 10, 1.0, 13):
            r.update(block="root", root=str(root_state),
                     neighbors=30, n_dcs=10, frac=1.0, seed=13)
            runs.append(r)

    df = pd.DataFrame(runs)
    df.to_csv(args.outdir / "continuum_sensitivity_runs.tsv", sep="\t", index=False)

    # ---- root audit: what each requested root actually resolved to ---------
    audit = (df[["block", "root", "root_matched_on", "root_cell_program",
                 "root_cell_index"]]
             .drop_duplicates()
             .sort_values(["block", "root"]))
    audit["is_canonical"] = audit["root"] == args.canonical_root
    audit.to_csv(args.outdir / "continuum_sensitivity_root_audit.tsv",
                 sep="\t", index=False)
    print("\n== Root audit (requested -> cell actually used) ==")
    print(audit.to_string(index=False))

    # ---- summary over canonical-root parameter sweep -----------------------
    can = df[df["block"] == "param"]
    summ = (can.groupby(["feature", "level"])["spearman_rho"]
            .agg(n_settings="count", mean_rho="mean", sd_rho="std",
                 min_rho="min", max_rho="max")
            .reset_index())
    # fraction of canonical-root settings with the same sign as the mean
    sign = (can.assign(sgn=np.sign(can["spearman_rho"]))
            .groupby(["feature", "level"])["sgn"]
            .apply(lambda s: (s == np.sign(s.mean())).mean())
            .rename("sign_consistency").reset_index())
    summ = summ.merge(sign, on=["feature", "level"])
    # Carry the root's provenance into the summary itself: this table is what gets
    # quoted, and quoting it without its root is how the old numbers travelled.
    summ["root"] = args.canonical_root
    summ["root_matched_on"] = can_col
    summ["root_cell_program"] = can_prog
    summ.to_csv(args.outdir / "continuum_sensitivity_summary.tsv", sep="\t", index=False)
    print("\n== Canonical-root parameter robustness (Spearman rho vs DPT) ==")
    print(summ.to_string(index=False))

    # ---- figure: rho across all settings, faceted by block -----------------
    order = [c for c in TREND_COLS if c in set(df["feature"])]
    g = sns.catplot(data=df, x="spearman_rho", y="feature", hue="block",
                    col="level", col_order=["cell", "donor"],
                    order=order, kind="strip", dodge=True, height=4.5, aspect=1.2,
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
