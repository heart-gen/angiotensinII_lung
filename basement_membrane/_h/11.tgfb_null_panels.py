"""
Detection-matched empirical null for the TGF-beta specificity arms.

WHY THIS EXISTS.  The specificity test splits TGFB_RESPONSE into a SMAD-proximal
arm (mean pericyte detection 9.2%) and an immediate-early/mechano arm (32.2%).
The arms are confounded with detection *by construction* -- that is the concern
restated, not a flaw in the split -- so a null result on the sparse SMAD arm is
uninterpretable on its own: we would not know whether the arm has no signal or
simply cannot carry one at that sparsity.

This script builds the missing yardstick.  For each arm it draws K random gene
panels matched gene-by-gene on pericyte detection rate, scores each with the
IDENTICAL sc.tl.score_genes call, and aggregates to the same donor x cluster
units the real models use.  12.tgfb_specificity.R then fits the same model to
every null panel, which yields (a) an empirical p for each observed arm against
panels of equal sparsity, and (b) a power statement -- what fraction of matched
panels could have reached the reported |beta| = 0.196 at all.

The design and its decision rule were fixed in _h/TGFB_SPECIFICITY_PLAN.md
before any of this was fitted.

Outputs (to --outdir):
  - tgfb_null_pseudobulk.tsv.gz  donor x cluster mean of every null panel score
  - tgfb_null_panels.tsv         which genes each null panel drew, for audit
  - tgfb_arm_detection.tsv       per-gene detection for the two real arms
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path

import bm_panels


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad (same input as 00.bm_score.py)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--n-null", type=int, default=1000,
                   help="null panels PER ARM")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def load_anndata(path: Path):
    adata = sc.read_h5ad(path)
    if "logcounts" not in adata.layers:
        raise KeyError("expected a 'logcounts' layer")
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        if not (symbols.values == adata.var_names.values).all():
            adata.var_names = symbols
            adata.var_names_make_unique()
    adata.X = adata.layers["logcounts"]
    return adata


def detection_frac(adata):
    """Per-gene fraction of pericytes with non-zero logcounts."""
    X = adata.layers["logcounts"]
    if hasattr(X, "getnnz"):
        n_pos = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        n_pos = (np.asarray(X) > 0).sum(axis=0)
    return pd.Series(n_pos / adata.n_obs, index=adata.var_names)


def excluded_genes():
    """Every gene that appears in any panel this analysis regresses against.

    A null panel must not accidentally draw a matrix or state gene: that would
    build the association we are testing for into the null itself, and the
    empirical p would be conservative for the wrong reason.
    """
    ex = set()
    for genes in bm_panels.PANELS.values():
        ex |= set(genes)
    for _, genes in bm_panels.GENE_BLOCKS:
        ex |= set(genes)
    return ex


def matched_panels(arm_genes, detect, pool, n_null, rng):
    """n_null panels, each matching arm_genes gene-by-gene on detection rate.

    Nearest neighbour on detect_frac, sampled without replacement WITHIN a
    panel so a panel never doubles a gene. Ties are broken at random rather
    than by gene order, so the null is not systematically drawn from one end
    of the alphabet.
    """
    pool_det = detect[pool].to_numpy()
    pool_arr = np.asarray(pool)
    targets = detect[arm_genes].to_numpy()
    panels = []
    for _ in range(n_null):
        used = set()
        picked = []
        for t in targets:
            d = np.abs(pool_det - t)
            # jitter breaks ties randomly; 40 nearest then a random draw keeps
            # the match tight without making every panel the same panel.
            order = np.argsort(d + rng.uniform(0, 1e-12, d.size))
            for ix in order:
                g = pool_arr[ix]
                if g not in used:
                    used.add(g)
                    picked.append(g)
                    break
        panels.append(picked)
    return panels


def main():
    configure_logging()
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    adata = load_anndata(args.adata)
    logging.info(f"pericytes: {adata.n_obs} cells x {adata.n_vars} genes")

    detect = detection_frac(adata)

    arms = {"smad": bm_panels.TGFB_SMAD, "ieg": bm_panels.TGFB_IEG}
    det_rows = []
    for arm, genes in arms.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            raise KeyError(f"[{arm}] genes absent from the object: {missing}. "
                           "The specificity arms are pre-specified; a missing "
                           "gene changes the panel and must be handled "
                           "explicitly, not dropped silently.")
        for g in genes:
            det_rows.append({"arm": arm, "gene": g,
                             "detect_frac": float(detect[g])})
        logging.info(f"[{arm}] {len(genes)} genes, mean detection "
                     f"{detect[genes].mean():.3f}")
    pd.DataFrame(det_rows).to_csv(
        args.outdir / "tgfb_arm_detection.tsv", sep="\t", index=False)

    # Candidate pool: expressed at all, and in no panel we test against.
    ex = excluded_genes()
    pool = [g for g in adata.var_names if g not in ex and detect[g] > 0]
    logging.info(f"null candidate pool: {len(pool)} genes "
                 f"({len(ex)} panel genes excluded)")

    panel_rows, score_cols = [], {}
    for arm, genes in arms.items():
        nulls = matched_panels(genes, detect, pool, args.n_null, rng)
        # Report how tight the match actually is -- if it is loose, the null
        # is not answering the question and the reader must be told.
        tgt = detect[genes].mean()
        got = np.mean([detect[p].mean() for p in nulls])
        logging.info(f"[{arm}] {args.n_null} null panels; mean detection "
                     f"target {tgt:.4f}, achieved {got:.4f}")
        for k, panel in enumerate(nulls, start=1):
            col = f"null_{arm}_{k:04d}"
            panel_rows.append({"arm": arm, "panel": col,
                               "genes": ",".join(panel),
                               "mean_detect": float(detect[panel].mean())})
            sc.tl.score_genes(adata, panel, score_name=col,
                              random_state=args.seed, use_raw=False)
            score_cols[col] = arm
            if k % 100 == 0:
                logging.info(f"[{arm}] scored {k}/{args.n_null}")
    pd.DataFrame(panel_rows).to_csv(
        args.outdir / "tgfb_null_panels.tsv", sep="\t", index=False)

    # Leave-one-gene-out on the full 17-gene panel. A panel score can be carried
    # by one dominant gene -- JUNB is detected in 56.5% of pericytes, more than
    # the entire SMAD arm combined -- and re-scoring is the only honest way to
    # test that, because sc.tl.score_genes subtracts an expression-matched
    # control set and is not a plain mean over genes.
    full = [g for g in bm_panels.TGFB_RESPONSE if g in adata.var_names]
    for g in full:
        col = f"logo_{g}"
        sc.tl.score_genes(adata, [x for x in full if x != g], score_name=col,
                          random_state=args.seed, use_raw=False)
        score_cols[col] = "logo"
    logging.info(f"leave-one-out: {len(full)} refits of the full panel")

    # Aggregate to the donor x cluster units the real models use. The n_cells
    # floor is applied in R against its own pseudobulk, so every unit is
    # emitted here and the merge decides.
    obs = adata.obs
    for need in ("donor_id", "pericyte_state"):
        if need not in obs.columns:
            raise KeyError(f"obs is missing '{need}'")
    cols = list(score_cols)
    df = pd.DataFrame(obs[["donor_id", "pericyte_state"]]).copy()
    df[cols] = adata.obs[cols].to_numpy()
    pbn = df.groupby(["donor_id", "pericyte_state"], observed=True)[cols].mean()
    pbn = pbn.reset_index()
    logging.info(f"null pseudobulk: {pbn.shape[0]} units x {len(cols)} panels")
    pbn.to_csv(args.outdir / "tgfb_null_pseudobulk.tsv.gz",
               sep="\t", index=False)

    session_info.show()


if __name__ == "__main__":
    main()
