"""
Detection-matched empirical null for the AGTR1-versus-matrix models.

WHY.  The TGF-beta specificity run established that a gene-set score regressed on
the BM score has a strongly POSITIVE baseline: detection-matched random panels
give beta = +0.45 against basement_membrane_score_z, because any score shares a
general-expression component with the BM score that the depth covariate does not
absorb. Against BM - fibrillar the same null sits at +0.026.

`04.bm_state_stats.R` runs the identical model with AGTR1 as the predictor, and
reports AGTR1 -> BM as significantly positive (+0.170 expr, +0.209 detect). If
the single-gene null is also centred above zero, those estimates are being tested
against the wrong reference and the sign of the CLAIM, not just its size, is in
question.

AGTR1 is a single gene, not a panel, so the panel result does not transfer -- it
has to be measured. This builds the matching null: K random single genes matched
on pericyte detection rate to AGTR1, each supplying both an expression and a
detection predictor, aggregated to the donor x cluster units the models use.

The denoised lens (AGTR1_scvi) gets no null here: reproducing it would require
retraining scVI per null gene. Stated as a limitation rather than approximated.

Outputs (to --outdir):
  - agtr1_null_pseudobulk.tsv.gz  donor x cluster means of every null gene
  - agtr1_null_genes.tsv          which gene each draw used, for audit
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path

import bm_panels


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--gene", default="AGTR1")
    p.add_argument("--tol", type=float, default=0.03,
                   help="max |detection - target| for a gene to enter the null")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def main():
    configure_logging()
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    adata = sc.read_h5ad(args.adata)
    if "feature_name" in adata.var.columns:
        sym = adata.var["feature_name"].astype(str)
        if not (sym.values == adata.var_names.values).all():
            adata.var_names = sym
            adata.var_names_make_unique()
    X = adata.layers["logcounts"]
    dense = lambda m: np.asarray(m.todense()) if hasattr(m, "todense") else np.asarray(m)

    n_pos = (np.asarray((X > 0).sum(axis=0)).ravel() if hasattr(X, "getnnz")
             else (np.asarray(X) > 0).sum(axis=0))
    detect = pd.Series(n_pos / adata.n_obs, index=adata.var_names)

    if args.gene not in adata.var_names:
        raise KeyError(f"{args.gene} absent from the object")
    target = float(detect[args.gene])
    logging.info(f"{args.gene} pericyte detection: {target:.4f}")

    # Exclude every panel gene AND the target itself: a null gene must not be a
    # matrix gene, or the association being tested is built into the null.
    ex = {args.gene}
    for genes in bm_panels.PANELS.values():
        ex |= set(genes)
    for _, genes in bm_panels.GENE_BLOCKS:
        ex |= set(genes)
    pool = np.asarray([g for g in adata.var_names
                       if g not in ex and detect[g] > 0])
    pool_det = detect[pool].to_numpy()
    logging.info(f"pool: {len(pool)} genes ({len(ex)} excluded)")

    # EVERY gene inside a detection tolerance, not the K nearest.
    #
    # A single gene cannot be matched the way a 7-gene panel can: only 110 genes
    # sit within +/-0.02 of AGTR1's 0.374 detection, so asking for 1,000
    # distinct neighbours reaches far down the tail and drags the null's mean
    # detection to 0.326 -- a 13% mismatch on the one variable that has to
    # match. Taking the whole window instead keeps the match tight and lets n
    # be whatever the data supports, which is the honest constraint.
    keep = np.abs(pool_det - target) <= args.tol
    picked = pool[keep]
    if picked.size < 50:
        raise RuntimeError(
            f"only {picked.size} genes within +/-{args.tol} of {args.gene}'s "
            f"detection ({target:.4f}); widen --tol rather than proceeding on "
            "a null too small to place an estimate against.")
    got = detect[picked].mean()
    logging.info(f"{picked.size} null genes within +/-{args.tol} of target; "
                 f"detection target {target:.4f}, achieved mean {got:.4f} "
                 f"[{detect[picked].min():.4f}, {detect[picked].max():.4f}]")
    pd.DataFrame({"gene": picked,
                  "detect_frac": detect[picked].to_numpy(),
                  "target_gene": args.gene, "target_detect": target}
                 ).to_csv(args.outdir / "agtr1_null_genes.tsv",
                          sep="\t", index=False)

    ix = [adata.var_names.get_loc(g) for g in picked]
    mat = dense(X[:, ix])

    obs = adata.obs
    for need in ("donor_id", "pericyte_state"):
        if need not in obs.columns:
            raise KeyError(f"obs lacks '{need}'")

    # Both lenses the R models use: mean log-expression, and detection rate.
    cols = {}
    for j, g in enumerate(picked):
        cols[f"nullexpr_{j:04d}"] = mat[:, j].astype(np.float32)
        cols[f"nulldet_{j:04d}"] = (mat[:, j] > 0).astype(np.float32)
    df = pd.concat([obs[["donor_id", "pericyte_state"]].reset_index(drop=True),
                    pd.DataFrame(cols)], axis=1)
    pbn = df.groupby(["donor_id", "pericyte_state"], observed=True).mean().reset_index()
    logging.info(f"null pseudobulk: {pbn.shape[0]} units x {len(cols)} predictors")
    pbn.to_csv(args.outdir / "agtr1_null_pseudobulk.tsv.gz", sep="\t", index=False)
    session_info.show()


if __name__ == "__main__":
    main()
