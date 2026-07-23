"""
Donor x cell-type pseudobulk over the CCC niche, for an arbitrary gene panel.

Generic on purpose: `basement_membrane/step_2.sh` calls it with the BM +
fibrillar-contrast panel, and `agt_axis/step_0.sh` calls it with the local-RAS
panel. One implementation, two gene lists -- the alternative was duplicating the
aggregation logic across modules, which is how the marker panels in this repo
ended up copy-pasted across six files.

NORMALIZATION (verified, load-bearing): ccc_niche.h5ad `X` is log1p(CP10K) --
expm1(X) sums to exactly 10,000 per cell. Pseudobulk therefore has to be
expm1 -> mean -> log1p. Averaging log values instead would introduce a Jensen
bias that scales with sequencing depth, which is precisely the cross-cell-type
confound this analysis is trying to control for. The assertion below fails loudly
if that ever stops being true.

Raw `counts` / `soupX` layers exist but are populated for only ~5% of cells
(Meyer_2021 only), so no count-based model is possible atlas-wide; everything
downstream works from the normalized values.

Output: one row per (donor_id x cell_type group), with per-gene
log1p(mean CP10K) and detection fraction, plus the covariates the LMMs need.
"""
import numpy as np
import pandas as pd
import anndata as ad
import session_info
import logging, argparse
from pathlib import Path


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="ccc_niche.h5ad")
    p.add_argument("--genes", required=True, type=Path,
                   help="TSV with a `gene` column (extra columns ignored)")
    p.add_argument("--outfile", required=True, type=Path)
    p.add_argument("--group-key", default="ccc_group")
    p.add_argument("--chunk", type=int, default=20000)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def check_normalization(adata, seed, n=300, tol=0.01):
    """expm1(X) must sum to 1e4 per cell; the pseudobulk recipe depends on it."""
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(adata.n_obs, min(n, adata.n_obs), replace=False))
    X = adata[idx].to_memory().X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    sums = np.expm1(X).sum(axis=1)
    med = float(np.median(sums))
    if abs(med - 1e4) / 1e4 > tol:
        raise ValueError(
            f"X is not log1p(CP10K): median expm1 row-sum = {med:.1f}, "
            "expected 10000. The expm1->mean->log1p pseudobulk recipe is invalid "
            "for this object; do not proceed.")
    logging.info(f"Normalization check passed (median expm1 row-sum {med:.2f})")


def main():
    configure_logging()
    args = parse_args()
    args.outfile.parent.mkdir(parents=True, exist_ok=True)

    panel = pd.read_csv(args.genes, sep="\t")
    want = sorted(set(panel["gene"].astype(str)))

    adata = ad.read_h5ad(args.adata, backed="r")
    logging.info(f"Opened {args.adata} ({adata.n_obs} x {adata.n_vars})")
    check_normalization(adata, args.seed)

    present = [g for g in want if g in adata.var_names]
    missing = sorted(set(want) - set(present))
    if missing:
        logging.warning(f"absent from object: {missing}")
    if not present:
        raise ValueError("no panel genes present in the object")
    logging.info(f"Aggregating {len(present)}/{len(want)} genes")

    obs_keys = [k for k in [args.group_key, "donor_id", "study", "dataset",
                            "disease_group", "assay", "suspension_type",
                            "log10_total_counts", "total_counts"]
                if k in adata.obs.columns]
    obs = adata.obs[obs_keys].copy()
    obs["_row"] = np.arange(adata.n_obs)

    gcols = [adata.var_names.get_loc(g) for g in present]
    # Accumulate sums per (donor x group) rather than holding the matrix:
    # 423k x 55k backed, sliced to the panel columns, stays under ~1 GB.
    unit = (obs[args.group_key].astype(str) + "||" + obs["donor_id"].astype(str))
    codes, levels = pd.factorize(unit)
    n_units = len(levels)
    sum_cp10k = np.zeros((n_units, len(present)), dtype=np.float64)
    sum_det = np.zeros((n_units, len(present)), dtype=np.float64)
    n_cells = np.zeros(n_units, dtype=np.int64)

    for start in range(0, adata.n_obs, args.chunk):
        stop = min(start + args.chunk, adata.n_obs)
        block = adata[start:stop].to_memory()
        Xb = block.X[:, gcols]
        Xb = Xb.toarray() if hasattr(Xb, "toarray") else np.asarray(Xb)
        cp10k = np.expm1(Xb)
        cb = codes[start:stop]
        np.add.at(sum_cp10k, cb, cp10k)
        np.add.at(sum_det, cb, (Xb > 0).astype(np.float64))
        np.add.at(n_cells, cb, 1)
        logging.info(f"  aggregated {stop}/{adata.n_obs} cells")

    mean_cp10k = sum_cp10k / n_cells[:, None]
    out = pd.DataFrame(np.log1p(mean_cp10k),
                       columns=[f"{g}__expr" for g in present])
    for j, g in enumerate(present):
        out[f"{g}__detect"] = sum_det[:, j] / n_cells

    meta = pd.DataFrame({
        args.group_key: [s.split("||")[0] for s in levels],
        "donor_id": [s.split("||")[1] for s in levels],
        "n_cells": n_cells,
    })
    # Unit-level covariates: constant within donor, averaged within unit for depth.
    agg = obs.groupby(unit.values, observed=True).agg(
        study=("study", "first") if "study" in obs.columns else ("donor_id", "first"),
        dataset=("dataset", "first") if "dataset" in obs.columns else ("donor_id", "first"),
        disease_group=("disease_group", "first") if "disease_group" in obs.columns else ("donor_id", "first"),
        mean_log10_total_counts=("log10_total_counts", "mean"),
    ).reindex(levels)
    meta = pd.concat([meta, agg.reset_index(drop=True)], axis=1)

    res = pd.concat([meta, out], axis=1)
    res.to_csv(args.outfile, sep="\t", index=False)
    logging.info(f"Wrote {args.outfile}: {res.shape[0]} units x {res.shape[1]} cols")
    logging.info("Units per group (n_cells >= 5):\n" +
                 res[res.n_cells >= 5][args.group_key].value_counts().to_string())

    session_info.show()


if __name__ == "__main__":
    main()
