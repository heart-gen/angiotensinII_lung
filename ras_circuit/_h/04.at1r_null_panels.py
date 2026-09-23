"""
Detection-matched null for the AT1R-response signature (Figures 5D, 5E, S18).

A panel score regressed against another panel score has a strongly positive null
(+0.45 against the BM score; memory: gene-set-score-null-not-zero), and even the
count-model arbiter must be asked whether the observed signature is unusual among
panels of the same sparsity. This builds that yardstick with the functions of
basement_membrane/_h/11.tgfb_null_panels.py, unchanged: each null panel matches the
kept signature gene-by-gene on pericyte detection, is scored with the IDENTICAL
sc.tl.score_genes call (up minus down, as in 03), and is aggregated to the same
donor x pericyte_state units the models use.

The candidate pool excludes every panel gene 03 prunes against, the signature
itself, and the legacy AT1R_PROGRAM, so no null panel can carry the association
under test by construction.
"""
import argparse
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "basement_membrane" / "_h"))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", type=Path, required=True)
    p.add_argument("--audit", type=Path, required=True,
                   help="stats_data/at1r_signature_audit.tsv from 03")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--n-null", type=int, default=1000)
    p.add_argument("--seed", type=int, default=13)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    nm = load(ROOT / "basement_membrane" / "_h" / "11.tgfb_null_panels.py", "tgfb_null")
    sc3 = load(Path(__file__).resolve().parent / "03.at1r_response_score.py", "at1r03")
    rng = np.random.default_rng(a.seed)

    adata = nm.load_anndata(a.adata)
    det = nm.detection_frac(adata)
    aud = pd.read_csv(a.audit, sep="\t")
    up = aud.query("kept and direction == 'up'")["gene"].tolist()
    dn = aud.query("kept and direction == 'down'")["gene"].tolist()
    use_dn = len(dn) >= 5

    ex = set(sc3.excluded()) | set(aud["gene"]) | set(sc3.AT1R_PROGRAM) | nm.excluded_genes()
    pool = [g for g in adata.var_names if g not in ex and det[g] > 0]
    logging.info(f"null pool: {len(pool)} genes; matching {len(up)} up"
                 + (f" + {len(dn)} down" if use_dn else ""))

    up_null = nm.matched_panels(up, det, pool, a.n_null, rng)
    dn_null = nm.matched_panels(dn, det, pool, a.n_null, rng) if use_dn else [None] * a.n_null
    logging.info(f"target detection up {det[up].mean():.4f}, achieved "
                 f"{np.mean([det[p].mean() for p in up_null]):.4f}")

    store, rows = {}, []
    for k in range(a.n_null):
        col = f"null_at1r_{k + 1:04d}"
        nm.score_into(adata, up_null[k], "__u", a.seed, store)
        s = store.pop("__u")
        if use_dn:
            nm.score_into(adata, dn_null[k], "__d", a.seed, store)
            s = s - store.pop("__d")
        store[col] = s
        rows.append({"panel": col, "up_genes": ",".join(up_null[k]),
                     "down_genes": ",".join(dn_null[k]) if use_dn else "",
                     "mean_detect_up": float(det[up_null[k]].mean())})
        if (k + 1) % 100 == 0:
            logging.info(f"scored {k + 1}/{a.n_null}")
    pd.DataFrame(rows).to_csv(a.outdir / "at1r_null_panels.tsv", sep="\t", index=False)

    cols = list(store)
    df = pd.concat([adata.obs[["donor_id", "pericyte_state"]].reset_index(drop=True),
                    pd.DataFrame(store)], axis=1)
    pbn = df.groupby(["donor_id", "pericyte_state"], observed=True)[cols].mean().reset_index()
    pbn.to_csv(a.outdir / "at1r_null_pseudobulk.tsv.gz", sep="\t", index=False)
    # Cell-level null scores are needed for the within-donor Spearman null in 06.
    cell = pd.DataFrame(store, index=adata.obs_names)[cols[: min(200, len(cols))]]
    cell.index.name = "index"
    cell.to_csv(a.outdir / "at1r_null_cells.tsv.gz", sep="\t")
    logging.info(f"null pseudobulk: {pbn.shape[0]} units x {len(cols)} panels")


if __name__ == "__main__":
    main()
