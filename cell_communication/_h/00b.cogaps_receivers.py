"""
Add an ORTHOGONAL CoGAPS-derived pericyte receiver scheme to the CCC niche.

The primary CCC state scheme (00.prepare_ccc_input.py -> ccc_group_state) splits
pericytes by the curated NVU-pattern programs (state_program). This module adds a
fully data-driven alternative: split pericytes by their dominant CoGAPS pattern
(unsupervised Bayesian NMF, pericyte_cogaps/), so the headline conclusion -- which
niche ligands target the injury-leaning pericytes -- can be shown to NOT depend on
how the pericyte states were defined.

Pipeline:
  - Read the already-built ccc_niche.h5ad (cells capped per donor; pericytes carry
    the same barcodes as the CoGAPS run, both descend from pericyte_states).
  - Read patterns_npN.tsv.gz (cell x pattern weights). Dominant pattern per cell =
    argmax of the non-negative pattern weights (standard CoGAPS assignment).
  - Annotate each pattern with its best-matching curated program from the validation
    cell-level Spearman table (pattern_score_spearman.tsv.gz): the program with the
    highest positive rho (>= --min-rho), else "mixed". This makes the CoGAPS receivers
    directly comparable to the curated-state receivers.
  - Build ccc_group_cogaps: pericytes -> "Pericyte_cg_<program>" (patterns sharing a
    program collapse); unmatched pericytes -> "Pericyte_cg_unassigned"; every other
    cell keeps its ccc_group (main) label, so senders / AT2 receivers are unchanged.
  - Write the annotation table and overwrite ccc_niche.h5ad with the new column.

Outputs:
  - cogaps_receiver_annotation_npN.tsv  (pattern -> program, rho per program, n cells)
  - ccc_group_cogaps_counts.tsv         (receiver x disease)
  - ccc_niche.h5ad                       (with ccc_group_cogaps added)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--niche", required=True, type=Path, help="ccc_niche.h5ad")
    p.add_argument("--patterns", required=True, type=Path,
                   help="patterns_npN.tsv.gz (barcode x Pattern_*)")
    p.add_argument("--score-spearman", required=True, type=Path,
                   help="validation_npN/pattern_score_spearman.tsv.gz")
    p.add_argument("--npatterns", type=int, required=True)
    p.add_argument("--min-rho", type=float, default=0.15,
                   help="Min positive Spearman rho to name a pattern by a program")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def annotate_patterns(score_path, pattern_names, min_rho):
    """pattern -> best curated program (max positive rho >= min_rho), else 'mixed'."""
    spear = pd.read_csv(score_path, sep="\t")
    spear = spear.dropna(subset=["rho"])
    # ignore non-program score columns (e.g. anatomical_region_ccf_score is all-NA)
    spear = spear[spear["score"].str.endswith("_score")]
    spear["program"] = spear["score"].str.replace("_score$", "", regex=True)
    ann = {}
    rho_tbl = {}
    for pat in pattern_names:
        sub = spear[spear["pattern"] == pat]
        rho_tbl[pat] = sub.set_index("program")["rho"].to_dict()
        if sub.empty:
            ann[pat] = "mixed"
            continue
        best = sub.loc[sub["rho"].idxmax()]
        ann[pat] = best["program"] if best["rho"] >= min_rho else "mixed"
    return ann, rho_tbl


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    pat = pd.read_csv(args.patterns, sep="\t", index_col=0)
    pattern_names = list(pat.columns)
    logging.info(f"CoGAPS patterns: {len(pat)} cells x {len(pattern_names)} patterns")

    # Dominant pattern per cell = argmax of non-negative pattern weights
    dom = pat.to_numpy().argmax(axis=1)
    dom_pat = pd.Series([pattern_names[i] for i in dom], index=pat.index, name="dom_pattern")

    ann, rho_tbl = annotate_patterns(args.score_spearman, pattern_names, args.min_rho)
    logging.info(f"Pattern -> program: {ann}")

    # barcode -> Pericyte_cg_<program>
    cg_label = dom_pat.map(lambda p: f"Pericyte_cg_{ann[p]}")

    # Annotation table (auditable pattern -> program with per-program rho + cell counts)
    dom_counts = dom_pat.value_counts().to_dict()
    rows = []
    for p in pattern_names:
        row = {"pattern": p, "assigned_program": ann[p],
               "n_dominant_cells": int(dom_counts.get(p, 0))}
        row.update({f"rho_{k}": v for k, v in rho_tbl.get(p, {}).items()})
        rows.append(row)
    pd.DataFrame(rows).to_csv(
        args.outdir / f"cogaps_receiver_annotation_np{args.npatterns}.tsv",
        sep="\t", index=False)

    adata = sc.read_h5ad(args.niche)
    base = adata.obs["ccc_group"].astype(str).to_numpy().copy()
    is_peri = (adata.obs["cell_type"].astype(str) == "Pericytes").to_numpy()

    mapped = pd.Series(adata.obs_names, index=adata.obs_names).map(cg_label)
    n_match = int((is_peri & mapped.notna().to_numpy()).sum())
    logging.info(f"CoGAPS label matched {n_match}/{is_peri.sum()} niche pericytes "
                 f"({n_match / max(is_peri.sum(), 1):.1%})")

    cg = np.where(mapped.isna().to_numpy(), "Pericyte_cg_unassigned",
                  mapped.astype(str).to_numpy())
    g = base.copy()
    g[is_peri] = cg[is_peri]
    adata.obs["ccc_group_cogaps"] = pd.Categorical(g)
    adata.obs["cogaps_dominant_pattern"] = pd.Series(adata.obs_names, index=adata.obs_names)\
        .map(dom_pat).to_numpy()

    pd.crosstab(adata.obs["ccc_group_cogaps"], adata.obs["disease_group"]).to_csv(
        args.outdir / "ccc_group_cogaps_counts.tsv", sep="\t")

    adata.write(args.niche)
    logging.info(f"Wrote ccc_group_cogaps into {args.niche.name}; "
                 f"receivers: {sorted(set(g[is_peri]))}")
    session_info.show()


if __name__ == "__main__":
    main()
