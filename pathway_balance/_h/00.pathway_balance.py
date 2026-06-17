"""
AGTR1 (AT1R) vs AGTR2 (AT2R) downstream pathway-balance scoring.

Strengthens the therapeutic rationale for AT1 blockade (losartan) by quantifying
the balance between AT1R-driven (pro-fibrotic / proliferative / contractile /
inflammatory) and AT2R-driven (anti-proliferative / vasodilatory / apoptotic)
transcriptional programs per cell, then per donor and per pericyte state.

    AT1R_AT2R_balance = AT1R_program_score - AT2R_program_score

A positive, disease-amplified balance provides a transcriptomic rationale for AT1R
blockade. IMPORTANT INTERPRETIVE CAVEAT: the AT1R downstream program is built from
generic pro-fibrotic/inflammatory effector genes (TGFB1, CCN2, SERPINE1, COL1A1/3A1,
FN1, ACTA2, IL6, CCL2, ...) that OVERLAP heavily with the injury-stromal program and
the NicheNet target gene set. The balance is therefore a profibrotic-vs-protective
transcriptional axis that partly RE-MEASURES injury intensity; it is NOT a
receptor-specific readout of AT1R signaling (AGTR1 itself is not disease-associated;
see pericyte_states/03.agtr1_lenses). The disease shift is also driven at the
disease level, not cleanly by pericyte program (program contrasts are NS; smallest
pairwise p = 0.067). Read it as consistent with, not proof of, AT1R dominance.

Inputs: pericyte_states.h5ad (Module 2). Memory-light (small object).
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse

# AT1R downstream: Gq/PLC -> pro-fibrotic, proliferative, contractile, inflammatory
AT1R_PROGRAM = [
    "TGFB1", "CCN2", "CTGF", "CCN1", "CYR61", "SERPINE1", "COL1A1", "COL3A1",
    "FN1", "ACTA2", "EDN1", "NOX4", "IL6", "CCL2", "NFKB1", "MYC", "EGR1",
    "FOS", "AGT", "ACE", "PLAU",
]
# AT2R downstream: anti-proliferative, vasodilatory, apoptotic, phosphatase
AT2R_PROGRAM = [
    "NOS3", "NOS1", "BDKRB2", "PTPN6", "DUSP1", "DUSP6", "CASP3", "BAX",
    "PPARG", "NPR1", "NPPB", "MAPK1",
]
RECEPTORS = ["AGTR1", "AGTR2"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def add_expr(adata, gene):
    if gene not in adata.var_names:
        return
    x = adata[:, gene].layers["logcounts"]
    x = x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
    adata.obs[f"{gene}_expr"] = x


def score(adata, genes, name, seed):
    present = [g for g in genes if g in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if missing:
        logging.warning(f"[{name}] absent: {missing}")
    sc.tl.score_genes(adata, present, score_name=name, random_state=seed, use_raw=False)
    return present


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    adata.X = adata.layers["logcounts"]

    at1 = score(adata, AT1R_PROGRAM, "AT1R_score", args.seed)
    at2 = score(adata, AT2R_PROGRAM, "AT2R_score", args.seed)
    adata.obs["AT1R_AT2R_balance"] = adata.obs["AT1R_score"] - adata.obs["AT2R_score"]
    for g in RECEPTORS:
        add_expr(adata, g)

    # Record which genes were used
    with open(args.outdir / "balance_signature_genes.txt", "w") as fh:
        fh.write("AT1R_program: " + ", ".join(at1) + "\n")
        fh.write("AT2R_program: " + ", ".join(at2) + "\n")

    keep = [c for c in [
        "donor_id", "disease", "lung_condition", "smoking_status", "sex",
        "self_reported_ethnicity", "age_or_mean_of_age_range", "leiden_pericytes",
        "dataset", "study", "pericyte_state", "state_program",
        "AT1R_score", "AT2R_score",
        "AT1R_AT2R_balance", "AGTR1_expr", "AGTR2_expr"] if c in adata.obs.columns]
    adata.obs[keep].to_csv(args.outdir / "pathway_balance_metadata.tsv.gz", sep="\t")
    logging.info(f"Scored balance for {adata.n_obs} cells")
    session_info.show()


if __name__ == "__main__":
    main()
