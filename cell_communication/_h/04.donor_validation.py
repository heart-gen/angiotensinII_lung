"""
Donor-level validation of the NicheNet prediction.

NicheNet prioritizes ligands by a network ranking (not a donor-replicated test).
This script extracts, per donor, (i) the expression of the top prioritized
ligands in their niche sender cells and (ii) the pericyte injury/activation
target-program expression in ALL their pericytes, so the prediction
"fibroblast-derived TGF-beta drives the pericyte injury program" can be tested at
the donor level (05.donor_validation_stats.R). Using all pericytes (not an AGTR1+
subset) keeps the test independent of binarizing a dropout-prone GPCR transcript.

Memory-light: reads only the required gene columns from ccc_niche.h5ad.
"""
import numpy as np
import pandas as pd
import anndata as ad
import logging, argparse
from pathlib import Path
from scipy import sparse

# Top NicheNet-prioritized ligands into pericytes (TGF-beta family first)
LIGANDS = ["TGFB1", "TGFB2", "TGFB3", "CCN2", "TIMP2", "ICAM1"]
# Pericyte injury / activation target program (matches the NicheNet geneset)
TARGET_PROGRAM = ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "SPARC",
                  "MMP2", "TIMP1", "POSTN", "ACTA2", "TAGLN", "ADAMTS1",
                  "SERPINE1", "IL6", "CXCL2", "ICAM1", "VCAM1"]
# TGF-beta PATHWAY ACTIVITY in the RECEIVER, added 2026-09-07 (defect P1-6).
#
# WHY. The ligand-specific claim does not validate: TGFB2 ranks first by NicheNet
# AUPR and has essentially no donor-level association with the receiver program
# (rho 0.064, p 0.52), while TGFB1 sits on the significance boundary and moves
# across it depending on whether dataset is random or fixed. A ligand TRANSCRIPT
# in the sender is a poor proxy for signalling anyway -- TGF-beta is secreted
# latent and activated post-translationally, so sender mRNA need not track
# receiver exposure. The defensible version of the question is whether TGF-beta
# PATHWAY ACTIVITY IN THE PERICYTE tracks the injury program.
#
# The two arms come from basement_membrane/_h/bm_panels.py and are kept split for
# the reason established there (TGFB_SPECIFICITY_PLAN.md, pre-registered): the
# full TGFB_RESPONSE panel is dominated by an immediate-early/dissociation
# program, so a result on the whole panel cannot be attributed to SMAD signalling.
#   TGFB_SMAD -> canonical SMAD2/3 negative-feedback targets. THE READOUT.
#   TGFB_IEG  -> AP-1 / BMP / YAP-TAZ genes that respond to warm dissociation.
#                A DISCRIMINATING CONTROL: if IEG predicts and SMAD does not, the
#                association is a protocol artifact, not TGF-beta signalling.
TGFB_SMAD = ["SMAD7", "SKIL", "SKI", "PMEPA1", "KLF10", "BAMBI", "TGIF1"]
TGFB_IEG = ["JUNB", "ID1", "ID2", "ID3", "CCN1", "CCN2", "CDKN1A"]

# Fibroblast / myofibroblast / mural senders (dominant TGF-beta sources)
FIBRO_SENDERS = ["Alveolar fibroblasts", "Adventitial fibroblasts",
                 "Peribronchial fibroblasts", "Subpleural fibroblasts",
                 "Myofibroblasts", "Vascular smooth muscle"]
RECEIVER = "Pericytes"


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="ccc_niche.h5ad")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-receiver", type=int, default=5,
                   help="Min pericytes per donor")
    p.add_argument("--min-sender", type=int, default=10,
                   help="Min fibroblast/mural senders per donor")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Circularity guard, asserted rather than documented: the predictor panels must
    # share no gene with the outcome panel, or the correlation is partly definitional.
    for name, panel in (("TGFB_SMAD", TGFB_SMAD), ("TGFB_IEG", TGFB_IEG)):
        overlap = set(panel) & set(TARGET_PROGRAM)
        if overlap:
            raise SystemExit(
                f"{name} overlaps TARGET_PROGRAM on {sorted(overlap)}: the donor-level "
                "test would be partly circular. Remove the shared genes from one side.")

    backed = ad.read_h5ad(args.adata, backed="r")
    want = sorted(set(LIGANDS + TARGET_PROGRAM + TGFB_SMAD + TGFB_IEG))
    present = [g for g in want if g in backed.var_names]
    logging.info(f"genes present: {len(present)}/{len(want)}")
    col_idx = [backed.var_names.get_loc(g) for g in present]
    obs = backed.obs.copy()
    n = backed.n_obs

    # Chunked row reads (CSR row-slicing is efficient and memory-light); avoids
    # materializing the full 12 GB matrix that backed column-subset would force.
    chunk = 50000
    blocks = []
    for i in range(0, n, chunk):
        Xc = backed[i:min(i + chunk, n)].to_memory().X
        Xc = Xc[:, col_idx]
        blocks.append(Xc.toarray() if sparse.issparse(Xc) else np.asarray(Xc))
        logging.info(f"  read rows {i}-{min(i + chunk, n)}")
    X = np.vstack(blocks)
    del backed, blocks
    expr = pd.DataFrame(X, columns=present, index=obs.index)

    for c in ["dataset", "study"]:
        if c in obs.columns:
            obs[c] = obs[c].astype(str)
    lig = [g for g in LIGANDS if g in present]
    tgt = [g for g in TARGET_PROGRAM if g in present]

    # per-cell summaries
    obs["_target_mean"] = expr[tgt].mean(axis=1).to_numpy()
    obs["_tgfb1"] = expr["TGFB1"].to_numpy() if "TGFB1" in present else np.nan
    # TGFB2 is the top-ranked ligand by corrected AUPR (0.208, just above TGFB1's
    # 0.203), so the donor-level check must be able to test it on its own rather
    # than only inside the composite -- otherwise the leading ligand is the one
    # ligand with no individual donor-replicated estimate.
    obs["_tgfb2"] = expr["TGFB2"].to_numpy() if "TGFB2" in present else np.nan
    obs["_ligand_mean"] = expr[lig].mean(axis=1).to_numpy()

    # Receiver-side pathway activity. Mean over the genes actually present; the
    # count is logged so a shrunken panel is visible rather than silent.
    smad = [g for g in TGFB_SMAD if g in present]
    ieg = [g for g in TGFB_IEG if g in present]
    logging.info(f"TGFB_SMAD present: {len(smad)}/{len(TGFB_SMAD)} {smad}")
    logging.info(f"TGFB_IEG present: {len(ieg)}/{len(TGFB_IEG)} {ieg}")
    obs["_tgfb_smad"] = expr[smad].mean(axis=1).to_numpy() if smad else np.nan
    obs["_tgfb_ieg"] = expr[ieg].mean(axis=1).to_numpy() if ieg else np.nan

    is_recv = obs["ccc_group"].astype(str) == RECEIVER
    is_fibro = obs["ccc_group"].astype(str).isin(FIBRO_SENDERS)

    # receiver target-program expression per donor (all pericytes; RECEIVER above)
    recv = obs[is_recv].groupby("donor_id", observed=True).agg(
        receiver_target_expr=("_target_mean", "mean"),
        # Pathway activity is measured in the RECEIVER, on the same pericytes as
        # the outcome, so it is a within-cell exposure proxy rather than a
        # cross-compartment one.
        receiver_TGFB_SMAD=("_tgfb_smad", "mean"),
        receiver_TGFB_IEG=("_tgfb_ieg", "mean"),
        n_receiver=("_target_mean", "size"),
        disease_group=("disease_group", "first"))
    if "dataset" in obs.columns:
        recv["dataset"] = obs[is_recv].groupby("donor_id", observed=True)["dataset"].first()

    # sender ligand expression per donor (fibroblast/mural senders)
    send = obs[is_fibro].groupby("donor_id", observed=True).agg(
        sender_TGFB1=("_tgfb1", "mean"),
        sender_TGFB2=("_tgfb2", "mean"),
        sender_ligand_mean=("_ligand_mean", "mean"),
        n_sender=("_ligand_mean", "size"))

    donor = recv.join(send, how="inner").reset_index()
    donor = donor[(donor["n_receiver"] >= args.min_receiver) &
                  (donor["n_sender"] >= args.min_sender)].copy()
    donor.to_csv(args.outdir / "donor_validation_table.tsv.gz", sep="\t", index=False)
    logging.info(f"donors with paired sender+receiver: {donor.shape[0]}")


if __name__ == "__main__":
    main()
