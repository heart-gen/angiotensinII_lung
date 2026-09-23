"""
Receiver programs for the OUTGOING pericyte -> neighbour analysis (Figure 5F),
and the endothelial readout node of the donor-level RAS network (Figure 5C).

Single source of truth: 09.nichenet_outgoing.R and 07.ras_network.R read the
table this module writes (`python receiver_programs.py --out receiver_programs.tsv`)
rather than keeping their own copies.

Design rules, enforced below or at run time:

  1. A receiver program must not contain a gene the SENDER could be scored for.
     NicheNet's lr_network lists collagens, laminins, SPARC, FN1 and others as
     ligands, so a naive "fibroblast matrix" program would contain genes that
     pericytes also express as ligands; a pericyte->fibroblast donor correlation
     on such a gene is partly ambient RNA, not signalling. 09.nichenet_outgoing.R
     drops any program gene that is also a pericyte-expressed potential ligand
     and writes the dropped genes to disk, so the exclusion is visible.
  2. No RAS gene may sit in a receiver program or in EC_READOUT: the network in
     07 regresses these nodes on RAS nodes, and a shared gene would make an edge
     arithmetic. Asserted here.
  3. Programs are RESPONSE genes (junctional, angiogenic, matrix-processing,
     epithelial injury-repair), not identity markers, so a donor-level shift is
     interpretable as a response rather than as a change in cell-type purity.

Provenance of the lists (curated 2026-09-22):
  EC_BARRIER_ANGIOGENIC -- adherens/tight junction and barrier-maintenance genes
    (CDH5, CLDN5, OCLN, TJP1, ESAM, JAM2, S1PR1, TEK) plus the VEGF/Notch
    angiogenic-response set (KDR, FLT1, ANGPT2, ESM1, APLNR, DLL4, HEY1, NRP1).
    Pericyte-endothelial crosstalk regulates exactly these two axes
    (Armulik et al., Dev Cell 2011; Daneman et al., Nature 2010).
  FIBRO_MATRIX -- collagen-processing and fibre-assembly enzymes and matrisome
    organisers rather than the structural collagens themselves (see rule 1):
    LOX/LOXL1/LOXL2 cross-linking, P4HA1/P4HA2/PLOD1/PLOD2 hydroxylation,
    SERPINH1 (HSP47) chaperone, PCOLCE/ADAMTS2/BMP1 procollagen processing,
    ELN/FBLN1/FBLN2/MFAP5 elastic-fibre assembly, CTHRC1.
  EPI_REPAIR -- the alveolar transitional / damage-associated programme
    (KRT8, KRT18, KRT19, KRT7, CLDN4, SFN, SOX4, CDKN1A, TP53I3, LCN2):
    Strunz et al., Nat Commun 2020; Kobayashi et al., Nat Cell Biol 2020;
    Choi et al., Cell Stem Cell 2020.
  EC_READOUT -- EC_BARRIER only (the angiogenic half is a response that could
    also be driven by VEGF from epithelium, so it is kept out of the network node).
"""
import argparse
from pathlib import Path

EC_BARRIER = ["CDH5", "CLDN5", "OCLN", "TJP1", "ESAM", "JAM2", "S1PR1", "TEK"]
EC_ANGIOGENIC = ["KDR", "FLT1", "ANGPT2", "ESM1", "APLNR", "DLL4", "HEY1", "NRP1"]
EC_BARRIER_ANGIOGENIC = EC_BARRIER + EC_ANGIOGENIC

FIBRO_MATRIX = ["LOX", "LOXL1", "LOXL2", "P4HA1", "P4HA2", "PLOD1", "PLOD2",
                "SERPINH1", "PCOLCE", "ADAMTS2", "BMP1", "ELN", "FBLN1", "FBLN2",
                "MFAP5", "CTHRC1"]

EPI_REPAIR = ["KRT8", "KRT18", "KRT19", "KRT7", "CLDN4", "SFN", "SOX4",
              "CDKN1A", "TP53I3", "LCN2"]

EC_READOUT = list(EC_BARRIER)

PROGRAMS = {
    "EC_barrier_angiogenic": EC_BARRIER_ANGIOGENIC,
    "fibroblast_matrix": FIBRO_MATRIX,
    "epithelial_repair": EPI_REPAIR,
    "EC_readout": EC_READOUT,
}

# Which program each outgoing target is scored on. Macrophages and VSMC are
# kept as LIANA targets (edges are reported) but have no pre-specified response
# programme, so NicheNet is not run for them -- inventing one post hoc would
# turn a descriptive edge list into an untested claim.
TARGET_PROGRAM = {
    "EC aerocyte capillary": "EC_barrier_angiogenic",
    "EC general capillary": "EC_barrier_angiogenic",
    "Alveolar fibroblasts": "fibroblast_matrix",
    "Adventitial fibroblasts": "fibroblast_matrix",
    "AT1": "epithelial_repair",
    "AT2": "epithelial_repair",   # AT2_AGTR2det + AT2_AGTR2undet, recombined
}
OUTGOING_TARGETS = [
    "EC aerocyte capillary", "EC general capillary",
    "AT1", "AT2_AGTR2det", "AT2_AGTR2undet",
    "Alveolar fibroblasts", "Adventitial fibroblasts",
    "Alveolar macrophages", "Interstitial macrophages",
    "Vascular smooth muscle",
]

# The RAS genes of agt_axis/_h/ras_panel.tsv (substrate, proteases, receptors).
RAS_GENES = {"AGT", "REN", "ACE", "ACE2", "CMA1", "CTSG", "CTSD", "ENPEP", "MME",
             "AGTR1", "AGTR2", "LRP2", "MAS1"}


def _assert_programs():
    for name, genes in PROGRAMS.items():
        if len(genes) != len(set(genes)):
            raise AssertionError(f"{name} has duplicate genes")
        bad = sorted(set(genes) & RAS_GENES)
        if bad:
            raise AssertionError(
                f"{name} contains RAS genes {bad}; the network in 07 regresses "
                "response nodes on RAS nodes, so a shared gene makes the edge "
                "arithmetic.")
    for t, prog in TARGET_PROGRAM.items():
        if prog not in PROGRAMS:
            raise AssertionError(f"target {t} maps to unknown program {prog}")


_assert_programs()


def table():
    import pandas as pd
    rows = [{"program": p, "gene": g} for p, genes in PROGRAMS.items() for g in genes]
    return pd.DataFrame(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--targets-out", type=Path, default=None)
    a = ap.parse_args()
    table().to_csv(a.out, sep="\t", index=False)
    if a.targets_out is not None:
        import pandas as pd
        pd.DataFrame({"target": OUTGOING_TARGETS,
                      "program": [TARGET_PROGRAM.get(
                          "AT2" if t.startswith("AT2_") else t, "")
                          for t in OUTGOING_TARGETS]}).to_csv(
            a.targets_out, sep="\t", index=False)
    print(f"wrote {a.out}")
