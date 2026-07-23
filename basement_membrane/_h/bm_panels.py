"""
Basement-membrane (BM) gene panels -- single source of truth for this module.

The classic lung basement-membrane components, as specified by the collaborator.
Both the Python and R steps read these from `_m/bm_panel_genes.tsv`, which
00.bm_score.py writes, so the panel is defined exactly once (same pattern as
pathway_balance/_h/00.pathway_balance.py -> balance_signature_genes.txt).

Why this module exists: `pericyte_states` folds COL4A1 -- a basement-membrane
collagen -- into the `fibroblast_like` panel alongside fibrillar collagens
(COL1A1/COL1A2/COL3A1), so basement-membrane deposition and fibrillar-ECM
deposition are currently inseparable. The FIBRILLAR_CONTRAST panel below exists
to test that dissociation rather than assume it.

All 13 BM genes were verified present in pericyte_states/_m/pericyte_states.h5ad
(55,329 genes) and disease_association/ipf_analysis/_m/ipf_dataset.h5ad (45,947).
"""

# Full panel. Collaborator's list, mapped to HGNC symbols:
#   collagen IV a1/a2, collagen XVIII, laminin a3/a4/a5, laminin b1/b2,
#   laminin g1, nidogen (entactin), perlecan, agrin. NID2 added for completeness.
BM_PANEL = [
    "COL4A1", "COL4A2", "COL18A1",
    "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMC1",
    "NID1", "NID2", "HSPG2", "AGRN",
]

# Structural sub-panels: network collagens, the laminin heterotrimer arms, and
# the linker/proteoglycan components that bridge the two networks.
BM_SUBPANELS = {
    "bm_collagen_iv": ["COL4A1", "COL4A2"],
    "bm_laminin": ["LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1"],
    "bm_linker": ["NID1", "NID2", "HSPG2", "AGRN", "COL18A1"],
}

# Fibrillar/interstitial ECM. NOT a basement-membrane panel -- this is the
# contrast used to test whether BM deposition is separable from the existing
# fibroblast_like axis. Deliberately excludes every BM_PANEL gene.
FIBRILLAR_CONTRAST = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1",
    "POSTN", "LUM", "DCN", "FBN1", "BGN",
]

# The existing fibroblast_like panel from pericyte_states/_h/00.state_discovery.py,
# minus COL4A1. Quantifies how much of that panel's behaviour was BM contamination.
FIBROBLAST_LIKE_NO_COL4A1 = [
    "COL1A1", "COL1A2", "COL3A1", "FN1", "LUM", "DCN", "PDGFA", "FBLN1",
]

# Pre-specified expectation used as a method sanity check in 03.bm_selectivity_stats.R.
# LAMA3 (laminin-332, with LAMB3/LAMC2) is epithelial and LAMA5 is broadly
# expressed; neither should come out pericyte-selective. The mural set should.
# If the selectivity metric does not recover this layout, the method is wrong.
EXPECTED_NOT_PERICYTE_SELECTIVE = ["LAMA3", "LAMA5"]
EXPECTED_MURAL_VASCULAR = ["LAMA4", "COL4A1", "COL4A2", "HSPG2", "NID1"]

PANELS = {
    "basement_membrane": BM_PANEL,
    **BM_SUBPANELS,
    "fibrillar_ecm": FIBRILLAR_CONTRAST,
    "fibroblast_like_noCOL4A1": FIBROBLAST_LIKE_NO_COL4A1,
}


def panel_table():
    """Long-format panel -> gene table, written to _m/bm_panel_genes.tsv."""
    import pandas as pd
    rows = [{"panel": panel, "gene": gene}
            for panel, genes in PANELS.items() for gene in genes]
    return pd.DataFrame(rows)
