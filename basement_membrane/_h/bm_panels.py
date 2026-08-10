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

2026-08-10 expansion. The fibrillar side was resolved into the two blocks a
collaborator asked for -- fibril-forming types I/III versus the regulatory
types V/XI -- and an ambient-tracer panel was added so the fibrillar signal in
pericytes can be separated from fibroblast-derived soup. The collagen V/XI
chains (COL5A2, COL5A3, COL11A1, COL11A2) and all five tracers were verified
present in both pericyte_states.h5ad and cell_communication/_m/ccc_niche.h5ad.
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
#
# FROZEN 2026-08-10. `fibrillar_ecm_score` feeds the published primary endpoint
# (bm_primary_endpoint_emmeans.tsv), the orthogonalization, and the continuum
# rho values quoted in writings/BM_SUMMARY.md and MECHANISM_ANALYSES.md. Adding
# the collagen V/XI chains here would silently move every one of those numbers.
# New collagen blocks go in FIBRILLAR_CORE / FIBRILLAR_MINOR below; this list
# does not change.
FIBRILLAR_CONTRAST = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1",
    "POSTN", "LUM", "DCN", "FBN1", "BGN",
]

# Fibrillar collagens proper, split into the two blocks the collaborator asked
# for (2026-08-10). The split is structural, not statistical: types I and III are
# the load-bearing fibril-forming collagens that define a fibrotic matrix, while
# types V and XI are low-abundance regulatory chains that nucleate fibrils and
# set their diameter. Types V/XI are therefore corroborating evidence for a
# fibril-forming program, not independent evidence of one -- a cell can carry
# COL5A2 without producing interstitial collagen, but it cannot produce
# interstitial collagen without types I/III.
#
# The matrix-identity claim is made on FIBRILLAR_CORE. FIBRILLAR_MINOR is
# reported alongside it and never carries a claim on its own.
FIBRILLAR_CORE = ["COL1A1", "COL1A2", "COL3A1"]
FIBRILLAR_MINOR = ["COL5A1", "COL5A2", "COL5A3", "COL11A1", "COL11A2"]

# Off-lineage transcripts used to quantify each unit's ambient-RNA burden.
# None of these is transcribed by pericytes, so whatever is detected in a
# pericyte is soup (or a doublet). They are the denominator for the ambient
# controls in 08.fibrillar_ambient.R: if pericyte COL1A2 is ambient, it must
# track these across donors; if it is transcribed, it need not.
#
# Chosen to span the three lineages that dominate lung soup: alveolar epithelium
# (SFTPC/SFTPB), airway epithelium (SCGB1A1/SCGB3A2), and leukocytes (PTPRC).
# Deliberately NOT fibroblast genes -- the fibroblast-derived component is
# measured directly, by regressing on donor fibroblast burden.
AMBIENT_TRACER = ["SFTPC", "SFTPB", "SCGB1A1", "SCGB3A2", "PTPRC"]

# Pre-specified negative control on the fibrillar side, mirroring LAMA3/LAMA5 on
# the BM side. COL11A1 is a fibrotic-fibroblast / CAF gene that is near-absent
# from mural cells; if it comes out pericyte-selective, the metric is wrong.
EXPECTED_NOT_PERICYTE_FIBRILLAR = ["COL11A1"]

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
    "fibrillar_core": FIBRILLAR_CORE,
    "fibrillar_minor": FIBRILLAR_MINOR,
    "ambient_tracer": AMBIENT_TRACER,
    "fibroblast_like_noCOL4A1": FIBROBLAST_LIKE_NO_COL4A1,
}

# Display blocks: a gene belongs to exactly one, unlike PANELS where genes recur
# across overlapping panels. Drives gene ordering and the block separators in
# figures/_h/basement_membrane_figure.R, so the BM -> fibrillar inversion reads
# as a block structure rather than an alphabetical jumble. Order is meaningful.
GENE_BLOCKS = [
    ("basement_membrane", BM_PANEL),
    ("fibrillar_core", FIBRILLAR_CORE),
    ("fibrillar_minor", FIBRILLAR_MINOR),
    ("ambient_tracer", AMBIENT_TRACER),
]


def gene_block(gene):
    """Display block for one gene; `fibrillar_other` for the score-only genes."""
    for block, genes in GENE_BLOCKS:
        if gene in genes:
            return block
    return "fibrillar_other"


def block_order():
    """Genes in display order, block by block."""
    return [g for _, genes in GENE_BLOCKS for g in genes]


def panel_table():
    """Long-format panel -> gene table, written to _m/bm_panel_genes.tsv.

    `block` is constant per gene and carries the display grouping; `block_index`
    fixes the within-block order so downstream R does not re-derive it.
    """
    import pandas as pd
    order = {g: i for i, g in enumerate(block_order())}
    rows = [{"panel": panel, "gene": gene, "block": gene_block(gene),
             "block_index": order.get(gene, 10_000)}
            for panel, genes in PANELS.items() for gene in genes]
    return pd.DataFrame(rows)
