"""
Basement-membrane (BM) gene panels -- single source of truth for this module.

The classic lung basement-membrane components, as specified by the collaborator.
Both the Python and R steps read these from `_m/bm_panel_genes.tsv`, which
00.bm_score.py writes, so the panel is defined exactly once (same pattern as
pathway_balance/_h/00.pathway_balance.py -> balance_signature_genes.txt).

Why this module exists: `pericyte_states` folded COL4A1 -- a basement-membrane
collagen -- into the `fibroblast_like` panel alongside fibrillar collagens
(COL1A1/COL1A2/COL3A1), so basement-membrane deposition and fibrillar-ECM
deposition were inseparable. The FIBRILLAR_* panels below exist to test that
dissociation rather than assume it.

2026-08-10 expansion. The fibrillar side was resolved into the two blocks a
collaborator asked for -- fibril-forming types I/III versus the regulatory
types V/XI -- and an ambient-tracer panel was added so the fibrillar signal in
pericytes can be separated from fibroblast-derived soup.

2026-09-01 expansion (collaborator request). BM_PANEL grew from 13 to 20 genes:
the collaborator's full list (collagen IV a1/a2, collagen XV, all twelve laminin
chains, nidogen 1/2, perlecan, agrin) PLUS COL18A1, which their list omitted but
which is detected in 52.6% of pericytes and for which pericytes rank #1 of 22
cell types. A combined FIBRILLAR_COLLAGEN score was added so the two matrix
categories are directly comparable as two scores, and a TGF-beta response panel
was added so both categories can be tested against TGF-beta signalling.

BEFORE/AFTER AUDIT. Growing BM_PANEL moves `basement_membrane_score` and
therefore the state gate, the BM x cluster estimates, the BM - fibrillar
endpoint, the continuum rho and the COPD contrasts. That movement is intended,
but it must be auditable: BM_PANEL_V1 preserves the 13-gene panel every
published BM number was computed on, is scored alongside as `bm_v1_score`, and
`bh_family()` keeps the original 13 genes in their own multiple-testing family
so their published per-gene adjusted p-values do not move. `bm_v1_score` exists
ONLY for that audit and must never carry a claim.

Pericyte detection (11,680 cells, logcounts) for the seven added genes, measured
2026-09-01 before the change: LAMC3 57.1%, LAMA2 34.5%, COL15A1 3.0%, LAMB3
0.8%, LAMB4 0.4%, LAMC2 0.3%, LAMA1 0.2%. Only LAMC3 and LAMA2 carry real
signal; the other five are near-absent and dilute the panel score. That is the
known cost of adopting the full structural list, and it is why the gate in
01.state_gate.py is re-run as a hard checkpoint rather than assumed to hold.

All genes below were verified present in pericyte_states/_m/pericyte_states.h5ad
(55,329 genes), cell_communication/_m/ccc_niche.h5ad (55,329) and
disease_association/ipf_analysis/_m/ipf_dataset.h5ad (45,947) on 2026-09-01.
Note CCN2/CCN1 are the current HGNC symbols; the CTGF/CYR61 aliases are absent
from all three objects.
"""

# ---------------------------------------------------------------- BM panel ----
# Collaborator's list mapped to HGNC symbols, in structural order: network
# collagens, the multiplexins, the laminin alpha/beta/gamma arms, then the
# linker/proteoglycan components that bridge the two networks.
BM_PANEL = [
    "COL4A1", "COL4A2", "COL15A1", "COL18A1",
    "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMB3", "LAMB4",
    "LAMC1", "LAMC2", "LAMC3",
    "NID1", "NID2", "HSPG2", "AGRN",
]

# The panel every published BM number was computed on (2026-07-21 .. 2026-08-10).
# Scored as `bm_v1_score` for the before/after concordance table only.
BM_PANEL_V1 = [
    "COL4A1", "COL4A2", "COL18A1",
    "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMC1",
    "NID1", "NID2", "HSPG2", "AGRN",
]

# Structural sub-panels: network collagens, the laminin heterotrimer arms, and
# the linker/proteoglycan components that bridge the two networks.
BM_SUBPANELS = {
    "bm_collagen_iv": ["COL4A1", "COL4A2"],
    "bm_laminin": ["LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
                   "LAMB1", "LAMB2", "LAMB3", "LAMB4",
                   "LAMC1", "LAMC2", "LAMC3"],
    "bm_linker": ["NID1", "NID2", "HSPG2", "AGRN", "COL15A1", "COL18A1"],
}

# --------------------------------------------------------- fibrillar side ----
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

# The collaborator's fibrillar list, which is exactly CORE + MINOR. Scored as one
# panel so "basement membrane" and "fibrillar collagen" can be compared as two
# scores rather than as one score against two sub-blocks. Unlike
# FIBRILLAR_CONTRAST this is collagen only -- no proteoglycans or glycoproteins.
FIBRILLAR_COLLAGEN = FIBRILLAR_CORE + FIBRILLAR_MINOR

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

# ------------------------------------------------------------- TGF-beta ----
# Canonical SMAD2/3 transcriptional response. This is a SIGNALLING readout, and
# it is only interpretable as one if it shares no gene with the things it is
# regressed against -- otherwise "TGF-beta predicts matrix" is arithmetic, not
# biology. Excluded for exactly that reason, and asserted below:
#   SERPINE1, THBS1, POSTN  -- in pericyte_states `activated_migratory`
#   TAGLN                   -- in pericyte_states `synthetic_contractile`
#   every COL*/LAM*/NID*/HSPG2/AGRN and FN1/LUM/DCN/FBN1/BGN -- matrix panels
# TGFBI is retained: it is a canonical TGF-beta-induced gene and is in no panel.
TGFB_RESPONSE = [
    "SMAD7", "SKIL", "SKI", "JUNB", "PMEPA1", "TGFBI",
    "ID1", "ID2", "ID3", "CDKN1A", "SERPINE2",
    "CCN2", "CCN1", "KLF10", "BAMBI", "TGIF1", "SNAI1",
]

# Leave-one-out sensitivity. TGFBI (TGF-beta-induced protein ig-h3) is a canonical
# SMAD target AND a secreted extracellular-matrix protein, so it is the one member
# of the panel whose correlation with a matrix score could be partly definitional.
# The claim is made on TGFB_RESPONSE and must hold on this panel too; if the two
# disagree, the association is carried by an ECM gene and cannot be reported as a
# signalling result.
TGFB_RESPONSE_NO_ECM = [g for g in TGFB_RESPONSE if g != "TGFBI"]

# Receptor/transducer availability. Descriptive only -- reported so a reader can
# see that pericytes can receive TGF-beta at all; never used as a response score.
TGFB_RECEPTOR = ["TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD3", "SMAD4"]

# ---------------------------------------------------------------- controls ----
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
    "bm_v1": BM_PANEL_V1,
    **BM_SUBPANELS,
    "fibrillar_ecm": FIBRILLAR_CONTRAST,
    "fibrillar_collagen": FIBRILLAR_COLLAGEN,
    "fibrillar_core": FIBRILLAR_CORE,
    "fibrillar_minor": FIBRILLAR_MINOR,
    "ambient_tracer": AMBIENT_TRACER,
    "tgfb_response": TGFB_RESPONSE,
    "tgfb_response_noECM": TGFB_RESPONSE_NO_ECM,
    "tgfb_receptor": TGFB_RECEPTOR,
    "fibroblast_like_noCOL4A1": FIBROBLAST_LIKE_NO_COL4A1,
}

# Display blocks: a gene belongs to exactly one, unlike PANELS where genes recur
# across overlapping panels. Drives gene ordering and the block separators in
# figures/_h/basement_membrane_figure.R, so the BM -> fibrillar inversion reads
# as a block structure rather than an alphabetical jumble. Order is meaningful.
#
# `tgfb_response`/`tgfb_receptor` are here so 00.bm_score.py emits their per-gene
# values; the matrix figure filters to the four matrix blocks explicitly and does
# not render them.
GENE_BLOCKS = [
    ("basement_membrane", BM_PANEL),
    ("fibrillar_core", FIBRILLAR_CORE),
    ("fibrillar_minor", FIBRILLAR_MINOR),
    ("ambient_tracer", AMBIENT_TRACER),
    ("tgfb_response", TGFB_RESPONSE),
    ("tgfb_receptor", TGFB_RECEPTOR),
]

# Multiple-testing families, PRE-SPECIFIED and frozen at their original
# membership. If the 2026-09-01 BM additions joined the original family, every
# published per-gene p_BH in bm_selectivity_emmeans.tsv would shift purely
# because the panel grew -- a silent revision of results the manuscript already
# quotes. The R steps read this from the `bh_family` column of
# bm_panel_genes.tsv rather than re-deriving it, so the partition is declared
# once here. Strings are stable; downstream code matches on them.
BH_FAMILIES = [
    ("original_panel", list(BM_PANEL_V1) + list(FIBRILLAR_CONTRAST)),
    ("ambient_tracer", AMBIENT_TRACER),
    ("tgfb_response", list(TGFB_RESPONSE) + list(TGFB_RECEPTOR)),
    ("bm_expansion", [g for g in BM_PANEL if g not in BM_PANEL_V1]),
]


def _assert_tgfb_disjoint():
    """The TGF-beta response score must share no gene with anything it is
    regressed against, or the association is circular. Enforced, not just
    documented -- the whole module exists because a matrix gene once sat in the
    wrong panel unnoticed."""
    matrix = set(BM_PANEL) | set(FIBRILLAR_CONTRAST) | set(FIBRILLAR_COLLAGEN)
    # pericyte_states STATE_PANELS, mirrored here so this file has no import
    # cycle with 00.state_discovery.py (which imports THIS module).
    state_panel_genes = {
        "RGS5", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "HIGD1B", "COX4I2",
        "NDUFA4L2", "GJA4", "CSPG4", "PLXDC1",
        "IL6", "CCL2", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CXCL6", "CXCL8",
        "CXCL10", "IL1A", "IL1B", "MIF", "ICAM1", "VCAM1", "SELE", "NFKBIA",
        "ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "DES", "VIM",
        "ADAMTS1", "THBS1", "TIMP1", "MMP2", "MMP3", "MMP9", "SERPINE1", "POSTN",
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "LUM", "DCN", "PDGFA",
        "FBLN1",
    }
    bad_matrix = sorted(set(TGFB_RESPONSE) & matrix)
    bad_state = sorted(set(TGFB_RESPONSE) & state_panel_genes)
    if bad_matrix or bad_state:
        raise AssertionError(
            "TGFB_RESPONSE overlaps the panels it is tested against: "
            f"matrix={bad_matrix}, state={bad_state}. Remove them -- an "
            "overlapping gene makes the TGF-beta association arithmetic.")


def _assert_panels_consistent():
    if len(BM_PANEL) != len(set(BM_PANEL)):
        raise AssertionError("BM_PANEL has duplicate genes")
    missing_v1 = sorted(set(BM_PANEL_V1) - set(BM_PANEL))
    if missing_v1:
        raise AssertionError(
            f"BM_PANEL_V1 genes absent from BM_PANEL: {missing_v1}. The v1 panel "
            "must be a subset, or the concordance table compares two panels that "
            "differ in both directions and the audit is uninterpretable.")
    sub = {g for genes in BM_SUBPANELS.values() for g in genes}
    if sub != set(BM_PANEL):
        raise AssertionError(
            "BM_SUBPANELS do not partition BM_PANEL; symmetric difference: "
            f"{sorted(sub ^ set(BM_PANEL))}")
    if set(FIBRILLAR_COLLAGEN) & set(BM_PANEL):
        raise AssertionError("fibrillar and BM panels overlap")
    fam = [g for _, genes in BH_FAMILIES for g in genes]
    if len(fam) != len(set(fam)):
        raise AssertionError("a gene appears in >1 BH family")


_assert_panels_consistent()
_assert_tgfb_disjoint()


def gene_block(gene):
    """Display block for one gene; `fibrillar_other` for the score-only genes."""
    for block, genes in GENE_BLOCKS:
        if gene in genes:
            return block
    return "fibrillar_other"


def bh_family(gene):
    """Multiple-testing family for one gene (see BH_FAMILIES)."""
    for family, genes in BH_FAMILIES:
        if gene in genes:
            return family
    return "fibrillar_expansion"


def block_order():
    """Genes in display order, block by block."""
    return [g for _, genes in GENE_BLOCKS for g in genes]


def panel_table():
    """Long-format panel -> gene table, written to _m/bm_panel_genes.tsv.

    `block` and `bh_family` are constant per gene and carry, respectively, the
    display grouping and the multiple-testing partition; `block_index` fixes the
    within-block order so downstream R does not re-derive it.
    """
    import pandas as pd
    order = {g: i for i, g in enumerate(block_order())}
    rows = [{"panel": panel, "gene": gene, "block": gene_block(gene),
             "block_index": order.get(gene, 10_000),
             "bh_family": bh_family(gene)}
            for panel, genes in PANELS.items() for gene in genes]
    return pd.DataFrame(rows)
