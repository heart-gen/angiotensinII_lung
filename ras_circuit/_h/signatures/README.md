# AngII response signature — provenance and pre-specified rule

Written **2026-09-22, before any differential expression was run**. The rule below
fixes what the signature is, so it cannot be tuned to the lung result it is later
scored against.

## Source

McLellan MA, et al. *High-resolution transcriptomic profiling of the heart during
chronic stress reveals cellular drivers of cardiac fibrosis and hypertrophy.*
Circulation 2020;142:1448–1463. doi:10.1161/CIRCULATIONAHA.119.045115.
Data: ArrayExpress **E-MTAB-8810**, processed file `full_count_matrix.tsv`
(genes × 29,615 cells; barcode suffix `_N` = sample `AP1800N`).

| Sample | Treatment (sdrf `Factor Value[compound]`) |
| --- | --- |
| AP18003, AP18005, AP18007, AP18009 | angiotensin II, 1.5 mg/kg/day, 2 weeks, osmotic pump |
| AP18004, AP18008 | saline (vehicle pump) |
| AP18002, AP18006 | none (untreated) |

Why this source: it is the only public dataset found that gives an **in vivo**,
AngII-versus-vehicle, **single-cell** response in which **pericytes** are resolved
as their own population (the paper lists pericytes, smooth muscle and fibroblasts
among the recovered cell types). Bulk VSMC AngII RNA-seq (e.g. Das et al., Nat
Commun 2017) is a different lineage in culture; the human podocyte AngII ±
losartan study (IJMS 2023) is the wrong lineage.

Known limits, stated up front: mouse, heart, 2-week infusion; AngII infusion
also raises blood pressure, so the response includes haemodynamic load; controls
are pooled vehicle + untreated (below). The signature is a **reference response
programme**, not a lung measurement.

## Pre-specified rule (`00b.derive_angii_signature.py`, `00c.map_orthologs.R`)

1. **QC**: cells with ≥ 200 genes and < 30% mitochondrial reads.
2. **Annotation**: Leiden clusters (resolution 1.0 on 30 PCs of 2,000 HVGs) are
   labelled by the highest mean marker score among pericyte (*Kcnj8, Abcc9, Rgs5,
   Vtn, Higd1b, Pdgfrb, Cspg4, Notch3*), smooth muscle (*Myh11, Acta2, Cnn1,
   Tagln, Des*), fibroblast (*Col1a1, Pdgfra, Dcn, Gsn, Tcf21*), endothelial,
   macrophage and cardiomyocyte panels. Labels are logged per sample.
3. **Unit**: raw-count pseudobulk per sample × cell type, ≥ 20 cells per unit.
4. **Contrast**: DESeq2 (pydeseq2), `~ condition`, **AngII (n = 4) vs control
   (saline + untreated, n = 4)**. Saline alone is n = 2, below the ≥ 3-per-arm
   floor, so it is the sensitivity contrast (AngII vs saline only), reported for
   concordance of sign, never used to build the signature.
5. **Selection**: padj < 0.05 and |log2FC| ≥ 0.5; direction kept; at most 200
   genes per direction, by padj.
6. **Cell type, in order**: pericytes, then mural (pericyte + smooth muscle
   pooled), then fibroblasts. A cell type needs ≥ 3 pseudobulk units per arm, and
   is used if ≥ 30 genes pass step 5 (twice the step-8 floor, to leave room for
   ortholog loss and pruning); otherwise the next is tried. If none reaches 30, the
   first cell type with ≥ 3 units per arm is used. If step 8 then leaves < 15 genes,
   `03` stops and the frozen file must be deleted to re-derive. The cell type used
   is written into the signature file.
   *(Clarified 2026-09-22 before any DE ran, to state the implemented rule exactly.)*
7. **Orthologs**: mouse → human with `babelgene`, strict one-to-one only.
8. **Disjointness** (applied in `03.at1r_response_score.py`): genes in any
   pericyte state panel, any basement-membrane / fibrillar / tracer / TGF-β panel,
   or the RAS panel (incl. *AGTR1*) are removed, and each removal is logged in
   `stats_data/at1r_signature_audit.tsv`. ≥ 15 genes must remain.

The first derivation is **frozen** into `angii_response_signature.tsv` in this
directory; later runs read the frozen copy and never overwrite it. Delete it
deliberately to re-derive.
