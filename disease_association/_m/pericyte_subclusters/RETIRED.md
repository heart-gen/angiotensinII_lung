# RETIRED — do not cite

Everything in this directory was produced by
`disease_association/_h/02.pericytes_disease_analysis.R` on **2026-01-05** and is
retired as of **2026-09-07** (`writings/TODO.md` P2-24).

## Why

The analysis is *AGTR1* by `leiden_pericytes` subcluster. It has been superseded
twice:

- **`pericyte_states/`** replaced the leiden subclusters with the functional state
  model (vascular-stabilizing, inflammatory, synthetic-contractile,
  activated-migratory, fibroblast-like, basement-membrane). Subcluster identity
  is not the unit any current claim is made on.
- **`disease_association/_h/05.agtr1_celltype_disease.R`** replaced the disease
  contrast with a per-cell-type analysis on a donor-aware model.

## The defect it carries

`02.pericytes_disease_analysis.R:57` still has `filter(age > 20)`. `dplyr::filter`
drops `NA`, and HLCA reports no age for **89 %** of Fibrotic/ILD stroma cells, so
this is a cohort filter rather than an age restriction — the same defect family as
P1-2, fixed in `01.disease_association.R` on 2026-07-30 and never backported here.
It is left in place on purpose: repairing a retired script would suggest its
outputs are usable.

## Its input is also retired

The script reads
`localization/pericyte_analysis/_m/pericyte_with_embeddings.h5ad`, from the
2024-preprint-era module that P1-17 marked **do not cite** (superseded labels,
un-logged `X`, 670 pericytes with 24 Control cells).

## What to use instead

| Question | Live source |
| --- | --- |
| *AGTR1* by pericyte state | `pericyte_states/_m/stats_data/` |
| *AGTR1* by cell type and disease | `disease_association/_m/mean_expr/agtr1_celltype_disease_effects.tsv` |
| Injury-program disease effects | `disease_association/_m/mixed_model_forest/` |

`step_2.sh` is deliberately absent from `submit_pipeline.sh` (P2-1).
`DISEASE_SUMMARY.md` is the only document that still cites this directory and
should stop.
