# `copd/` — RETIRED

**Status:** fully retired. **No exemptions** — nothing in this directory is read
by any current module. (`grep` for `copd/_m` and `COPD_dataset` across every
`.R`/`.py`/`.sh` returns only this module's own files.)

Written 2026-09-07 as part of defect **P1-18** (`writings/TODO.md`). Contrast
with `../ipf_analysis/RETIRED.md`, which *does* carry live exemptions.

---

## What this was

`_h/01.angiotensinII_analysis.R`, launched by `step_1.sh`. Ran on **JHPCE**
(Rocky Linux 9.2, R 4.3.1, user `jbenjami`) on **2024-04-25**, for the preprint.
It cannot run here: `--partition=bluejay,shared` and `module load R` are
JHPCE-only.

**This is a different cohort from everything else in the project** — and that is
the most important thing to know about it. SmartSeq2 `lungsoup` from
`inputs/copd/_m/COPD_dataset.RData`, with **39 cell-type labels** including
`Fibroblast CTHRC1+` and `Aberrant Basaloid`, and **no Myofibroblast class** at
all. It is *not* GSE136831, which is what `disease_association/agtr1_copd_ipf/`
uses. Donor and cell counts are recorded nowhere.

## Results — all null

Cell-type panel *P*: **0.65, 0.33, 0.78, 0.092**. Combined with `ipf_analysis/`,
the minimum over ~30 uncorrected 2024 tests is **0.092**. There was no
significant 2024 disease effect, so the current analyses never contradicted a
2024 *statistical* claim.

Venn diagram: **317 *AGTR1*⁺, 67 *AGTR2*⁺, 0 co-expressing.**

## Known defects

- **No tabular output survives.** `_h/01:179` writes `normalized_expression.tsv`;
  that file exists nowhere in the repository. Every number here is locked in PDF
  vector text and **cannot be regenerated**. (P2-30.)
- **Two hand-picked four-type panels ship, out of 39 cell types**, with the
  all-cell-type call **commented out at line 169**. It changed no conclusion —
  nothing was significant — but the shipped figure is a subset chosen after the
  fact. (P3-18.)

## Credit

Donor-level aggregation before every test — no pseudoreplication. Both
`all_cells/` and `filter_cells/` variants emitted and kept. Complete
`session_info()` in `_m/summary.log`.

## If you are here to clean up

Retire the whole directory, but do not **delete** it: `_m/all_cells/`,
`_m/filter_cells/`, `_m/angiotensinii_venn_diagram.pdf` and `_m/summary.log` are
the only surviving record of the preprint's COPD numbers.
