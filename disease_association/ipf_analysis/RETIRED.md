# `ipf_analysis/` — RETIRED, with two explicit exemptions

**Status:** the *analysis* in this directory is retired. **The data build is not.**
Written 2026-09-07 as part of defect **P1-18** (`writings/TODO.md`).

---

## ⚠️ DO NOT DELETE — two files here are live dependencies

| Exempt file | Why |
| ----------- | --- |
| **`_h/02.generate_ipf_data.R`** | The only builder of `_m/ipf_dataset.h5ad` (8.4 GB, 312,928 × 45,947, **raw counts**) |
| **`_h/step_2.sh`** | Its only launcher — and the only one in this directory ever ported to Bridges-2 |
| **`_h/sample_demo.csv`** | The GEO cohort demographic table |

**Four current modules read them:**

| Consumer | Reads |
| -------- | ----- |
| `basement_membrane/_h/step_3.sh` (→ `05.bm_copd.py`) | `ipf_dataset.h5ad`, `sample_demo.csv` |
| `disease_association/agtr1_copd_ipf/_h/{step_1.sh, step_1b.sh, 02.unit_decomposition.py}` | `ipf_dataset.h5ad`, `sample_demo.csv` |
| `disease_association/pericyte_analysis/_h/{01.subset_data.py, 03.transfer_labels.py}` | `ipf_dataset.h5ad` |
| `tables/_h/01.cohort_mouse.R` (Table S1) | `sample_demo.csv` |

So the natural housekeeping action — retiring the directory as a unit — would
remove a build step four live modules depend on. That is the whole reason this
file exists rather than a plain `RETIRED.md`.

**`_h/02` applying no QC is the correct call**, not an omission: it is a format
conversion, and the QC belongs to each consumer. This also means the live h5ad is
clean of both bugs recorded below.

### The rebuild path now works

As of 2026-09-07, `step_2.sh` can actually run under `sbatch`. Two things were
wrong with it:

- It pointed at `bio260021p/shared/opt/env/R_env`, **which does not exist**
  (fixed 2026-09-06 in the 22-script sweep, commit `08071bf`).
- It had **no `source /etc/profile.d/modules.sh` / `conda shell.bash hook`
  bootstrap**, so `conda activate` could not resolve under a batch shell
  (fixed 2026-09-07).

It also reported `"Error: Python execution failed"` when an `Rscript` call failed
— copied boilerplate, now corrected. Same class as the `pericyte_analysis`
launcher bug in P1-17.

---

## What *is* retired

`_h/01.angiotensinII_analysis.R` and `_h/03.AGTR1_analysis_IPF.R`, and their
launchers `step_1.sh` / `step_3.sh`. Both ran on **JHPCE** (Rocky Linux 9.2,
R 4.3.1, user `jbenjami`) in **April–May 2024**, for the preprint. Neither can
run here: `--partition=bluejay,shared` and `module load R` are JHPCE-only.

**`_h/03` cannot be re-run anywhere.** It reads `normalized_expression.tsv.gz`,
written by `_h/01`. That file does not exist anywhere in the repository, so every
2024 number is locked in the vector text of the shipped PDFs. (Tracked as P2-30.)

### The 2024 results were null, and that matters

Extracted from the shipped PDFs — there are **no result tables**:

- cell-type panel *P*: 0.82, 0.43, 0.20, 0.32, 0.12, 0.92, 0.69, 0.45, 0.32,
  0.78, 0.33, 0.45, 0.20, 0.95, 0.19
- `_h/03` whole-lung ANOVA: *F*₂,₇₅ = 2.498, ***P* = 0.089**

**Minimum over ~30 uncorrected tests, across this module and `copd/`: 0.092.**

So **2024 and the 2026 analyses were never in *statistical* conflict** — there
was no significant 2024 disease effect for the current work to contradict. Any
disagreement was always about *descriptive direction* read off dot plots and the
detection-conditioned `filter_cells/` boxplots.

`_h/03`'s "Positive" proportions (Control 0.750 / COPD 0.778 / IPF 0.938) are
`mean(AGTR1 across all cells of a donor) > 0.001` — a threshold on a **whole-lung
average**, so donors with more sequenced cells clear it more easily and the IPF
arm contributes the most cells. **This is not a prevalence** and must never be
quoted as one.

### Two real bugs, confined to `_h/01`

1. **The QC filter is an `OR` where it must be an `AND`, and `discard` holds the
   KEEP mask.** `qc.lib <- nUMI > 1000` and `qc.mito <- mito% < 25` are both
   *pass* masks; the code does `discard <- qc.mito | qc.lib; sce <- sce[, discard]`.
   The log announces 56 + 1,082 removals; the code removes only cells failing
   **both** — at most 56 of 312,928. (P2-29.)
2. **`grep("^MT", rownames)` is unanchored.** 430 genes match, only 34 are `MT-`:
   **396 false positives**, including `MTOR`, `MTHFR`, the `MTMR` family and the
   highly expressed metallothioneins `MT1A`/`MT1X`/`MT2A`. `subsets_Mito_percent`
   is therefore not a mitochondrial fraction.

Also unaddressed: `computeSumFactors` warned `non-positive size factor estimates`
and continued; those cells yield NaN log-counts that propagate through
`na.rm = TRUE`.

**Blast radius is zero for current work** — `_h/02` applies no QC, so the live
h5ad and all four consumers are unaffected.

### Credit where it is due

Donor-level aggregation preceded every human test — **no pseudoreplication
anywhere**, in 2024. Both `all_cells/` and `filter_cells/` variants were emitted
and kept, which is the only reason the later sign flip could be decomposed at all.
All three logs carry a complete `session_info()`.

---

## If you are here to clean up

- **Keep** `_h/02.generate_ipf_data.R`, `_h/step_2.sh`, `_h/sample_demo.csv`,
  `_h/sample_demo.md`, `_m/ipf_dataset.h5ad`, `_m/conversion.log`.
- **Retire** `_h/01`, `_h/03`, `step_1.sh`, `step_3.sh`, `_m/all_cells/`,
  `_m/filter_cells/`, `_m/angiotensinii_venn_diagram.pdf`, `_m/analysis.log`,
  `_m/summary.log` — but do not *delete* them: they are the only surviving record
  of the preprint numbers.
- See also `../copd/RETIRED.md`, which has no exemptions.
