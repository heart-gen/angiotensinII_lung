# `pericyte_analysis/` — RETIRED. Do not cite any output in this directory.

**Status:** retired 2026-09-07 (defects **P2-22** / **P2-23**, `writings/TODO.md`;
commit `827769d`). This file written 2026-09-08 to make that decision
discoverable from the directory itself rather than only from a commit message,
and to close **P2-21** in the same place.

Nothing here is a live dependency of any other module. Unlike
`../ipf_analysis/RETIRED.md`, there are **no exemptions** — this directory can be
ignored in full.

---

## Why it is retired

Three independent grounds, each measured rather than asserted. Full evidence in
`_m/transfer_diagnostics/` (written by `_h/05.transfer_diagnostics.py`).

### 1. The matrices are inverted from what the pipeline intended (P2-23)

`01.subset_data.py` called `normalize_total(adata, layer="counts")`, which
normalises the counts **layer in place** and leaves `X` alone — the opposite of
the intent on both matrices:

- `layers["counts"]` became CP10K-normalised and non-integer, **and that is the
  layer scANVI was set up on** — hence "does not contain unnormalized count
  data" and a negative-binomial likelihood evaluated off its support.
- `X` was `log1p`'d while still **raw**, so every HVG call, PCA and neighbour
  graph in this module ran on `log1p(raw counts)` with no depth normalisation.
  Shipped `X` reaches **1290.76**.

The second consequence appears in no defect entry and is the larger one.

### 2. The transfer does not reproduce the query's own structure (P2-22)

`_m/transfer_diagnostics/concordance.json`, on 670 cells:

| Metric | Value |
| --- | ---: |
| Adjusted Rand index (Leiden × transferred) | **0.369** |
| Normalised mutual information | **0.466** |
| Median per-cell confidence | **0.628** |
| Fraction below 0.5 confidence | **31.3%** |

Two of six transferred labels are 8 and 3 cells, at confidence 0.385 and 0.222.

### 3. It cannot be repaired from what is on disk

The query side is fixable and **was** fixed in `827769d` — upstream `X` holds raw
integers, so the counts were destroyed, not missing. The reference side is not:
`ref_hvg` `layers["counts"]` is bit-identical in range to `layers["soupX"]`
(0.0021–908.71, non-integer). **The reference's "counts" is the ambient-corrected
matrix, and no raw count matrix for it exists in this repository.** A rebuild is
blocked on regenerating the reference from original counts.

The code was corrected anyway so that a future rebuild starts from a correct
state, and a seed and training-history capture were added — no fit here was ever
reproducible, and the shipped run stopped at **epoch 11 of 500** with no learning
curve on file.

---

## P2-21 — the confidence threshold that was computed and never applied

Recorded here rather than fixed, because the module it applies to is retired.

`03.transfer_labels.py` calls `generate_filtering_summary(query_full, 0.75, ...)`,
writes per-label retention to `results/annotation_filtering_summary.tsv`, and then
writes and uses `query_full` **unfiltered**:

| Predicted label | Cells | ≥ 0.75 | % |
| --- | ---: | ---: | ---: |
| 0 | 237 | 100 | 42.2 |
| 1 | 141 | 65 | 46.1 |
| 2 | 211 | 45 | 21.3 |
| 3 | 8 | 0 | **0.0** |
| 4 | 70 | 6 | **8.6** |
| 6 | 3 | 0 | **0.0** |
| **Total** | **670** | **216** | **32.2** |

Retention collapses for the small reference clusters, and clusters 5 and 7
receive **zero** query cells — the signature of a classifier defaulting to
majority classes where the query holds no clean analogue.

**Also missing, and it is the more important gap:** no held-out reference cells
were ever scored, so **no accuracy number exists anywhere for this transfer**.
Confidence is the model's self-assessment, not its correctness. That is why
applying the 0.75 filter would not have rescued the module — it would have
produced a smaller set of labels with the same unknown error rate.

**Resolution: stated, not applied.** The labels are unfiltered, the module is
retired, and the reason not to filter is that filtering would imply the remaining
labels were validated. They were not.

---

## If this is ever revived

1. Regenerate the reference from raw counts (the current blocker).
2. Keep the `01.subset_data.py` fix: `normalize_total(adata)` then `log1p(adata)`
   acting on `X`, leaving `layers["counts"]` integral.
3. Score held-out reference cells and report an accuracy, before any label is
   used for anything.
4. Re-read `_m/transfer_diagnostics/` and require ARI well above 0.369 before
   treating the transferred labels as informative.
