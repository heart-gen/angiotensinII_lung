# Superseded outputs — quarantined 2026-07-22

These files are **not** produced by the current `_h/01.state_stats.R` and must not
be cited. They were left behind by earlier versions of that script and by the
pre-relabel run. Nothing was deleted; move a file back if you need it.

Two distinct reasons:

## 1. Superseded by the basement-membrane relabel

- `composition_program_fibroblast_like_{emmeans,posthoc}.tsv` (2026-06-16)

`00.state_discovery.py` was re-run on 2026-07-21 with the sixth (basement-membrane)
panel, which renamed the program covering clusters 1/3/5. `01.state_stats.R` was
re-run on 2026-07-22 and now writes
`composition_program_basement_membrane_{emmeans,posthoc}.tsv` in their place. The
underlying cells and numbers are unchanged — cluster assignment agrees for
11,680/11,680 cells and all five original program scores are bit-identical
(max |difference| = 0) against `_m_backup_5panel/` — so these two files are
correct arithmetic under a label that no longer exists.

## 2. Orphaned by an earlier refactor of `01.state_stats.R`

- `composition_<program>_{emmeans,posthoc}.tsv` (no `state`/`program` tag)
- `composition_disease_anova_all.tsv`, `composition_by_disease.{pdf,png}`
- `state_score_*.tsv`
- `injury_score_*.{tsv,pdf,png}`

The current script writes only `composition_state_*`, `composition_program_*` and
`injury_fraction_*`. It contains no per-program *score*-by-disease analysis at all,
so `state_score_*` has no generating code in the repository. These files date from
2026-06-15/16 and cannot be regenerated as-is.

**Open question for the manuscript:** `MECHANISM_ANALYSES.md` states that "the
per-program state scores are likewise flat (all BH-adjusted p ≥ 0.79)". The
closest surviving file, `state_score_disease_anova_all.tsv`, reports BH ≥ 0.98 —
so the quoted figure comes from a still-earlier run than any output on disk.
Either restore the score-by-disease analysis to `01.state_stats.R` and re-run, or
drop the sentence. Do not update the number from these files.
