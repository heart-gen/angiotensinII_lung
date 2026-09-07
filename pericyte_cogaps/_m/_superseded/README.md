# Superseded outputs — quarantined 2026-09-07

Files here are **not** produced by any script in `pericyte_cogaps/_h/` and must
not be cited or read by downstream code. Nothing was deleted; move a file back if
you need it.

## `cogaps_seed_stability_summary.tsv` (2026-07-22 09:35)

Orphaned by commit `68846fd` ("Removed orphaned files and updated CoGAPS"). The
current `02.select_rank.R` writes `cogaps_seed_stability.tsv` (151 rows, one per
`ref_pattern` × `seed`, **with an `np` column**) and `cogaps_nP_selection.tsv`.
Nothing writes this summary.

**It is the nP = 5 aggregation, not the selected rank.** Three independent
tells:

- 5 rows, `Pattern_1`…`Pattern_5`, and **no `np` column** — so nothing on its
  face says which factorization it describes;
- its `Pattern_5` mean, `0.700339288223707`, is *exactly*, to 15 digits, the
  nP = 5 `min_r` in `cogaps_nP_selection.tsv`;
- its mtime (2026-07-22 09:35) predates the live outputs (2026-07-23 07:54).

**Why this mattered (defect P1-12).** `tables/_h/04.validation.R` read this file
into supplementary Table **S06B2**, titled *"CoGAPS per-pattern cross-seed
reproducibility at the selected rank"*. The selected rank is **nP = 8**, whose
weakest pattern reproduces at **r = 0.978**; the shipped supplement attached
**0.402** to it — the nP = 5 weakest pattern. The published CoGAPS validation
claim was never in doubt, but the table that documented it described a rank the
manuscript does not use.

`04.validation.R` now reads `cogaps_seed_stability.tsv`, filters to
`np == SELECTED_NP`, aggregates per `ref_pattern`, and **refuses to write the
part** unless that aggregation reproduces the nP = 8 row of S06B1 to 1e-9 and
S06B1 flags nP = 8 as `SELECTED (main)`. Both checks would have failed on this
file, which is the point.
