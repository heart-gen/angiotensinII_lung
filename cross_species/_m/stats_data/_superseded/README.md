# Superseded `cross_species` outputs

**Do not cite anything in this directory.** These nine files are the complete
output of `cross_species/_h/03.conservation_stats.R` (run 2026-07-21 22:11–22:12),
which was retired on 2026-07-22 and is no longer invoked by `step_3.sh`. They are
kept for provenance only.

## Why they were retired

`04.species_comparability.py` showed that neither test in `03` measures what its
name says:

- The n = 1,144 "mural" set it operates on is **96.4 % smooth muscle**
  (1,016 vSMC + 87 PA-SMC + **41 pericytes**), and the mouse `basement_membrane`
  state contains **zero** pericytes. So the injury-vs-stabilizing contrast is a
  smooth-muscle contrast wearing pericyte-state labels.
- Cell type is perfectly aliased with dataset, so `dataset` as a fixed effect
  absorbs the pericyte-vs-vSMC comparison entirely.
- The continuous `Agtr1a ~ vascular_stabilizing_score` test runs on the dense
  `scvi_corrected` layer, where every cell is non-zero, so it largely restates
  cell-type identity rather than measuring conservation.

**No mouse state-level comparison is possible from these data.** The defensible
cross-species claim is compartment-level, lives in `04.species_comparability.py`,
and rests on raw counts: `Agtr1a` detected in 31 of 41 mouse pericytes versus
0 of 87 PA-SMC.

## What is live instead

The ten `species_comparability_*.tsv` files in the parent directory. The
manuscript supplement figure reads exactly two of them —
`species_comparability_mural_cells.tsv` and `species_comparability_agtr1a_tests.tsv`
(`figures/_h/manuscript_mechanism_figure.R:512–513`).

`mouse_Agtr1a_by_state.pdf` is **not** in any figure; it was retired at the same
time (`figures/_h/assemble_mechanism_figures.R:119`).

## Contents

| File | Was |
| ---- | --- |
| `mouse_Agtr1a_by_state.tsv` / `.pdf` / `.png` | Per-state `Agtr1a` summary. The script annotated this "consumed by the manuscript supplement figure" — it never was; that comment was removed 2026-09-07 (P2-15). |
| `mouse_Agtr1a_injury_vs_stable.tsv` | Categorical injury-vs-stabilizing contrast (the smooth-muscle contrast above). |
| `mouse_Agtr1a_vs_stabilizing_continuous.tsv` | Continuous score test on the dense scVI layer. |
| `mouse_Agtr1a_per_dataset_consistency.tsv` | Per-dataset breakdown of the above. |
| `mouse_injury_grouping_used.tsv` | Which states were counted as injury. |
| `mouse_mural_state_composition.tsv` | Mural state fractions — the table that shows the 96.4 % problem. |
| `mouse_state_composition_by_dataset.tsv` | The same, split by dataset. |

Moved here 2026-09-07 (`writings/TODO.md` P2-15). The convention follows
`pericyte_states/_m/stats_data/_superseded/` and `pericyte_cogaps/_m/_superseded/`.
