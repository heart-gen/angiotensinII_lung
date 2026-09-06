# Figure regeneration after the P1 fixes (opened 2026-09-02)

Seven P1 defects were fixed on 2026-09-02 and each changed tables that figures
read. **No figure in `figures/mechanism/` has been regenerated since.** Until
that happens, several figures plot superseded numbers, and two plot a
*direction* that has since inverted.

This file is the task list. It is scoped to the figure layer only; the analysis
side is closed (see `writings/TODO.md` P1-1, P1-2, P1-3, P1-5, P1-7, P1-13,
P1-14).

---

## What changed upstream, and why it reaches figures

| Fix | Tables that moved | Nature of the change |
| --- | --- | --- |
| **P1-3** DPT root | `pseudotime_trend_correlations.tsv`, `continuum_metadata.tsv.gz`, `bm_continuum_summary.tsv` | **Direction inverted.** All 16 trend correlations flip sign. Two endpoints lose significance entirely. |
| **P1-1 / P1-2** age filter + study RE | `composition_*`, `injury_fraction_*`, `AGTR1_by_*` | Fitted N 47 → **89**; model gains `(1 \| study)`; **COPD (*n* = 1) reappears as a fourth group**. |
| **P1-7** + airspace | `niche_index_per_donor.tsv.gz`, `airspace_donor_summary.csv` | Stability arm goes from 1 component (46 donors) to **2 components (89 donors)**. |
| **P1-5 / P1-14** contrast signs | `bm_*_posthoc`, `bm_selectivity_*`, `agt_source_posthoc`, `agtr1_copd_*` | `t.ratio` and labels corrected in 9 tables / 2,655 rows. |
| **P1-13** reference group | `injury_pattern_disease_np{8,9}.tsv` | Baseline COPD → **Healthy**; BH columns added. |
| TGF-β / AGTR1 nulls | `tgfb_specificity_*`, `agtr1_null_*`, `agtr1_count_null_*` | New panel already added (S17E); legends rewritten. |

---

## Figures to regenerate

Ordered by how badly the current version misleads.

### Tier 1 — currently plots an inverted direction

- [ ] **`figure_pericyte_layer`** (Figure 2) — panel F is the continuum lollipop.
      Every feature's ρ has flipped sign. **The panel currently shows the
      continuum running the wrong way.**
      *Also check:* the panel-F guard added for `AGTR1_scvi` still fires.
- [ ] **`figureS_pericyte_layer`** — same script, same continuum inputs.
- [ ] **`figure_basement_membrane`** — panel G is per-donor continuum ρ by
      metric. `rho_bm` is now a **null** (−0.014, BH = 0.715, was −0.177 at
      *P* = 0.013) and `rho_switch` is **n.s.** (BH = 0.061, was 0.013). Panels
      that annotate significance must be re-checked, not just redrawn.

### Tier 2 — data moved substantially

- [ ] **`figureS_state_composition`** (S11) — was fitted on 47 donors and
      reported as 93; now 89 donors with `(1 | study)`. **A fourth disease group
      (COPD, *n* = 1) now exists** and drives cluster 5 to BH = 0.0005 on its
      own. The figure must show per-group N and must not present the COPD
      contrast as a finding — use `p_excl_small_groups`.
      *Its legend in `mechanism/README.md` is already flagged as superseded.*
- [ ] **`figure_mechanism_main`** — reads `pseudotime_trend_*` and
      `niche_index_per_donor`. Both changed. (Panel C separately has open defect
      P1-8: computed on 5 donors via the deprecated label-based selection.)
- [ ] **`figureS_sensitivity`** — reads `niche_index_per_donor`; the stability
      arm is now a 2-component composite over 89 donors instead of 1 over 46.
- [ ] **`figureS_balance_by_state`** — same script as `figure_mechanism_main`.

### Tier 3 — regenerate for consistency; content likely unchanged

- [ ] **`figureS_bm_associations`** (S17) — already regenerated with panel E on
      2026-09-02, but *before* the continuum re-root. Re-run to pick up
      `bm_continuum_summary`.
- [ ] **`figureS_receiver_robustness`** — reads `bm_vs_*`.
- [ ] **`figure_disease_main`**, **`figureS_disease_robustness`** — read
      `agtr1_copd_*`, re-exported with corrected signs.
- [ ] Emitted as side effects of the above scripts, no known content change:
      `figureS_acta2_control`, `figureS_alluvial`, `figureS_crossspecies_mouse`,
      `figureS_program_category`, `figure_ccc_nichenet`, `figureS_bm_copd`,
      `figureS_ras_landscape`.

### Not affected
`figureS_state_annotation`, `figureS_agtr1_dropout`, `figureS_cogaps_validation`,
`figureS_continuum_stability`, `figureS_nichenet_specificity` — their inputs did
not change. (`figureS_continuum_stability` reads the sensitivity sweep, which was
**not** re-run under the new root — see Open questions.)

---

## How to run

```bash
cd figures/_m && sbatch ../_h/step_figures.sh
```

One job runs all twelve scripts in dependency order, ending with
`assemble_mechanism_figures.R` (the manifest step). Everything upstream is
already re-run, so no analysis step is needed first.

---

## After regeneration — the manual layer

These are hand-curated and do **not** come out of `step_figures.sh`:

- [ ] `figures/pericyte-programs.02.svg` — Inkscape re-assembly (md5 currently
      `50541edf`, differs from the generated `figure_pericyte_layer.svg`).
- [ ] `figures/basement-membrane.03.svg` — currently md5-identical to the
      generated file, so a straight copy suffices *if* it stays that way.
- [ ] `figures/supplementary/*.SNN.{svg,png}` — the curated export set. It stops
      at **S13** while `mechanism/README.md` references up to **S17**, and
      `acta2-control.S05.svg` is from 2026-07-28. S01/S07/S08 also disagree with
      the manifest in `assemble_mechanism_figures.R`.
      **Blocked on a decision: is the manifest or the folder authoritative?**

---

## Verification checklist

- [ ] `figure_pericyte_layer` panel F: `basement_membrane_score` ρ is now
      **positive**, `vascular_stabilizing_score` **negative**.
- [ ] `figureS_state_composition`: shows 4 disease groups, per-group N visible,
      COPD marked non-estimable.
- [ ] `figure_basement_membrane` panel G: `rho_bm` no longer starred.
- [ ] `figureS_sensitivity`: stability arm labelled as 2-component.
- [ ] Legends in `mechanism/README.md` match the redrawn panels — S11 and S17
      legends were rewritten on 2026-09-02 and describe the *new* numbers.

---

## Open questions that affect what the figures should say

1. **Is this still an "injury continuum"?** Rooted correctly, the injury
   programs *fall* along the axis (inflammatory −0.452, activated/migratory
   −0.392) and basement membrane is the only program that *rises* (+0.322). But
   **every** program except BM falls, which is equally consistent with an overall
   score-magnitude gradient. Any axis label asserting "injury" is currently
   unsupported. **Resolve before finalising Figure 2F's axis title.**
2. **The continuum sensitivity sweep (`02b`) has not been re-run under the new
   root**, so `figureS_continuum_stability` reports stability across 8 roots that
   were all evaluated against the old rooting logic.
3. **`niche_stability_score` is now a real 2-component composite** — any legend
   describing it as vascular-stabilizing alone is stale.
