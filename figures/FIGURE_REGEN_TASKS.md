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
| **P1-16** (2026-09-07) effect size | `agtr1_celltype_disease_{omnibus,ranking}.tsv` | Panel C's statistic **changed**: partial η² → Δ marginal *R*². Ordering survives (ρ = 0.893) but the **top cell type moves** alveolar → peribronchial fibroblasts. |
| **P1-15** (2026-09-07) min-cells | `mincells_sensitivity.tsv` | The ≥5 rung is real for the first time: +0.526, *P* = 0.033 on 79 donors, against a duplicated +0.728 before. |
| **P1-9** (2026-09-07) study guard | `balance_disease_injury_adjusted.tsv`, `balance_arm_decomposition.tsv` | Each row is now **duplicated** across a `study_guard` column. Any figure reading these must filter to the PRIMARY rows or it will double-plot. |
| **P1-11** (2026-09-07) LOSO caption | *(no table moved)* | Axis label only: "study left out" → "dataset left out". |

---

## Figures to regenerate

Ordered by how badly the current version misleads.

### Tier 1 — currently plots an inverted direction

- [ ] **`figure_pericyte_layer`** (Figure 2) — panel F is the continuum lollipop.
      Every feature's ρ has flipped sign. **The panel currently shows the
      continuum running the wrong way.**
      *Also check:* the panel-F guard added for `AGTR1_scvi` still fires.
- [ ] **`figureS_pericyte_layer`** — same script, same continuum inputs.
- [ ] **`figureS_continuum_stability`** (S7) — the `02b` sweep was **re-run under
      the corrected root on 2026-09-07**, so every ρ in panels A–C has flipped and
      the input now carries two extra roots (the `basement_membrane` and
      `activated_migratory` programs) and an eighth feature (`AGTR1_scvi`). The
      script was taught both on 2026-09-07; the figure itself has not been redrawn.
      *Also check:* `AGTR1_scvi` lands on OKABE position 8 (grey) — the denoised
      lens is the readout, so grey is probably the wrong colour for it.
- [ ] **`figure_basement_membrane`** — panel G is per-donor continuum ρ by
      metric. `rho_switch` is now **n.s.** (−0.055 at BH = 0.014 → +0.024 at
      BH = 0.061); the tracer (−0.196 → +0.203) and TGF-β (+0.133 → −0.131)
      invert. `rho_bm` is a null (−0.014, BH = 0.715) but **was already one
      before the re-root** (−0.021, BH = 0.624) — do not annotate it as a change.
      Panels that annotate significance must be re-checked, not just redrawn.

### Tier 2 — data moved substantially

- [ ] **`figure_disease_main`** — **panel C's statistic changed** (P1-16,
      2026-09-07). The script now plots `delta_r2_marginal` and **`stop()`s if
      that column is missing**, so a stale
      `agtr1_celltype_disease_ranking.tsv` fails the job loudly rather than
      falling back. The x-axis label changed to "Δ marginal *R*²" and the row
      order changes: **peribronchial fibroblasts move to the top**, ahead of
      alveolar. Moved up from Tier 3, where it sat only for the `agtr1_copd_*`
      sign corrections.
      *Also check:* the panel's *P* labels are unchanged (the omnibus test did
      not move), so a redraw that changes the labels means something else broke.

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
      **Panel D's x-axis label also changed** (P1-11, 2026-09-07): the refits drop
      *datasets*, not studies. The script now prints refits / positive /
      significant / estimate range / most-influential drop for all four
      responses — **write the caption from that block**, which is what the
      13-of-16 defect existed for.
- [ ] **`figureS_balance_by_state`** — same script as `figure_mechanism_main`.

### Tier 3 — regenerate for consistency; content likely unchanged

- [ ] **`figureS_bm_associations`** (S17) — already regenerated with panel E on
      2026-09-02, but *before* the continuum re-root. Re-run to pick up
      `bm_continuum_summary`.
- [ ] **`figureS_receiver_robustness`** — reads `bm_vs_*`.
- [ ] **`figureS_disease_robustness`** — reads `agtr1_copd_*`, re-exported with
      corrected signs.
- [ ] Emitted as side effects of the above scripts, no known content change:
      `figureS_acta2_control`, `figureS_alluvial`, `figureS_crossspecies_mouse`,
      `figureS_program_category`, `figure_ccc_nichenet`, `figureS_bm_copd`,
      `figureS_ras_landscape`.

### Not affected
`figureS_state_annotation`, `figureS_agtr1_dropout`, `figureS_cogaps_validation`,
`figureS_nichenet_specificity` — their inputs did not change.

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
- [ ] `figure_disease_main` panel C: x-axis reads **Δ marginal *R*²**, and the
      top row of the fibroblast block is **peribronchial**, not alveolar.
- [ ] `figureS_sensitivity` panel D: x-axis reads "dataset left out"; exactly
      **one** point is coloured `p < 0.05` in the injury-stromal facet and one in
      the `AGTR1_pos_frac` facet (both are the `Lafyatis_2019` refit in the
      latter case).
- [ ] Legends in `mechanism/README.md` match the redrawn panels — S11 and S17
      legends were rewritten on 2026-09-02 and describe the *new* numbers; the
      S12 panel-D, S16 min-cells and main-figure panel-C legends were rewritten
      on 2026-09-07.

---

## Open questions that affect what the figures should say

1. ~~**Is this still an "injury continuum"?**~~ **RESOLVED 2026-09-07 — it is not,
   and the name has been dropped repo-wide.** Rooted correctly, the injury programs
   *fall* along the axis (inflammatory −0.452, activated/migratory −0.392) and
   basement membrane is the only program that *rises* (+0.322); the re-run `02b`
   sweep (question 2) then showed the poles are basement-membrane versus everything
   else, with `activated_migratory` on the *same* end as `vascular_stabilizing`.
   **Axis titles should read "vascular-stabilizing ↔ basement-membrane", never
   "injury".** Applies to Figure 2E/2F, `figure_mechanism_main` panel E, and S7.

   Still open, and unaffected: **every** program except BM falls, which is equally
   consistent with an overall score-magnitude gradient. Depth-adjusted partials
   preserve the pattern (BM +0.400, inflammatory −0.511), so it is not a pure depth
   artifact — but "BM rises while the rest fall" is the claim to make, not a
   mechanism.
2. ~~**The continuum sensitivity sweep (`02b`) has not been re-run under the new
   root.**~~ **RESOLVED 2026-09-07** — re-run, and it bears directly on question 1
   above. The axis is bipolar with **`vascular_stabilizing` and
   `activated_migratory` at the SAME end** and `basement_membrane` at the other:
   rooting at `activated_migratory` reproduces the canonical signs, rooting at
   `basement_membrane` reverses all 8 features together. So the axis is not
   stabilizing ↔ injury; the injury pole is not a pole. Parameter robustness is
   clean (sign consistency 1.00 for all 8 features over 18 settings, sd ρ ≤ 0.035)
   and |ρ| survives every re-rooting (Spearman of |ρ| ≥ 0.88, = 1.00 for the
   cluster-0 and PC1-max roots).
3. **`niche_stability_score` is now a real 2-component composite** — any legend
   describing it as vascular-stabilizing alone is stale.
