# Figure regeneration after the P1 fixes (opened 2026-09-02)

> **Superseded in part, 2026-09-10 — the disease analyses left this repository.**
> `disease_association/`, `sensitivity/`, `niche_index/_h/01.niche_disease_stats.R`
> and the disease figures and tables moved to **heart-gen/lung-pericyte-analysis**.
> Entries below that name `figure_disease_main`, `figureS_sensitivity`,
> `figureS_disease_robustness`, S12, S16, Tables S13A/S13E/S13F or S14 refer to
> work that is no longer built here. They are left in place because this is a
> record of what was done, not a list of what to do; do not act on them here.


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

## Status: ✅ CLOSED 2026-09-07

All three tiers regenerated and the verification checklist run and passed (below).
The tier checkboxes are left as written for the historical record; the
authoritative statement is the checklist.

Two things came out of the verification that were not on any list:

1. **`figureS_state_composition` was drawing the COPD donor as an unlabelled
   "NA" column**, in two places. Found by the checklist, fixed, re-run.
2. **One checklist item was itself stale** and expected the pre-P1-10 LOSO
   result. Rewritten.

Still open, and **not** part of P1-21: the `figures/supplementary/` curated set
(see the manual layer below) needs a decision, not a re-run.

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

## Verification checklist — ✅ RUN AND PASSED 2026-09-07

All figures regenerated by job **45458692** (16:52; earlier passes 45449432 and
45458501). Verified with `pdftotext` against the panels themselves where the
label is extractable, and against the source table where the check is a sign or
a significance the script derives.

- [x] `figure_pericyte_layer` panel F: `basement_membrane_score` ρ **positive**
      (+0.138 donor), `vascular_stabilizing_score` **negative** (−0.573).
      ✅ `pseudotime_trend_correlations.tsv`, figure redrawn from it.
- [x] `figureS_state_composition`: 4 disease groups, per-group N in the axis
      labels, COPD flagged. ✅ **This one failed on the first pass and was fixed.**
      `GRP_LEVELS` omitted COPD, so the donor P1-2 restored was drawn as an
      unlabelled **"NA"** column; panel D's legend carried a second NA series from
      the same cause. `dx()` now *errors* on an unmatched group rather than
      drawing it, axis labels read `Healthy\n(n = 42)` … `COPD\n(n = 1)\nnot
      estimable`, and panel D drops the 10 Healthy-vs-COPD contrasts with the
      reason in its axis title. Zero NA columns remain.
- [x] `figure_basement_membrane` panel G: `rho_bm` (BH **0.715**) and
      `rho_switch` (BH **0.0613**) both render **n.s.** ✅ The starring rule at
      `basement_membrane_figure.R:355–359` reads `p_BH` directly; 3 `n.s.` labels
      present in the PDF.
- [x] `figureS_sensitivity`: stability arm is 2-component
      (`n_niche_stability_score_components` = 2 in **89/89** donors after the
      airspace export fix).
- [x] `figure_disease_main` panel C: x-axis reads **"Disease-attributable AGTR1
      variance (Δ marginal *R*²)"**, and the fibroblast block runs
      **Peribronchial (P = 0.08) → Adventitial (0.10) → Alveolar (0.14)**. ✅
- [x] `figureS_sensitivity` panel D: x-axis reads **"Fibrotic/ILD effect (dataset
      left out)"**. ✅
      ⚠️ **The second half of this item was itself stale and has been rewritten.**
      It expected *"exactly one point coloured `p < 0.05` in the injury-stromal
      facet and one in the `AGTR1_pos_frac` facet"*. That described the 47-donor
      fit. On the corrected 89-donor fit the correct expectation is
      **23 of 23 significant in the injury-stromal facet, 22 of 23 in
      `niche_index`, and 0 of 23 in both `injury_frac` and `AGTR1_pos_frac`**.
- [x] Legends in `mechanism/README.md` match the redrawn panels. ✅ Verified: S12
      panel D reads "significant in all 23 refits, estimates 0.643–0.920" and
      panel B reads "the stability arm is **flat** (*F*₃,₅₂.₃ = 0.15, *P* = 0.93)"
      — both current.
      ⚠️ **One legend edit remains:** the S12 legend still opens with
      *"LEGEND UPDATED 2026-09-07 (P1-10); **THE FIGURE HAS NOT BEEN REDRAWN**…
      Panels B, C and D on disk still plot the old one"*. It has been redrawn
      three times since; that warning block should be deleted.

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

---

## Added 2026-09-07 — the niche_index / sensitivity model change (P1-10)

`niche_index/_h/01.niche_disease_stats.R` and `sensitivity/_h/00.sensitivity.R`
moved from `lm(~ disease_group + age + sex)` (no study term, 47 of 89 donors) to
`lmer(~ disease_group + sex + (1 | study))` on all 89. Both modules were re-run
(jobs 45444388 / 45444507 / 45444589). **Every figure reading their tables is
stale**, and three panels now plot something that is no longer true:

| Figure / panel | What is on disk | What the tables now say |
| -------------- | --------------- | ----------------------- |
| **S12 B** — stability arm | drawn as **rising** in fibrotic/ILD | **flat**, *F*₃,₅₂.₃ = 0.15, *P* = 0.93 |
| **S12 C** — smoking availability | 14 of 32 Healthy, 0 of 37 diseased | **21 of 42**, 0 of 47 |
| **S12 D** — LOSO | **17** refits, **1** significant, est 0.440–0.548 | **23** refits, **23** significant, est 0.643–0.920 |
| **S12 A / E / F** | 36-donor `+ age` marginal means | 89-donor primary means |

Two things to verify after the redraw, because they are new failure modes:

1. **The COPD arm is one donor.** Dropping `+ age` restored it. Panels that draw
   disease groups must either omit it or mark it — and any *P* quoted next to it
   should be the `p_excl_small_groups` column, which is now emitted alongside
   `p.value` in every `niche_index` posthoc/anova table.
2. **`leave_one_study_out.tsv` gained columns** (`n_studies`, `singular`, `model`)
   and grew from 68 to 92 rows. `sensitivity_robustness_figure.R` should be checked
   for anything that assumed the old width or the old 17-row-per-response shape.

---

## Added 2026-09-07 (second batch) — P1-6, P1-8, P1-19

### `figure_mechanism_main` panel C — the donor set changes from 5 to 59

The panel used a label-based injury selection that `pathway_balance` abandoned
(P1-8). It now reads `pathway_balance/_m/stats_data/balance_donor_injury_selected.tsv`
and **stops** if that file is missing or if the donor count comes back under 20.

| | On disk | After regeneration |
| - | --- | --- |
| donors | **5** (4 Healthy, 1 IPF) | **59** (26 H / 19 F / 14 O) |
| Wilcoxon bracket | computed on 4 vs 1 | 26 vs 19 |

The underlying statistics also moved when the module's `+ age` filter was dropped:
Healthy vs Fibrotic/ILD balance +0.145 (*P* = 0.0038) and the **AT1R arm** +0.186
(*P* = 9.6 × 10⁻⁴), against AT2R flat. Verify the panel's message matches its
"(corollary)" label — the disease term now **does** collapse under
injury-adjustment (0.145 → 0.064, *P* = 0.194), which is what the corollary framing
claims.

### Figure S12 panel D legend — no wording change needed, but check the axis

Already covered in the first batch. No new action.

### Panel D of `figure_disease_main` — legend wording (P1-19)

Drop "replication" for "independent evaluation". Replacement text is drafted in
`writings/pi_briefings/README_AUDIT.md` Defect 11. Two facts the legend must carry
and currently does not:

- **The COPD arm has no HLCA comparator at all.** `disease_association/_h/05`
  excludes COPD, and its "Other" group is COVID-dominated. Half the panel is
  therefore not a cross-cohort comparison of any kind.
- Neither HLCA estimate is significant, so the honest phrase is *"no directional
  finding to replicate"*, never *"the replication failed"*.

New supporting file: `agtr1_cross_cohort_signs.tsv`.

### NicheNet ligand ordering (P1-6) — any figure that ranks ligands

`nichenet_specificity_Pericytes.tsv` now carries `rank_by_z`, `rank_by_aupr`,
`p_emp_at_floor` and `rank_note`. **Order by `z`.** The two orderings disagree at
the top (TGFB1 leads on *z*, TGFB2 on AUPR), and the permutation *p* cannot break
the tie because it is floored.

Permutations were raised 1,000 → 10,000, so **every *z* and every empirical *p* in
the current figures is from the old run** and must be refreshed. Check
`nichenet_specificity_figure.R` for anything that plots `p_emp` or assumes the old
floor.

### `donor_validation_results.tsv` gained rows and a column

Two new predictors (`receiver_TGFB_SMAD`, `receiver_TGFB_IEG`) and a `role`
column. `donor_validation_scatter.pdf` is regenerated. If a supplementary panel
shows the donor-validation predictors as a fixed set, it will need the two new rows
— and the IEG row must be drawn as a **control**, not as a result.

---

## Done 2026-09-10 — `figure_pericyte_layer` panel D rebuilt on the count-model arbiter

Panel D drew three lenses and *cited* the count model in its legend. It now draws
four, and all four are fit at one unit.

**What was wrong.** The denoised series came from `03.agtr1_lenses.R`, a CELL-level
`lmer(y ~ state_program + (1|donor_id))` on 11,680 cells — no study term, no depth
covariate, no propagation of imputation uncertainty. The count-model arbiter it was
being compared against is 214 donor × cluster pseudobulks with `(1|study) +
(1|donor_id)` and a library-size offset. On a common `log(AGTR1 per 10⁴)` scale
that is SE 0.061–0.076 against 0.131–0.204 — a 2.0–2.7× gap that reads as the
denoiser being the more precise measurement. It is not: refit at the shared unit
the same denoised lens gives 0.110–0.221, within 5–43 % of the count model.
The denoised 95 % interval on basement-membrane − vascular-stabilizing
(1.112 [1.078, 1.147]) also **excludes** the no-imputation estimate (1.255), which
is what over-tight bars look like from the outside.

**What changed.**

- New producer `basement_membrane/_h/19.agtr1_lens_by_cluster.R` (+ `step_10.sh`,
  runs after `step_6.sh`). It refits raw / detection / denoised at the count
  model's own pseudobulk unit and hard-stops if the unit count ever stops matching
  `agtr1_count_by_cluster.tsv`.
- It also **replaces an orphan**: `agtr1_lens_by_cluster_emmeans.tsv` was on disk
  with no producing script anywhere in the repo and carried a superseded
  `denoised OLD (invalid)` series from the scVI model that failed its validity
  gate. An orphan table must not feed a main figure.
- Panel D's x-axis moved from `state_program` to the Leiden subclusters P0–P5.
  The subclusters exclude AGTR1 from the 2,000 HVGs that define them, so the
  grouping is independent of the readout; `state_program` is a marker-panel argmax.
  It is also where the evidence is — at program level only 1 of 5 count-model
  specs separates BM from VS, at subcluster level all twelve pseudobulk BM-vs-VS
  contrasts agree in sign with 10 of 12 significant.
- New source-data parts S05D1–S05D4. **S05D3/S05D4 are the count model's first
  appearance in the supplement at all** — before this it was quoted in legends and
  never tabulated. S05C1–C3 stay, because `figureS_acta2_control` panel A is still
  keyed on programs.

**Still open — the same defect in `figureS_acta2_control` panel A (S5).** That
panel draws the denoised series from the cell-level program fits and the count
series from pseudobulk, side by side on one axis, exactly the mismatch fixed here.
Closing it needs a program-level analogue of `19.agtr1_lens_by_cluster.R`
(`~ state_program + depth + (1|study) + (1|donor)` on the 154 donor × program
pseudobulks). The panel's conclusion is not at risk — ACTA2 and AGTR1 are still not
co-ordered — but the relative bar widths in it are not currently meaningful.
