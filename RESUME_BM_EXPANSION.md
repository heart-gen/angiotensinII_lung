# Resume notes — basement-membrane + AGT-axis expansion

Paused 2026-07-21. **Nothing has been committed.** Two background submitters are
still running and will keep feeding SLURM jobs as QOS slots free (10-job cap):

```
scratchpad/cascade2.sh        # cell_communication chain + pseudobulk refreshes + figures
scratchpad/cross_species.sh   # mouse scVI re-run with the BM panel
```
Scratchpad lives under
`/ocean/projects/bio250020p/shared/tmp/claude-96597/.../scratchpad/`.
To stop them: `pkill -f cascade2.sh; pkill -f cross_species.sh` (this does NOT
cancel already-submitted jobs; use `scancel` for those).

All jobs now run on account **bio250020p**. New module scripts carry it in the
`#SBATCH` line; older modules still say `bio260021p` on disk and are overridden at
submission with `sbatch --account=bio250020p`. Decide whether to make that
permanent in their tracked files.

## Is it safe to exit the session?

Yes. Verified at pause time:

- **SLURM jobs are fully independent of the session.** Once submitted they run to
  completion regardless of whether Claude Code, the terminal, or the SSH
  connection is alive.
- **Both background submitters have PPID 1** (already reparented to init) and were
  started with `nohup`, so they survive the session ending and will keep feeding
  jobs as QOS slots free. They are not guaranteed to survive a login-node reboot.
- **All work is written to disk.** Nothing is held in session state.
- **Nothing is committed to git.** The work is safe on disk but unversioned, so
  avoid destructive git operations (`git checkout .`, `git clean -fd`, `git reset
  --hard`) in this repo until it is committed.

If the submitters die before finishing, run `bash resubmit_remaining.sh` — it
submits everything still outstanding, in dependency order, and lists the job IDs
already submitted so you do not duplicate them. Check `squeue -u $USER` first;
the script does not deduplicate.

To stop everything instead:
`pkill -f cascade2.sh; pkill -f cross_species.sh; scancel -u $USER`

## State of the world

### Done and verified
- `pericyte_states` re-run with the 6th (basement-membrane) panel. **All 11,680
  cells kept identical cluster assignments**; same chosen solution
  (`leiden_n30_r0.5`, 6 clusters, median Jaccard 0.966). Only the label changed:
  clusters 1/3/5 `fibroblast_like` → `basement_membrane` (4,200 cells, 36%).
  5-panel backup preserved at `pericyte_states/_m_backup_5panel/`.
- All five original program **scores are bit-identical** before/after
  (max |diff| = 0), so score-based analyses are numerically unaffected.
- `basement_membrane/` steps 0,1,2,3,5 complete. `agt_axis/` step 0 complete.
- `niche_index`, `sensitivity` re-ran clean. `pathway_balance` re-submitted after
  the injury-grouping fix.
- Figures generated and visually checked: `figure_basement_membrane`,
  `figureS_ras_landscape`, `figureS_bm_copd`.
- All 18 edited Python/R scripts parse cleanly.

### Running / queued at pause
| Job | Status | Blocks |
|---|---|---|
| `bm_nichenet` | **DONE** — summarized | — |
| `agt_vs_ligands` (300-draw rank bootstrap) | running ~40 min | AGT summary section |
| `cogaps_run` (nP×seed sweep) → `cogaps_select` (de-novo rank) → `cogaps_validate` | pipeline renumbered 00–05; sweep pending submit | set chosen nP in `step_3`/`step_5` from `cogaps_nP_selection.tsv` |
| `ccc_prepare` | **DONE** — re-ran on `RM-shared` 2026-07-21 (EM no longer needed); `ccc_niche.h5ad` refreshed 2026-07-22 | downstream cell_communication chain ran; `figures` refreshed |
| `cross_species` step_1→3 | **DONE** — see the cross-species section below | — |

`ccc_prepare` waiting on EM resources is the long pole. If it stays queued,
consider whether `cell_communication` truly needs a full re-run: only
`ccc_group_state` / `pericyte_program` depend on `state_program`; `ccc_group`
(used by the BM selectivity result) does not, so the completed selectivity numbers
are already valid.

If the AGT bootstrap is too slow, drop `--nboot` from 300 to ~100 in
`agt_axis/_h/step_1.sh`. The point estimate (AGT rank 11) is already established;
the bootstrap only adds the interval.

## Open decisions — ALL RESOLVED 2026-07-22

Resolutions are recorded below each item. Two follow-up jobs were still running at
the time of writing: `cogaps_run`→`validate`/`project`/`stability` (42529933–36,
the BM-panel re-run) and `bm_nichenet` at `--nperm 10000` (42529989).

1. **Account** — ~~make `bio250020p` permanent?~~ **DONE.** All 28 remaining
   `step_*.sh` files were rewritten; 41 of 41 now carry `--account=bio250020p`, so
   no submission needs a command-line override.
2. **Mouse `03.conservation_stats.R` is now superseded in part** — the audit shows
   its injury-vs-stabilizing contrast (n = 1,144, `score_median_split`) is 96%
   smooth muscle and its continuous `Agtr1a ~ vascular_stabilizing_score` result is
   a cell-type contrast on denoised data. Decide whether to (a) restrict that
   script to the 41 pericytes, (b) relabel its outputs as mural/SMC-level and drop
   them from the manuscript, or (c) retire it in favour of `04.species_comparability.py`.
   **RESOLVED — option (c), retired.** The script now carries a `RETIRED — DO NOT
   RUN, DO NOT CITE` header stating all three reasons, `step_3.sh` no longer
   invokes it (it logs the retirement instead), and the run guide in
   `MECHANISM_ANALYSES.md` was updated. Script and outputs are kept for provenance.
3. **Downstream prose** — ~~"fibroblast-like pericyte program" needs rewriting.~~
   **DONE** in `MECHANISM_ANALYSES.md`: the panel list is now six not five, the
   cluster 1/3/5 program is basement-membrane in the Methods, Results, ACTA2
   control, three-lens paragraph, CCC receiver name and alluvial legend. Uses of
   `fibroblast_like` that refer to the *fibrillar-ECM panel or score* (continuum
   trends, BM-vs-fibrillar orthogonality, injury composite components, CoGAPS
   pattern 3) were deliberately left alone — those are still that panel.
4. **`pericyte_states/_m/stats_data/` labels were stale** — **FIXED.** `step_1`
   re-ran (42529924) and now writes `composition_program_basement_membrane_*`.
   Verified against `_m_backup_5panel/` first: cluster assignment agrees
   11,680/11,680 and all five original scores are bit-identical (max |diff| = 0),
   so this was a labelling fix, not a numerical one. 30 superseded/orphaned files
   were quarantined to `_m/stats_data/_superseded/` with a README.
   **One thing needs a decision there:** the current `01.state_stats.R` contains no
   per-program score-by-disease analysis at all, yet `MECHANISM_ANALYSES.md` quotes
   "per-program state scores are likewise flat (all BH ≥ 0.79)". The nearest
   surviving file reports BH ≥ 0.98, so the quoted number predates every output on
   disk. Either restore that analysis to the script and re-run, or drop the
   sentence — do not patch the number from the quarantined files.
5. **`pericyte_cogaps` steps 3 and 4 were stale, and the BM panel was missing** —
   **FIXED, re-running.** `00.prepare_cogaps_input.py` now forces in the
   basement-membrane panel; the feature space grew from 2,030 to **2,038 genes**
   (2,000 HVG + 38 forced), so all 13 BM genes are now present instead of 5. The
   whole chain was resubmitted with dependencies (42529932–36), which also
   regenerated the stale projection and stability outputs. **All CoGAPS text has
   been refreshed** — `pericyte_cogaps/_h/ANALYSIS_SUMMARY.md` and the three
   CoGAPS passages in `MECHANISM_ANALYSES.md` now carry the 2,038-gene numbers.
   What changed: the mural programs **merge** into one pattern (contractile
   ρ = 0.61, stabilizing ρ = 0.59, the only two panels enriched by marker overlap);
   **basement membrane now resolves into two positively-correlated patterns**
   (ρ = 0.45 and 0.28, both pericyte-first in the projection) instead of one
   uniformly-negative one; the disease-associated pattern is the fibrillar-matrix
   factor (BH = 5.8e-3), which is fibroblast-dominant, not pericyte-specific;
   inflammation still has no marker-supported pattern. Seed stability improved for
   the outlier (0.30 → 0.700) but that outlier **is now the basement-membrane
   pattern**, which is a caveat that must travel with it.
   **A latent trap was found and is being fixed.** `05.project_niche.R` labels its
   outputs from `cell_communication/_m/cogaps_receiver_annotation_np5.tsv`, which
   had not been regenerated — so after the re-run all five hardcoded labels in
   `projected_*_np5.tsv` and `injury_pattern_disease_np5.tsv` named the *wrong*
   biology (e.g. `Pattern_4 (vascular_stabilizing)` is in fact the fibrillar-matrix
   disease pattern). `cell_communication/step_5.sh` (42538664) regenerates the
   annotation from the current correlations and cogaps `step_3` (42538665) re-runs
   the projection behind it. **Do not quote per-cell-type disease contrasts until
   those land**; the manuscript paragraph carries a blockquote saying so.
6. **`pathway_balance` statistics had never successfully run** — **FIXED.**
   `01.balance_stats.R:94` used `df[, inj_cols]` on a data.table, which resolves
   the bare symbol as a column name; both 2026-07-21 attempts died there, leaving
   the disease-level tables at their 2026-06-16 values from the *old label-based*
   injury selection. Fixed to `..inj_cols`; the job now completes (42529923). A
   hardcoded `cat()` asserting that the disease effect "collapses" under
   injury adjustment was also removed — under the continuous selection it does not.
   **The conclusions changed and `MECHANISM_ANALYSES.md` was rewritten accordingly**
   (by-disease contrasts now n.s.; the adjustment no longer absorbs the disease
   term; the AT1R-arm-only result stands).

## Summaries (per the per-analysis-summarization skill)

- `basement_membrane/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** Covers selectivity,
  the gate, the BM-vs-fibrillar axis, the continuum, COPD, and the BM NicheNet
  target ranking.
- `agt_axis/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** The bootstrap CI (rank median
  18, 95% interval 6–63), ligand co-expression and target-convergence sections
  were added after `agt_vs_ligands` finished; numbers verified against
  `agt_ligand_rank_bootstrap.tsv` and `agt_target_overlap.tsv`.
- `cross_species/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** Covers the comparability
  audit and the cell-level *Agtr1a* compartment result, and states explicitly
  that `03.conservation_stats.R` outputs should not be cited.
- `pericyte_states/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** State discovery,
  continuum + sensitivity, the three-lens AGTR1 negative result, the ACTA2
  control, and the null disease-composition results.
- `pericyte_cogaps/_h/ANALYSIS_SUMMARY.md` — **COMPLETE**, with two caveats
  recorded in it (see open decisions 4 and 5 below).
- `cell_communication/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** Niche construction,
  four receiver schemes, the NicheNet ranking with its permutation control, and
  the 105-donor validation.
- `niche_index/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** Records that the
  injury-stromal *arm*, not the net index, is the interpretable readout.
- `pathway_balance/_h/ANALYSIS_SUMMARY.md` — **COMPLETE**, on the corrected
  2026-07-22 numbers.
- `sensitivity/_h/ANALYSIS_SUMMARY.md` — **COMPLETE.** Includes the exact
  leave-one-study-out breakdown (13 of 16 significant) and the two metadata
  limitations (no smoking status for any diseased donor; no medication data at all).

**All nine mechanism modules now have a per-analysis summary.**
`disease_association/` and `localization/` do not — they predate this expansion.

### Follow-up flagged by the NicheNet result — ADDRESSED, running
`07.bm_nichenet_targets.R` ran with `--nperm 200`, which floors the empirical
p-value at 1/201 and means **no ligand can survive BH correction across several
hundred candidates** (all `perm_p_BH` ≈ 0.18). The z-scores were interpretable
(MMP14 z = 87, TIMP2 z = 37.5, TGFB2 z = 5.9, TGFB1 z = 2.4) but the claim was not
publishable at that resolution.

Re-ran 2026-07-22 at `--nperm 10000` (42529989, **COMPLETED** in 2 h 48 m). The
serial `vapply` loop would not fit in any reasonable walltime at 50× the work, so
it was parallelised with `parallel::mclapply` over `SLURM_CPUS_PER_TASK` using
`RNGkind("L'Ecuyer-CMRG")` + `mc.set.seed` (still deterministic given SEED and
core count), with an explicit check that no permutation silently returned an error
object. `step_4.sh` now requests 16 cores and 24 h.

**Result: exactly one ligand of 321 survives BH < 0.05 — `MMP14`**
(p = 9.999e-5, BH = 0.032). Everything else, including TIMP2 and both TGF-β
ligands, sits at BH ~= 0.189. **The 200-perm z-scores were badly inflated** by the
underestimated null SD and must not be quoted anywhere: MMP14 87.0 -> 27.1,
TIMP2 37.5 -> 14.4, TGFB2 5.9 -> 3.2, TGFB1 2.4 -> 1.7. TGFB2 stays 5th and TGFB1
12th; TGFB3 274th. `basement_membrane/_h/ANALYSIS_SUMMARY.md` was updated.
Still to do: `MECHANISM_ANALYSES.md` does not yet state the MMP14 result.

## Cross-species result (mouse) — RESOLVED, section rewritten

`mouse_integrated.h5ad` was rebuilt and retains **13/13 BM orthologs** (previously
1/13). All six mouse panels retain 100% of their genes, so ortholog attrition is
NOT an issue.

A comparability audit (`cross_species/_h/04.species_comparability.py`, new; run via
`step_4.sh`) found that **mouse and human are not comparable at the state level**,
for reasons that predate the BM panel:

1. **The mouse "mural" set is 96.4% smooth muscle** — 1,016 vSMC + 87 PA-SMC +
   **41 pericytes**. Human is 11,680 pericytes. The mouse `basement_membrane`
   "state" is 84 of the 87 PA-SMCs and contains **zero pericytes**;
   `fibroblast_like` is 52 vSMC, also zero pericytes.
2. **Cell type is aliased with dataset** — 2 datasets are vSMC-only, 2 are
   pericyte+PA-SMC-only. Only pericyte-vs-PA-SMC is estimable within dataset.
3. **Different labelling rules** — human uses relative enrichment (z-score →
   cluster mean → argmax); mouse used raw magnitude, which the human code
   explicitly rejects. Applying the human rule moves BM 89→5, inflammatory 8→94.
   Mouse also clusters the whole lung, so 87% of mural cells sit in one cluster.
4. **`scvi_corrected` is dense** — detection of 1.000 is vacuous, and the
   rho = 0.68-0.90 "conservation" correlations mostly re-express cell-type identity.

### What survives — and it is solid
Cell-level inspection of all 41 pericytes on the raw `lognorm`/`counts` layers:

- **31/41 pericytes (75.6%) carry >=1 `Agtr1a` UMI**, across **15/18 donors** and
  both datasets; **0/87 PA-SMC** carry a single UMI.
  Fisher p = 1.4e-9 and 3.0e-13 per dataset; Mantel-Haenszel stratified p = 8.4e-20.
- The 24% negatives are **depth dropout**: datasets differ ~360x in depth
  (median 291 vs 104,973 UMI/pericyte); detection tracks depth (rho = 0.43,
  p = 0.005); deep dataset detection is 15/16 with median 928 UMI/positive cell.
- **`Agtr1b` and `Agtr2` are absent** from mouse pericytes (0/41, max 0 UMI), so
  `Agtr1a` is the sole AngII receptor there.
- **All 41 pericytes are from `normal` animals** (only vSMC have influenza), so the
  injury-loss/losartan-rescue phenotype is **not testable** in this dataset and must
  be cited to primary literature.

`MECHANISM_ANALYSES.md` Methods + Results were rewritten accordingly; the
NEEDS REVISION blockquote is gone. New outputs are the seven
`stats_data/species_comparability_*.tsv` tables.

## Headline results so far

- Pericytes rank **#1 of 22 cell types** for COL4A1, COL4A2, COL18A1, LAMB1, NID1,
  NID2 — but not the laminin α-chains. Pre-specified negative controls behaved
  (LAMA3 17th, LAMA5 9th). LAMA4, the textbook "pericyte laminin", is higher in
  alveolar fibroblasts.
- On BM-minus-fibrillar, pericytes separate from every fibroblast population
  (P down to 5×10⁻³⁰¹) but are **indistinguishable from capillary endothelium**.
- The BM cluster pattern **survives orthogonalization** against fibrillar ECM.
- Along the injury continuum, fibrillar rises (median ρ = 0.040, BH = 0.014) while
  BM does not (ρ = 0.025, P = 0.84); the BM-minus-fibrillar switch index falls
  (ρ = −0.054, BH = 0.014).
- COPD: null on LAMB1/LAMA4 across 7 powered compartments (all BH ≥ 0.81), with a
  working IPF positive control (endothelial HSPG2 P = 1.5×10⁻⁵). **No pericyte
  COPD claim** — only one Control donor clears 5 pericytes.
- **Zero of 22 cell types** carry a complete AGT→AngII→AT1R circuit, and none
  carries more than one of the three requirements. Renin effectively absent
  (max 0.023).
