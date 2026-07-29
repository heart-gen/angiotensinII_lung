# Mechanism figures — index and figure legends

Publication-ready figures for the mechanistic revision (target: *Circulation Research*).
All panels are exported as `.pdf` (vector, cairo), `.svg`, and `.png` (350 dpi) by the
scripts in `figures/_h/`, all of which run from `figures/_h/step_figures.sh` (submitted
from `figures/_m`). Shared visual language lives in `figures/_h/_fig_common.R`, sourced
by every figure script: Okabe–Ito palette, `theme_ms` (8–9 pt base), `save_fig`, and a
fixed RNG seed (`geom_jitter` draws at render time, so without it a figure could not be
regenerated). No in-panel titles — interpretation lives in the legends below — and bold
uppercase panel tags.

**Non-ASCII characters do not survive `cairo_pdf` here**: a literal `ρ`, `≥`, `β` or `×`
is emitted as `..`. Use plotmath (`expression(Spearman~rho)`, `parse = TRUE`) for
standalone titles, and spell symbols out in ASCII inside multi-line annotations.

## How to use this directory

The narrative arc the figures deliver:

> niche signals (CCC/NicheNet) → drive pericyte target programs → reframed as interpretable
> functional **states/programs** → arranged on a vascular-stabilizing ↔ injury/fibrotic
> **continuum** → summarized as a donor-level **niche-stability index** → with an
> **AT1R/AT2R balance** axis rationalizing AGTR1 blockade → robust to smoking/cohort
> confounders → and **conserved in mouse**.

**Two main figures:**

| File | Role | Built by |
|---|---|---|
| `figure_pericyte_layer.{pdf,svg,png}` | **Main Fig — pericyte layer (where → what/why).** Ties the localization "where AGTR1 is" to the state/continuum "what/why" on one shared UMAP, with the three-lens reversal as the linchpin. Self-contained (no vector-editor assembly needed). | `pericyte_layer_figure.R` |
| `figure_ccc_nichenet.{pdf,svg,png}` | **Main Fig — niche signaling.** Who signals to pericytes and what programs those signals drive. | `manuscript_mechanism_figure.R` |
| `figure_mechanism_main.{pdf,svg,png}` | **Main Fig — disease phenotype.** Donor-level niche-stability/injury readouts, state composition, and the continuum. Assemble alongside the image panels (state UMAP, DPT UMAP, PAGA) listed in `figure_panel_manifest.tsv`. | `manuscript_mechanism_figure.R` |
| `figure_disease_main.{pdf,svg,png}` | **Main Fig — disease association (continuous injury).** Graded injury-program engagement across disease groups, the programs that carry it, *AGTR1* disease effects across stromal cell types, and the independent COPD/IPF evaluation. Study-level robustness moved to S16 on 2026-07-29. Self-contained. | `disease_main_figure.R` |

**Supplements — S1–S16.** Filenames stay semantic so scripts and cross-references do not
churn; the manuscript number is carried by the `supp_number` column of
`figure_panel_manifest.tsv` and by the legend headings below.

| # | File | Built by |
|---|---|---|
| S1 | `figureS_pericyte_layer` | `pericyte_layer_figure.R` |
| S2 | `figureS_crossspecies_mouse` | `manuscript_mechanism_figure.R` |
| S3 | `figureS_state_annotation` | `state_annotation_figure.R` |
| S4 | `figureS_agtr1_dropout` | `agtr1_dropout_figure.R` |
| S5 | `figureS_acta2_control` | `manuscript_mechanism_figure.R` |
| S6 | `figureS_cogaps_validation` | `cogaps_validation_figure.R` |
| S7 | `figureS_continuum_stability` | `continuum_stability_figure.R` |
| S8 | `figureS_nichenet_specificity` | `nichenet_specificity_figure.R` |
| S9 | `figureS_receiver_robustness` | `receiver_robustness_figure.R` |
| S10 | `figureS_ras_landscape` | `basement_membrane_figure.R` |
| S11 | `figureS_state_composition` | `state_composition_figure.R` |
| S12 | `figureS_sensitivity` | `sensitivity_robustness_figure.R` |
| S13 | `figureS_program_category` | `manuscript_mechanism_figure.R` |
| S14 | `figureS_balance_by_state` | `manuscript_mechanism_figure.R` |
| S15 | `figureS_bm_copd` | `basement_membrane_figure.R` |
| S16 | `figureS_disease_robustness` | `disease_main_figure.R` |

The tail of the list is where displaced figures land, so that the numbering above stays put
when a new figure claims a slot:

- `figureS_program_category` held S3 until 2026-07-27, when the stability/annotation figure
  took that slot; **appended as S13 rather than shifting S4–S12 up by one**.
- `figureS_balance_by_state` held S5 until 2026-07-28, when the dropout figure entered at S4
  and pushed the ACTA2 control into S5; **appended as S14 rather than shifting S6–S13**.
- `figureS_bm_copd` held S6 until 2026-07-28, when the CoGAPS validation figure took that
  slot; **appended as S15 rather than shifting S7–S14**. It is a disease-side figure and
  belongs with the disease supplements; S15 keeps it numbered until that section is laid out.
- `figureS_disease_forest` was **new at S16** (2026-07-28). It was panel A of
  `figure_disease_main` until the continuous-injury rebuild replaced that panel with
  leave-one-study-out refits; it was appended rather than inserted for the same reason as
  the entries above.
- **Renamed to `figureS_disease_robustness` on 2026-07-29**, still S16. The leave-one-study-out
  panel came out of the main figure and joined it, so the file now holds two panels and
  "forest" no longer describes it. Same slot, same script, same manuscript number — only the
  filename and the panel count changed. `figureS_disease_forest.{pdf,svg,png}` was deleted;
  nothing outside this README and the manifest referenced it.

**Deliberately NOT numbered supplements:**
- `figureS_alluvial` — a **grant** figure, not a manuscript supplement (2026-07-27). It is
  still built by `manuscript_mechanism_figure.R` and still has a legend below; do not
  renumber it back into the S-series.
- `figureS_cogaps_transfer` — **retired** 2026-07-27. Its permutation-null panel became S8
  (A–B) and its CoGAPS-projection panel became S9C; keeping it would have duplicated both
  panels across two supplements. Note that it is *not* the predecessor of S6: it showed the
  CoGAPS→niche **projection**, whereas S6 shows rank selection and the pattern↔program
  correspondence. The two do not overlap.

**Supporting files (not figures):**
- `tableS_acta2_control.tsv` — source-data values underlying `figureS_acta2_control` (S4).
- `figure_panel_manifest.tsv` — two kinds of row. `figure = A|B` are the per-module PDFs
  (liana dotplots, NicheNet heatmap, state/DPT/PAGA UMAPs) assembled into the main figures
  in a vector editor; `figure = Supp` are the finished supplements keyed by `supp_number`.
  `exists` flags whether each has been generated.
- `figureB_states_continuum_niche.{pdf,png}` — **working draft only** (legacy
  `assemble_mechanism_figures.R`, `theme_bw`, in-panel titles). Superseded for submission by
  `figure_mechanism_main` + the manifest image panels; kept for quick QC, not for the manuscript.

**State-model key (applies to every panel that mentions programs).** Six stable pericyte
Leiden subclusters collapse onto three discrete dominant programs
(`pericyte_states/_m/annotations/state_program_map.tsv`): vascular-stabilizing (clusters 0, 2),
**basement-membrane** (clusters 1, 3, 5; 4,200 cells, 36%), activated/migratory (cluster 4).
Clusters 1/3/5 were previously labelled *fibroblast-like*; a pre-specified basement-membrane gate
(`basement_membrane/_h/01.state_gate.py`) showed that label was a least-negative default and the
positive identity is basement-membrane. Three further program *scores* (inflammatory,
synthetic/contractile, and the fibrillar-ECM **fibroblast-like** score — distinct from the
former state label) are continuous axes that are never a cell's dominant discrete program, so they
appear in score/continuum panels but not in discrete program-assignment panels.

---

## Figure legends

### `figure_pericyte_layer` — Pericyte localization reframed as a functional continuum

**Figure.** Integration of AGTR1 localization with pericyte functional state on a single shared
UMAP embedding (the embedding learned in the localization analysis and reused unchanged for state
discovery and trajectory, so all overlays share coordinates and cells; *n* = 11,680 pericytes).
**(A)** Pericyte UMAP colored by the six stable subclusters (P0–P5), hue families grouped by
program. **(B)** *AGTR1* expression (log) on the same embedding; *AGTR1* is detected diffusely
across the pericyte compartment rather than confined to one subcluster, supporting a
compartment-label rather than sub-state-marker interpretation. **(C)** Same embedding colored by
dominant state-program (vascular-stabilizing, basement-membrane, activated/migratory). **(D)**
*AGTR1* across programs under three measurement lenses — raw expression, binary detection, and
scVI-denoised — as donor-aware marginal means centered within lens (±SE). The raw and detection
lenses share a vascular-stabilizing–high pattern that **reverses** under denoising (denoised
*AGTR1* is flat-to-injury), identifying the raw enrichment as a transcript-capture/dropout effect
and establishing that *AGTR1* is not a vascular-stabilizing state marker. **(E)** Diffusion
pseudotime on the same embedding, ordering cells along a vascular-stabilizing ↔ injury/fibrotic
continuum (ordering reflects transcriptional similarity, not time). **(F)** Donor-level Spearman
correlation between each of the six program scores / *AGTR1* and pseudotime; orange = *P* < 0.05.
The five injury/mural scores rise monotonically along the continuum (donor ρ ≈ 0.19–0.52), and
*AGTR1* shows a weak positive trend (ρ ≈ 0.25), whereas the **basement-membrane score falls**
(donor ρ = −0.17, *P* = 0.016) — pericytes shed the structural BM program as they move toward the
injury pole. Together: where *AGTR1* sits (A–C), why it is a compartment label not a state
marker (D), and how the compartment is organized as an injury continuum (E–F).

### `figure_ccc_nichenet` — Niche signaling drives pericyte target programs

**Figure.** Cell–cell communication and ligand–target inference identify the niche signals
received by lung pericytes. **(A)** NicheNet ligand-activity ranking (top 15 ligands) for
predicting the pericyte target-gene program; bars show AUPR (corrected). TGF-β–family ligands
(*TGFB1/2/3*, *CCN1*, *CCN2*) are highlighted (orange). **(B)** Ligand→target regulatory-potential
heatmap for the prioritized ligands (rows) and their top predicted pericyte target genes
(columns); fill encodes NicheNet regulatory potential. **(C)** Dot plot of the fraction of cells
expressing each prioritized ligand across the niche sender cell types (size and color = fraction
expressing), identifying which neighboring populations supply each signal. **(D)** Donor-level,
covariate-adjusted validation of the predicted signaling axis. Each point is a donor; the **sender
ligand-composite score** (mean expression of the prioritized ligands across the niche sender cells)
and the **whole-pericyte injury-response score** (the curated NicheNet target program; receiver =
all pericytes, not a dropout-prone *AGTR1*⁺ split) were each residualized on dataset and disease
group, and the residuals are plotted with a single fitted line and 95% CI (partial regression /
Frisch–Waugh–Lovell). The sender→pericyte association persists after adjustment for study and
disease (**adjusted β = 0.50, *P* = 0.001, *n* = 105 donors**; the estimate is the ligand-score
coefficient from `lm(pericyte response ~ ligand score + dataset + disease group)`, equal to the
plotted residual slope). Disease group and dataset are intentionally **not** encoded by color or
shape, so the panel isolates the covariate-adjusted association. Abbreviations: CCC, cell–cell
communication; AUPR, area under the precision–recall curve. The cross-program CoGAPS projection
formerly shown as panel D is now Figure S9C.

### Figure S8 — `figureS_nichenet_specificity` — Specificity and donor-level validation of predicted niche regulation of pericyte injury programs

**Figure S8.** The prioritized ligands act specifically on the pericyte injury program, and the
prediction holds donor by donor. **(A)** For each of the 11 prioritized ligands, the observed
corrected AUPR against the curated pericyte target program (orange) compared with a null built
from 1,000 random gene sets of matched size (grey square = null mean, whiskers = mean ± 2 SD).
Every observation lies far outside its null, so the ranking reflects regulation of this specific
target program rather than generic connectivity of the prior network. **(B)** The same result as
a standardized distance: all 11 ligands lie 23.8–40.3 SD beyond their nulls and every one reaches
the empirical floor of the permutation test (*P* = 9.99 × 10⁻⁴ = 1/1001), so the permutation
resolves only that they are extreme, not how they rank against one another. **(C)** Donor-level
validation: niche expression of the prioritized-ligand composite against pericyte target-program
expression, one point per donor (*n* = 105), colored by disease group. Raw Spearman ρ = 0.54
(*P* = 4.4 × 10⁻⁹); adjusted for study as a random intercept, β = 0.65 (*P* = 1.2 × 10⁻⁶).
Disease colors are retained here — unlike the disease-neutral adjusted view in the main figure —
to show the relationship is not carried by one group. **(D)** The same for TGFB1 alone, which is
weaker but concordant (raw ρ = 0.24, *P* = 0.014; adjusted β = 0.41, *P* = 0.019), i.e. the
composite is not a restatement of a single ligand. **(E)** The same for TGFB2 alone. TGFB2 is the
**top**-ranked ligand in (A) by corrected AUPR (0.208, just above TGFB1's 0.203), so omitting it
would leave the leading ligand as the only one with no individual donor-level estimate. Its
donor-level support does not match its rank: raw ρ = 0.06 (*P* = 0.52), adjusted β = 0.54
(*P* = 0.087).

*Reading C–E.* In all three panels the **drawn line is the raw, unadjusted fit** while the
**annotated β is the dataset- and disease-adjusted LMM coefficient**; *n* = 105 donors throughout.
In C and D the two agree in direction and the distinction is cosmetic. In **E they do not**: the
line is flat while the adjusted β is the largest of the three. The adjustment is doing the work
there — donors with high TGFB2 are unevenly distributed across studies, so the marginal
correlation is confounded away and the within-study slope is positive but does not reach
significance. TGFB2's *x*-range is also compressed, with most donors below 0.15, so the estimate
rests on relatively few high-expressing donors. **Do not describe TGFB2 as donor-validated.** The
defensible statement is that the composite (C) is donor-validated, TGFB1 (D) individually so, and
TGFB2 (E) is not — despite leading the AUPR ranking. That dissociation between prior-network rank
and donor-level replication is the reason the manuscript's claims rest on the composite.

### Figure S9 — `figureS_receiver_robustness` — Predicted signaling is robust to receiver definition and differs between injury and basement-membrane programs

**Figure S9.** The prioritization does not depend on how the pericyte receiver is defined, and the
basement-membrane program has a distinct predicted sender set. **(A)** Ligand rank across nine
receiver definitions: all pericytes, the three stable-cluster programs, and five CoGAPS
dominant-pattern groups. Rows are the union of each receiver's top 8 ligands, ordered by their
rank under the primary all-pericyte definition; darker = better rank, grey = the ligand is not in
that receiver's candidate set (its receptor is not expressed there). All 36 receiver pairs agree
at Spearman ρ ≥ 0.77, and the two independently derived basement-membrane receivers — supervised
Leiden clustering versus unsupervised NMF — agree at ρ = 0.99 with 18 of 20 top ligands shared.
**(B)** Fibroblast- and myeloid-to-pericyte TGF-β and ECM edges (LIANA), faceted by receiver
scheme; the same senders reach pericytes under either definition. MMP14 is absent because LIANA's
resource carries no MMP14 ligand edge, not because the edge was tested and rejected. **(C)**
Pericyte-learned CoGAPS patterns (de-novo-selected nP = 8, one row per pattern labelled by its
assigned program) projected onto a cell-type × donor pseudobulk across 21 cell populations; fill =
within-pattern *z*. The patterns that rank pericytes first are the vascular-stabilizing and the two
basement-membrane patterns, the contractile pattern is strongest in vascular smooth muscle, and the
fibrillar/inflammatory patterns extend into bona-fide fibroblasts and myeloid populations — the
injury axes a subset of pericytes adopt are the ones fibroblasts constitutively run. **(D)**
Re-running the prioritization against the **basement-membrane** target set instead of the injury
program gives a substantially different ranking (Spearman ρ = 0.223 against the injury ranking,
only 8 of the top 20 ligands shared). The target set is the **10 of 13** BM-panel genes expressed
in pericytes above the 0.10 receiver threshold (*LAMA3*, *LAMA5* and *AGRN* fall below it and were
never tested). Of 321 candidates only **MMP14** survives correction across 10,000 matched
permutations (rank 1, empirical *P* = 9.999 × 10⁻⁵, FDR = 0.032); TIMP2 ranks second but is not FDR-significant,
and TGFB2/TGFB1 carry regulatory potential without reaching significance. BM deposition is
therefore predicted to be governed by pericyte-proximal matrix turnover rather than by the
TGF-β-dominated axis that drives the injury program.

### `figure_mechanism_main` — Donor-level disease phenotype and the stabilizing↔injury continuum

**Figure.** Pericyte niche state across health and disease at the donor level. **(A)**
Niche-stability index and **(B)** injury-stromal score per donor, grouped by disease
(Healthy, COPD, Fibrotic/ILD); box = median and IQR, whiskers = 1.5×IQR, white diamond = mean,
points = individual donors; brackets show Wilcoxon rank-sum comparisons versus Healthy. **(C)**
AT1R–AT2R pathway balance per donor (shown as a **corollary** of injury intensity: the disease
effect is redundant with the injury-stromal score and collapses after adjustment — see
`MECHANISM_ANALYSES.md`). **(D)** Mean donor-level pericyte state-program composition by disease
(stacked fractions across the three discrete programs — vascular-stabilizing, basement-membrane,
activated/migratory — plus the continuous-program assignments).
**(E)** Continuum trends: donor-level Spearman correlation between diffusion-pseudotime ordering
and each of the six program scores / *AGTR1*; points colored by significance (orange, *P* < 0.05).
The ordering reflects transcriptional similarity along a vascular-stabilizing ↔ injury/fibrotic
**continuum**, not a temporal axis. Vascular-stabilizing, inflammatory, synthetic/contractile,
and activated/migratory scores rise monotonically along the continuum (donor ρ ≈ 0.46–0.52,
*P* < 1×10⁻¹¹); the **basement-membrane score falls** (donor ρ = −0.17, *P* = 0.016); *AGTR1*
shows a weak positive trend (ρ ≈ 0.25). Assemble with the state UMAP,
DPT-pseudotime UMAP, and PAGA panels in `figure_panel_manifest.tsv`. *n* = 32 donors (donor-level
mixed-model marginal means; df = 31).

### `figure_disease_main` — Disease is associated with graded pericyte injury-program engagement, while *AGTR1* dysregulation sits in fibroblasts

*Rebuilt 2026-07-28 for the continuous-injury revision, trimmed to four panels 2026-07-29.
Panel provenance relative to the pre-2026-07-28 three-panel version: old C → **A** (now three
groups), old B → **B** (now both contrasts); **C** and **D** are new. The two study-level
robustness panels — the within-study forest (old A) and the leave-one-study-out refits that
briefly replaced it at C — now sit together as the two panels of **S16**: they answer two
different objections to the same panel-A estimate, so they belong beside each other, and the
main figure carries the biological chain A → B → C → D uninterrupted.*

**Figure.** Donor-level association between lung disease and pericyte injury-program
engagement, and the receptor-level question that follows from it. The figure begins from an
established null — disease is **not** associated with replacement of discrete pericyte states
(`figureS_state_composition`, `figure_mechanism_main` D) — and shows instead that disease
tracks a **graded, continuous** shift in program engagement. The endpoint in **A–B** (and in
both panels of **S16**) is
the donor-level **injury-program score**: the mean of the *z*-standardized donor means of the
inflammatory, activated/migratory, and fibrillar fibroblast-like per-cell program scores
(continuous; the basement-membrane program and the *AGTR1*⁺ fraction are deliberately excluded —
basement membrane is vascular-support biology, and binary *AGTR1* detection is dropout-prone and
would fold the focal receptor into the outcome). Models are donor-level linear mixed models
`~ disease_group + sex + (1 | dataset)` (`disease_association/_h/03.disease_forest.R`).

**(A)** Injury-program score by disease group. Each point is a donor; the black point and
whisker are the **sex + study-adjusted marginal mean ± 95 % CI**; brackets give the contrast
against Healthy. Engagement is **graded**: Fibrotic/ILD **+0.79 SD (95 % CI 0.38–1.20,
*P* = 3.5×10⁻⁴)** and Other disease **+0.62 SD (95 % CI −0.07–1.31, *P* = 0.086)** (*n* = 88
donors; 42 Healthy, 24 Fibrotic/ILD, 22 Other). **(B)** Which continuous programs carry the
score — the same three-group model fit to each program separately (points = contrasts vs
Healthy, whiskers = 95 % CI, asterisk = *P* < 0.05). The composite is driven by the **fibrillar
fibroblast-like** program (Fibrotic/ILD **+1.13 SD, *P* = 4.4×10⁻⁴**; Other **+1.10 SD,
*P* = 0.016**) and the **activated/migratory** program (Fibrotic/ILD **+0.87 SD, *P* = 0.006**),
with a smaller inflammatory contribution (**+0.55 SD, *P* = 0.040**) and a **flat
vascular-stabilizing** program (+0.15 SD, *P* = 0.46) — a targeted fibrogenic shift, not a
global rescaling. The estimate in **A** is robust to study composition from both directions —
no single cohort creates it, and it is reproduced *within* studies that sampled both arms — see
**S16**. **(C)** *AGTR1* disease effects
across stromal cell types (`disease_association/_h/05.agtr1_celltype_disease.R`): the
**omnibus** (2 df) disease effect on donor-mean *AGTR1*, expressed as partial η² — the share of
donor-level *AGTR1* variance attributable to disease group, computed within cell type so that
cell types with ~5× different baseline *AGTR1* are comparable. Rows are **blocked by lineage**
(fibroblast above mural) and sorted by η² within each block, so the class-level comparison the
panel makes is a property of the layout rather than a coincidence of a single sort order; the
strip labels replace a colour legend and the colours are redundant encoding. All three
**fibroblast** populations rank above both **mural** populations (peribronchial 0.195,
adventitial 0.091, alveolar 0.071 vs pericytes 0.052, vascular smooth muscle 0.043; lineage
means 0.119 vs 0.048).
**(D)** Independent evaluation in **GSE136831** (Adams/Kaminski), the only dataset here with a
real COPD arm (`disease_association/agtr1_copd_ipf/`). Donor × compartment pseudobulk *AGTR1*,
disease minus Control, adjusted for depth, sex, age and **ever-smoker** status (which the HLCA
metadata does not permit). Filled points are the pre-specified fibroblast-lineage family, open
points exploratory. *AGTR1* is **lower in myofibroblasts in both arms** (COPD **−0.18,
*P* = 0.037**; IPF **−0.17, *P* = 0.028**) — the same compartment, the same direction, in two
independent disease arms — directionally consistent with the negative fibroblast contrasts in
**C**; bulk fibroblasts are null (COPD −0.09, *P* = 0.17; IPF +0.03, *P* = 0.61). The primary
family is reported at **nominal alpha**: these two compartments were fixed by the independent
HLCA analysis in **C** before GSE136831 was touched, so this is a directional replication rather
than a screen, and the two myofibroblast tests are the same compartment against a shared Control
arm (correlated, not independent looks). A BH value is carried in the source table
(`p_BH_reference` = 0.074) for readers who want it, but it is not the reporting gate; the guard
against fishing is the pre-specification plus the `family` column, which labels every other
compartment and gene exploratory.

**Statistical notes / limitations** (see `MECHANISM_ANALYSES.md`):
- **Scale.** A and B read off one three-group model whose programs are *z*-standardized once
  over the three-group donor set, so they are mutually comparable, as is **S16A**. They are
  **not** on the same scale as the two-group primary contrast in **S16B** and must not be quoted
  interchangeably.
- **COPD is excluded from A–B.** All 12 HLCA COPD donors come from a single study
  (Kaminski_2020), so a COPD estimate there is inseparable from that study. COPD is evaluated
  instead in **D**, where it is a real independent arm. "Other disease" in A–B is
  COVID/carcinoma/etc., not COPD.
- **Donor gates differ between the HLCA panels and D, because the datasets differ.** A–B:
  **≥ 10 pericytes per donor** (`03 --min-cells 10`), and ≥ 2 donors per group for a study to
  enter S16B; the endpoint is insensitive to this — the effect holds across min-cells 5–30
  (all *P* < 0.005), and 5 and 10 give identical donor sets because no donor has 5–9 pericytes.
  C: **≥ 10 cells per donor per cell type** (from `01.disease_association.R`), ≥ 3 donors per
  cell type, plus **≥ 3 donors in both Healthy and Fibrotic/ILD** for a cell type to be tested
  at all — 3 being the smallest group that admits a within-group SD, i.e. the point at which the
  contrast stops being estimable rather than a power preference. This gate is **not binding**:
  the three excluded stromal populations have 1, 1 and 0 fibrotic donors (mesothelium,
  myofibroblasts, subpleural fibroblasts), so the tested set is the same five cell types at any
  threshold from 2 to 6. D: **≥ 5 cells per donor per compartment** and ≥ 5 donors per arm — a lower floor
  because GSE136831 is a smaller, disease-enriched cohort, which is precisely why the pericyte
  compartment still fails it.
- **Panel D gives only qualified support**, by design and not by outcome. GSE136831 has exactly
  **one** Control donor with ≥ 5 pericytes (against 6 COPD and 15 IPF), so a pericyte-specific
  contrast is **not estimable in either direction**; that row is drawn explicitly rather than
  omitted. This dataset can corroborate the **fibroblast** half of the claim only, and it does
  **not** establish consistent pericyte-specific *AGTR1* upregulation.
- **Nothing in C reaches significance** (all BH ≥ 0.71), and the fibrotic groups are small
  (6 donors per cell type). C supports an **ordering** — disease-associated *AGTR1* variation is
  larger in fibroblasts than in pericytes — not a positive fibroblast finding. Directional
  contrast estimates and CIs are in
  `disease_association/_m/mean_expr/agtr1_celltype_disease_effects.tsv`.
- **C supersedes the earlier Kruskal–Wallis *AGTR1*-by-cell-type analysis**
  (`disease_agtr1_kruskal_by_celltype.tsv`, from `01.disease_association.R`), which reported
  peribronchial fibroblasts as significant (KW *P* = 0.0014, FDR = 0.0068). That test used the
  raw CELLxGENE `disease` vocabulary (4–6 unbalanced groups per cell type) and modelled **no
  study term**. It should not be quoted: in the stroma subset the fibrotic donors of every
  tested cell type come from only **2–3 datasets, and only one dataset per cell type contains
  both a Healthy and a Fibrotic/ILD donor**, so the old contrast is largely *between* studies.
  In peribronchial fibroblasts specifically, all 6 fibrotic donors have donor-mean *AGTR1*
  exactly 0 and come from Lafyatis_2019 and Schiller_2020, and Schiller_2020's single Healthy
  peribronchial donor is also 0 — so the comparison is effectively one study's cases against
  Healthy donors pooled from twelve others. Adding `(1 | dataset)` absorbs this and the
  significance does not survive; the **direction** does (fibrotic lower, −0.86 SD). This is the
  same failure mode already documented for the pericyte niche-index disease effect.
- **Near the detection floor.** *AGTR1* in fibroblast populations is low (Healthy donor means
  ≈ 0.10–0.18 log-normalized in peribronchial/alveolar fibroblasts vs ≈ 0.85 in pericytes), and
  several fibrotic donor means are exactly 0. C is therefore partly a **detection** comparison,
  and inherits the dropout caveat that applies to *AGTR1* throughout this project.
- **Age was not adjusted** in the primary model because it is missing for 18/24 fibrotic donors
  (forcing it collapses the case group to *n* = 6; direction preserved, +0.22 SD, *P* = 0.55).
  **Smoking-adjusted disease contrasts are not estimable** (smoking is unrecorded for all
  fibrotic donors), and within healthy donors smoking does not move the endpoint
  (ever − never = −0.12 SD, *P* = 0.48), so it is not a hidden driver.
- The association is reported as a cross-sectional, cohort-anchored transcriptional correlate of
  disease, not evidence of pericyte loss or causal progression. **This figure supersedes the
  donor-level disease panels (A, B) of `figure_mechanism_main`**, which used a pooled composite
  ANCOVA whose only significant contrast was the off-target Healthy-versus-"Other" comparison.

### `figureS_alluvial` (grant figure, unnumbered) — Stable subclusters → dominant program → effector class

**Figure S.** Alluvial flow linking the six stable pericyte Leiden subclusters (left axis) to
their three dominant programs (middle axis: vascular-stabilizing, basement-membrane,
activated/migratory) to six functional effector molecule classes (right axis: signaling ligands;
chemokines/cytokines; adhesion molecules; ECM structural; matrix-remodeling enzymes; fibrotic
mediators). Left→middle ribbon width is proportional to subcluster cell count; middle→right width
is each program's effector composition (donor-averaged mean fraction of cells expressing the
class, normalized within program so flow is conserved). Ribbons are colored by **stable
subcluster**, with hue families grouped by program (blues = vascular-stabilizing [P0, P2];
oranges = basement-membrane [P1, P3, P5]; pink = activated/migratory [P4]) so both the
subcluster→program merge and the program grouping are legible. Cluster→program assignments are
the canonical `state_program_map.tsv` used throughout the analysis.

### Figure S13 — `figureS_program_category` — Program × protein-category enrichment

**Figure S.** Dot heatmap relating the six pericyte programs (rows) to eight curated marker
protein categories (columns: signaling ligands; chemokines; cytokines; adhesion molecules; ECM
structural; matrix-remodeling; fibrotic mediators; mural identity). Dot size = prevalence (mean
fraction of the program's cells expressing the category); fill = relative enrichment (*z* of
mean detection across programs, computed within category to expose program-specific signal
otherwise swamped by category-level baseline detectability). Cells are assigned to their dominant
program by *z*-scored argmax of the six program scores — now including basement-membrane
(`cell_communication/_h/03b.program_category_enrichment.py`), retaining all six programs.

### Figure S3 — `figureS_state_annotation` — Stability and annotation of human lung pericyte states

**Figure S.** Audit trail for the six stable pericyte states used throughout the study
(*n* = 11,680 pericytes). States are data-driven Leiden clusters on the study-integrated
embedding (`X_pca_harmony`) chosen by a bootstrap stability sweep; curated marker panels and
Wilcoxon differential expression **annotate** the clusters and played no part in defining them.
**(A)** Stability sweep over Leiden neighbourhood size (15, 30) and resolution (0.3–0.9),
30 bootstraps per setting at 80% resampling. Points are the median best-match Jaccard between
the full-data clustering and each bootstrap replicate; labels give the number of clusters *k*;
the dashed line is the 0.6 pass gate (`--stability-threshold`). Every setting clears the gate by
a wide margin, so the selected solution (black ring: 30 neighbours, resolution 0.5, *k* = 6;
median Jaccard 0.966, mean ARI 0.946, mean NMI 0.932) was chosen among stable alternatives
rather than rescued by the threshold. **(B)** The selected solution broken out per cluster, with
cluster size on the axis: point = median over 30 bootstraps, bar = 2.5–97.5% range, cross = worst
replicate, dotted line = median across clusters (0.965), numbers = per-cluster medians. **This is
the panel that qualifies (A).** The four large clusters are highly reproducible (P0 0.99, P1 1.00,
P2 0.91, P3 0.98), but the two smallest are not: **P4** (220 cells, activated/migratory) has a
median of only **0.52** with a worst replicate of 0.11, and **P5** (113 cells,
basement-membrane) is bimodal — median 0.95 but a 2.5th percentile of 0.03, i.e. in a minority of
resamples it dissolves entirely. The pooled median in (A) is taken over all 180 cluster ×
bootstrap values, so the 120 values contributed by P0–P3 hold it at 0.97 and it does **not**
communicate this. **(C)** What the clusters express.
*Top row:* data-driven Wilcoxon markers per cluster (top four by score, `pval_adj` < 0.05 and
logFC ≥ 0.25); ribosomal-protein and mitochondrial genes are excluded **from display only**, and
the unfiltered ranking is in `pericyte_states/_m/annotations/state_markers.tsv.gz`. *Bottom row:*
the curated `STATE_PANELS` used to annotate the clusters, four genes per program, ranked by how
much each discriminates among P0–P5 (SD of per-cluster mean expression) — most of the
inflammatory panel is near zero in every pericyte cluster and would otherwise fill the block with
blanks. Panels overlap (e.g. *COL4A1* is on both the fibroblast-like and basement-membrane
panels), so a gene is kept only for the first program that claims it. In both rows, dot size =
fraction of the cluster's cells expressing the gene and fill = *z* of mean log-expression across
the six clusters, so a row reads for specificity. **(D)** Relative program enrichment per cluster
(mean of cell-level *z*-scored program scores); black outline marks the dominant program assigned
to that cluster, giving P0/P2 → vascular-stabilizing, P1/P3/P5 → basement-membrane,
P4 → activated/migratory.

*Caveats.* (0) **P4 and P5 are not stable clusters at the level P0–P3 are** (see B). Neither
supports a claim that rests on it being a discrete, reproducible cell state. This is consistent
with how they are already used: the disease analyses report the *continuous* injury-program
score, and **S11** shows discrete composition does not differ by disease group — so nothing in
the study's conclusions rests on P4 or P5 as entities. It does mean the activated/migratory
program label in (D) is carried by the least reproducible cluster in the solution, and any text
describing P4 as a distinct state should be softened to a program-score statement.
(i) P0 and P2 — the two vascular-stabilizing clusters — are not separated by
distinctive lineage markers in the top row of (C); their top Wilcoxon genes are largely
translational and housekeeping transcripts, i.e. the two differ substantially in transcriptional
load rather than in identity. That is a property of the data, not of the display filter, and it
is why the functional claims in this study rest on the continuous program scores rather than on
P0-versus-P2 as separate biological entities. The curated row of (C) is what carries the
cluster→program link for those two. (ii) There is deliberately no program-score UMAP panel here:
**S1** (`figureS_pericyte_layer`) panel A already shows those six overlays on this same
embedding, and a panel duplicated across two supplements is worse than a cross-reference.

*Provenance.* Panels A and D come straight from `00.state_discovery.py`
(`stability/cluster_stability_grid.tsv`, `annotations/state_program_map.tsv`). Panels B and C
need two quantities that script aggregates away — per-cluster (rather than pooled) bootstrap
Jaccard, and per-cluster mean/detection expression for the marker genes. Both are produced
**without re-running the clustering**, because `pericyte_state` labels are keyed on by four
downstream modules: (B) by `00.state_discovery.py --stability-only` (driver `step_0c.sh`), which
re-runs the bootstrap at the stored solution and refuses to write unless the recomputed labels
reproduce `pericyte_state` for every cell; (C) by `00b.annotation_support.py` (driver
`step_0b.sh`), which only reads the published object. Note that (B)'s Jaccards are a **fresh**
bootstrap, not a replay of the sweep — the sweep consumes one RNG stream across all eight
settings, so a single-setting run necessarily draws a different resample sequence. (A)'s 0.966
and (B)'s median are therefore two draws of the same quantity, agreeing within the 0.05 tolerance
the script enforces; do not quote them as the same number.

### Figure S1 — `figureS_pericyte_layer` — Pericyte-layer supporting detail

**Figure S.** Companion detail for `figure_pericyte_layer`, all on the same shared UMAP
(*n* = 11,680 pericytes). **(A)** Per-program module-score overlays (vascular-stabilizing,
synthetic/contractile, activated/migratory, inflammatory, fibroblast-like, basement-membrane),
showing the spatial organization of each of the six programs across the embedding. **(B)** *ACTA2* expression (log), the canonical
contractile-mural benchmark, which mirrors the raw *AGTR1* vascular-stabilizing–high pattern.
**(C)** Binary *AGTR1* detection, visualizing the dropout structure that drives the raw
enrichment reversed in `figure_pericyte_layer` panel D. **(D)** Donor-mean *AGTR1* versus *ACTA2*
expression by stable subcluster (donor × subcluster means, subclusters with ≥5 cells; box =
median/IQR, whiskers = 1.5×IQR), comparing the two markers' subcluster profiles.

### Figure S4 — `figureS_agtr1_dropout` — Apparent AGTR1-undetected pericytes reflect transcript dropout rather than a distinct spatial population

**Figure S.** *AGTR1* is detected in only 4,371 of 11,680 human lung pericytes, which invites the
reading that the remaining cells are an *AGTR1*-negative subpopulation residing elsewhere in the
lung. Five panels test that reading and reject it.
**(A)** Observed versus expected number of *AGTR1*-undetected pericytes under a matched-gene
dropout model. The null is the observed zero rate of the **200 genes closest to *AGTR1* in pooled
abundance** (2.25–3.07 counts per 10k versus *AGTR1*'s 2.63) in the **same** cells, with
mitochondrial and ribosomal genes excluded as ambient-dominated. Observed 7,314 of 11,680 cells
(62.6%) versus a matched expectation of 7,459 (63.9%) [95% matched-gene interval 6,562–9,279];
*z* = −0.19, empirical *P* = 0.85. *AGTR1* is, if anything, marginally **less** dropout-prone than
its abundance peers.
**(B)** Airspace-affinity distributions (mean cosine similarity to AT1/AT2/aerocyte/general-capillary
centroids in `X_pca_harmony`) for the two groups; violin = full distribution, box = median/IQR.
All 11,680 pericytes.
**(C)** Effect of the *AGTR1* readout on airspace affinity across five specifications, all fitted on
one common complete-case set (*n* = 7,157 cells, 56 donors) with the outcome z-scored, so rows are
directly comparable in SD units. Per-cell rows are LMM with a donor random intercept and a dataset
variance component; the donor row is OLS on donor means. Binary detectability: −0.136 unadjusted,
−0.117 with depth covariates, **−0.114 fully adjusted** (age, sex, `n_counts`, `n_genes`,
`pct_mito`; *P* = 1.1 × 10⁻⁵). Continuous scVI-denoised *AGTR1*: −0.087 SD per SD
(*P* = 6.1 × 10⁻⁹). **Donor fraction *AGTR1*-detectable: −0.064 [−0.555, 0.427], *P* = 0.79** — the
largest possible contrast on the ladder (a 0→1 predictor) and the row that carries the panel.
**(D)** Among the 7,309 cells with **no observed *AGTR1* transcript**, scVI-denoised *AGTR1* (log
scale, pericyte-trained model) against airspace affinity. Two readings at once: the undetected
cells are **not** *AGTR1*-low — their denoised values span two orders of magnitude and their median
(0.770) **exceeds** that of the detected cells (0.744, dashed rule) — and within them the denoised
level carries no airspace signal (Spearman ρ = −0.02, *P* = 0.052).
**(E)** *AGTR1* detection rate by sequencing-depth decile (orange, 95% Wilson interval) against the
matched genes' detection rate in the same cells (grey, median and 2.5–97.5 percentile). Detection
climbs monotonically from 12.4% in the shallowest decile (median 712 UMIs) to 61.6% in the deepest
(median 6,118 UMIs), tracking the matched-gene median (13.3% → 61.4%) throughout.

**Reading A and E together.** They are complementary, not redundant: **A** asks whether the *number*
of zeros is unusual, **E** asks whether the zeros *land where dropout predicts*. A gene with a
normal zero count but depth-independent zeros would pass A and fail E; *AGTR1* passes both.

**What this figure does not claim.** It is **not** a null result. The per-cell effect in C is
significant, and it is significant on the scVI-denoised *continuous* lens as well — so a real,
graded association between *AGTR1* level and distance from the airspace does exist, as expected for
a perivascular receptor. The claim is narrower and is about **detectability as a population
boundary**: the zero count is what dropout predicts (A), the two groups' distributions coincide (B),
the contrast does not survive donor aggregation (C), and the zeros are internally graded with no
spatial structure of their own (D).

**Caveats.**
1. The two bands in E mean different things and must not be read as comparable error bars: the
   orange interval is *sampling* error for one gene (binomial/Wilson), the grey band is
   *between-gene* spread across 200 genes. That is why grey is far wider.
2. The `counts` layer of these AnnData objects is **soupX ambient-corrected, not raw** (float,
   minimum 0.00455). The dropout model is fitted on `raw/X` of
   `pericyte.hlca_full.dataset.h5ad`, the only integer matrix in the project. A sampling model
   fitted to the soupX layer would be meaningless.
3. Two definitions of detection coexist. The project's `AGTR1_detect` is soupX-derived (4,371
   detected); raw counts give 4,366. They agree for **99.96%** of cells. Panels A/E use raw
   detection so they are internally consistent with the matched genes; B/C/D use the project
   definition so C stays comparable with the published estimate.
4. The fully adjusted row of C reproduces the previously published per-cell estimate exactly
   (−0.003587 airspace units ÷ SD 0.0316 = −0.1136 SD), which is the intended cross-check on the
   re-implementation.
5. A Poisson reference (6,469 expected zeros) is written to `dropout_model_summary.tsv` but is
   **not** plotted and is not a fair test: real counts are overdispersed and Poisson under-predicts
   zeros for essentially every gene. The empirical matched-gene null is the primary comparison.
6. The complete-case restriction to 7,157 of 11,680 cells is driven by missing donor age/sex in
   HLCA, and by dropping donors that contribute only one *AGTR1* class (no within-donor contrast).

**Provenance.** Upstream: `localization/airspace_analysis/_h/03.agtr1_dropout.py`
(driver `step_3.sh`, run from `localization/airspace_analysis/_m`), writing
`_m/agtr1_dropout/{dropout_model_summary,matched_gene_null,detection_by_depth,airspace_effect_ladder}.tsv`
and `agtr1_cells.tsv.gz`. Figure: `figures/_h/agtr1_dropout_figure.R`.

### Figure S5 — `figureS_acta2_control` — AGTR1 is not reducible to ACTA2⁺ contractile identity

**Figure S.** Control benchmarking *AGTR1* against canonical *ACTA2*⁺ contractile mural identity.
**(A)** Centered donor-aware program marginal means (±SE) for *ACTA2* (raw), *AGTR1* (raw), and
*AGTR1* (scVI-denoised) across the three discrete programs. The two **raw** markers share a
vascular-stabilizing–high pattern; the **denoised** *AGTR1* does not, indicating the raw pattern
is a shared transcript-capture/dropout effect. **(B)** Forest plot of donor × program pseudobulk
Pearson correlation (point) with Fisher-*z* 95% CI between *AGTR1* (raw vs denoised) and either
*ACTA2* or the leave-*ACTA2*-out contractile program; filled = *P* < 0.05, open = n.s. The raw
coupling (*r* = 0.28 vs *ACTA2*; *r* = 0.39 vs contractile program) collapses after denoising
(*r* = 0.04, *P* = 0.65; *r* = 0.16, *P* = 0.054), so *AGTR1* is **not** reducible to contractile
identity. *n* = 154 donor × program pseudobulks. Source data: `tableS_acta2_control.tsv`.
Framed as a benchmark/control, not a mechanistic pillar.

### Figure S14 — `figureS_balance_by_state` — AT1R–AT2R balance by pericyte program

**Figure S.** Donor-level AT1R–AT2R pathway balance across pericyte programs (box = median/IQR,
points = donors). Differences among programs are **not significant** (smallest pairwise
*P* = 0.067); this panel documents that the AT1R/AT2R skew tracks continuous injury intensity and
disease rather than discrete program identity, justifying its demotion from the main figure.

### Figure S6 — `figureS_cogaps_validation` — Unsupervised CoGAPS validation of pericyte programs

**Figure S6.** The curated six-program model is marker-panel driven and can therefore only recover
what its panels encode. CoGAPS is given no panels: it factorizes the same 2,038-HVG × 11,680-pericyte
matrix into *nP* patterns de novo (Bayesian NMF, distributed mode, 5,000 iterations, 8 subsets).
**(A)** Rank selection. Cross-seed reproducibility across *nP* = 4–10: for each rank, the canonical
(seed-13) fit is the reference and each of three non-canonical seeds (1, 42, 2024) is greedily
1:1-matched to it by |Pearson *r*| over shared gene loadings. The line is the mean over patterns, the
band the min–max; the **floor** is what separates the ranks, since the mean stays ≥ 0.92 everywhere.
Against a 0.80 threshold, *nP* = 5 (floor 0.700), 6 (0.577) and 10 (0.405) fail on a single collapsing
pattern, while 7 (0.951), 8 (**0.978**) and 9 (0.952) pass. **nP = 8 is the main rank** — it maximizes
the weakest pattern's reproducibility — and **nP = 9 is the sensitivity rank**, the largest that still
clears the threshold (the recommendation rule in `02.select_rank.R` reports 9 for exactly that reason;
both are carried forward, and E shows the choice does not matter). Reconstruction error is the other
standard rank signal and does not adjudicate: unexplained variance falls monotonically from 0.727 at
*nP* = 4 to 0.694 at *nP* = 10 with no elbow, so the rank is chosen on reproducibility alone.
**(B)** Do the de-novo patterns carry the curated programs? Cell-level Spearman of each pattern's
weight against each curated program score at *nP* = 8 (*n* = 11,680). Every curated axis has a leading
pattern: vascular-stabilizing → P7 (ρ = 0.58), synthetic/contractile → P5 (0.62), inflammatory → P3
(0.45) and P1 (0.40), fibroblast-like → P4 (0.42) and P6 (0.38), basement-membrane → P8 (0.39).
**(C)** Is that agreement visible in the *genes*, not only the scores? Overlap of each pattern's top-50
PatternMarker genes with the six curated panels, hypergeometric against the 2,038-gene HVG background,
BH-adjusted over the 48 tests. Three survive, and they are the three the score panel predicts:
P7 × vascular-stabilizing (5 of 11 panel genes, BH = 1.4 × 10⁻⁴), P5 × synthetic/contractile (4 of 8,
BH = 5.0 × 10⁻⁴), P3 × inflammatory (4 of 16, BH = 7.5 × 10⁻³). Crosses mark measured zeros.
**(D)** Is the *nP* = 8 solution itself stable? All 24 pattern × seed matches reach |*r*| ≥ 0.958, and
the weakest pattern (P6) averages 0.978 — the number quoted as the *nP* = 8 floor in A.
**(E)** Does any of this depend on picking 8 over 9? All **8 of 8** *nP* = 8 patterns match an *nP* = 9
pattern (median |*r*| = 0.998, min 0.867), the matched pairs keep the same six-program correlation
*profile* (median *r* = 0.997, min 0.853), and the leading curated program agrees for **8 of 8**. The
one extra pattern *nP* = 9 buys is unmatched and correlates weakly with every curated panel (max
ρ = 0.16), i.e. it adds no new program.

**Reading B and C together.** They are independent lines of evidence and a pattern can pass one and
fail the other, which is why both are shown. P6 does exactly that: fibroblast-like by cell-level score
(ρ = 0.38) with **zero** marker-gene overlap against any panel. C is also the harsher test — the panels
are 7–16 genes inside a 2,038-gene background, so an overlap of 3 is already unusual and the
non-significant cells are mostly small-*n*, not contradictory.

**What this figure does not claim.** It does not claim the CoGAPS patterns *are* the curated programs
one-for-one. The mapping is many-to-one (two patterns lead on inflammatory, two on fibroblast-like,
two on basement-membrane), and P2 and P8 correlate **negatively** with five of the six panels — they
behave as low-activity/baseline factors, and their basement-membrane argmax label rests on a single
modest positive ρ (0.25 and 0.39). The claim is narrower and is what A–E support: the curated program
*axes* fall out of an unsupervised factorization of the same cells, at a rank the data chose, and
neither the axes nor their gene evidence depend on that rank.

**Caveats.** (1) Argmax over ρ is a weak label wherever the winning ρ is small; the underlying values
are in `cogaps_np8_vs_np9_concordance.tsv` (`top_rho_*`) so a reader can see which labels are strong.
(2) Distributed CoGAPS matches patterns across data subsets to form the consensus and is free to
return a count other than the one requested — the *nP* = 6 fit came back with **7** patterns. It is
plotted in A at the rank requested, and it is not a rank this figure relies on.
(3) A pattern that collapses in one seed has no correlate to match against and yields no *r*; those
are dropped from the mean rather than scored as zero (1 of 15 at *nP* = 5, 3 of 21 at *nP* = 6, 1 of 27
at *nP* = 9, **none at *nP* = 8**). Dropping them is the convention in `02.select_rank.R` and is
conservative *against* the ranks that have them, since a collapse is itself instability.
(4) B and C are computed at *nP* = 8 only; the *nP* = 9 counterparts are in `validation_np9/` and E is
the bridge between them. (5) This figure is about the pattern ↔ program correspondence; the CoGAPS →
niche **projection** is a different analysis and lives in S9C.

Source data: `pericyte_cogaps/_m/cogaps_nP_selection.tsv`, `cogaps_seed_stability.tsv`,
`validation_np8/{pattern_score_spearman,pattern_panel_overlap}.tsv.gz`,
`cogaps_np8_vs_np9_concordance.tsv` (written by `_h/06.rank_concordance.R`, `step_6.sh`).

### Figure S7 — `figureS_continuum_stability` — Stability of the pericyte transcriptional continuum across parameter and root choices

**Figure S7.** The vascular-stabilizing ↔ injury continuum is not an artifact of the analyst's
choices. Diffusion pseudotime was recomputed across a grid of neighborhood sizes, diffusion-component
counts, cell subsamples and root cells; every correlation shown is **donor-level** (Spearman of
per-donor mean pseudotime against per-donor mean score), aggregated exactly as the headline analysis
does. **(A)** Donor-level ρ between pseudotime and each of the six program scores plus *AGTR1*
across the 18 parameter settings (columns: n_neighbors / n_DCs / subsample %). **(B)** The spread of
those 18 estimates per feature; the label is the fraction of settings sharing the sign of the mean,
which is **100% for all seven features**. The estimates are tight (s.d. ≤ 0.031, full range ≤ 0.08
for every feature), and the canonical setting reproduces the headline analysis exactly
(vascular-stabilizing ρ = 0.520, *AGTR1* ρ = 0.249, basement-membrane ρ = −0.172, identical to the
values in `pseudotime_trend_correlations.tsv`). **(C)** The same correlations after re-rooting pseudotime at each of
the six stable clusters in turn and at both extremes of the leading latent axis. **(D)** |ρ| under
the canonical vascular-stabilizing root against |ρ| under each alternative root (open symbols = the
sign flipped); points fall on the identity line at Spearman ρ(|ρ|) = 0.92. Re-rooting at the
opposite pole reverses the *direction* of pseudotime, as it must, but leaves the strength of every
feature's association intact — changing the root changes the orientation, not the underlying axis.

### Figure S11 — `figureS_state_composition` — Discrete pericyte-state composition does not differ across disease groups

**Figure S11.** A deliberate null that bounds what the discrete state model can be asked to do.
Donor-level ANCOVA (fraction ~ disease group + age + sex) on donors with ≥ 20 pericytes; boxes and
points are donors, diamonds are age/sex-adjusted marginal means ±95% CI. **(A)** Fractions of the
six stable pericyte clusters and **(B)** of the three dominant programs, by disease group.
**(C)** The grouped injury-associated state fraction (Healthy 0.019, Fibrotic/ILD 0.034, Other
0.034). **(D)** Forest of every disease contrast against Healthy across all six clusters, the three
programs and the grouped fraction. **No contrast is significant after Benjamini–Hochberg
correction** (all BH *P* ≥ 0.47 for the grouped fraction; ≥ 0.65 for programs; ≥ 0.87 for
clusters), and every interval crosses zero. Intervals are nominal 95% while *P* values are
BH-adjusted, so the two need not agree exactly at the 0.05 boundary. The disease signal in this
study is therefore carried by the **continuous** injury-program score (disease main figure), not by
a shift in how many pericytes fall into each discrete cluster — the continuous result is not a
relabelled abundance change.

### Figure S12 — `figureS_sensitivity` — Robustness and limitations of the disease-associated injury-stromal signal

**Figure S12.** What the donor-level disease association does and does not survive. **(A)** The
injury-stromal score by disease group under the primary composite (*AGTR1* excluded) and under a
sensitivity composite that adds the *AGTR1*-positive fraction; points are donors, diamonds are
marginal means ±95% CI. Excluding *AGTR1* avoids both dropout and circularity, and it does not
create the disease association — the two composites give the same picture. **(B)** The
vascular-stability component and the net niche-stability index. Stability **rises alongside**
injury in fibrotic/ILD lungs, which is precisely why the net index (stability minus injury) is a
weaker discriminator than its injury half; the net index is reported for completeness, not as the
primary endpoint. **(C)** Smoking status is recorded for 14 of 32 Healthy donors and for **no**
diseased donor, so a smoking-stratified disease contrast is not merely underpowered but
inestimable. **(D)** Leave-one-study-out: across all 16 refits the Fibrotic/ILD effect on the
injury-stromal score stays positive, and it remains significant in 13 of 16; significance is lost
only when one of three fibrosis-enriched cohorts is removed, which reflects loss of exposed donors
rather than dependence on a single cohort. The other three responses are non-significant across
every refit. **(E)** Among donors that do carry a smoking label, the injury/*AGTR1* readouts show
no smoking gradient (marginal means ±95% CI). **(F)** The Healthy-vs-Fibrotic/ILD effect is
essentially unchanged when smoking is added as a covariate. Donor counts here come from the
`sensitivity`/`niche_index` modules; the within-study meta-analysis in the disease main figure uses
a different donor set and endpoint and reports its own smoking availability (21 of 42 Healthy) —
the two should not be quoted for the same claim. HLCA carries no medication metadata; medication
sensitivity is noted as a limitation/future-cohort question (see Methods).

### Figure S2 — `figureS_crossspecies_mouse` — *Agtr1a* marks the mouse pericyte compartment

**Figure S.** The mouse lung data support a compartment-level, but not a state-level,
cross-species comparison; *Agtr1a* is genuinely transcribed by mouse pericytes in raw counts.
**(A)** Composition of the mouse mural compartment across the four integrated datasets (M1–M4);
point area is proportional to cell number, with counts labeled. Two structural features preclude
a state-level comparison with the human analysis: the compartment contains only **41 pericytes**
(3.6% of 1,144 mural cells; the remainder are 1,016 vascular-associated and 87 pulmonary-artery
smooth muscle cells, against 11,680 human pericytes), and cell type is almost perfectly aliased
with dataset — M3/M4 contribute vSMC only, M1/M2 pericytes and PA-SMC only — so the
pericyte-versus-vSMC contrast is absorbed by any dataset term and only pericyte-versus-PA-SMC is
estimable within dataset. **(B)** *Agtr1a* in raw log-normalized counts for every pericyte and
PA-SMC in the two datasets containing both (all cells shown; bars are medians; labels give the
detected fraction and *n*). Detection is 64% (16/25) and 94% (15/16) in pericytes versus **0/43
and 0/44** in PA-SMC (Fisher exact *P* = 1.4 × 10⁻⁹ and 3.0 × 10⁻¹³; Mantel–Haenszel stratified
by dataset *P* = 8.4 × 10⁻²⁰), so the contrast is not a batch effect. Across both datasets, 31 of
41 pericytes (75.6%) spanning 15 of 18 donors carry ≥1 *Agtr1a* UMI. **(C)** The undetected
pericytes reflect sequencing depth, not absence: the two datasets differ ~360-fold in depth
(median 291 versus 104,973 UMI per pericyte) and detection tracks depth across cells (Spearman
ρ = 0.43, *P* = 0.005); at full depth detection is 15/16 with a median of 928 *Agtr1a* UMI per
positive cell. Open circles denote cells with no *Agtr1a* UMI. **(D)** Fraction of *Agtr1a*-positive
cells computed on raw counts versus the dense scVI-denoised layer. Denoising reconstructs every
cell as non-zero, converting 0% detection in PA-SMC to 100%; detection fractions and correlations
computed on that layer are therefore uninformative, which is why all claims here use raw counts.
*Agtr1b* and *Agtr2* are absent from mouse pericytes entirely (0/41 cells, maximum 0 UMI), making
*Agtr1a* the sole angiotensin II receptor transcribed in this compartment. Note that all 41
pericytes derive from healthy animals, so these data do not test the injury-associated pericyte
loss or its losartan rescue. Source: `cross_species/_h/04.species_comparability.py`.

### `figure_basement_membrane` — Pericytes deposit a selective vascular basement membrane

**Figure.** Lung pericytes selectively express a collagen-IV/nidogen basement-membrane (BM)
module that is transcriptionally distinct from fibrillar extracellular matrix. All panels use
donor × cell-type pseudobulk as the unit of analysis (*n* = 2,329 units, 22 cell types, 220
healthy donors, ≥5 cells per unit); expression is log1p of the mean CP10K within each unit
(back-transformed with `expm1` before averaging, so no Jensen bias is carried into the
cross-cell-type comparison), standardized within dataset, and modeled as
`expr_z ~ cell type + mean log10 total counts + (1|donor) + (1|study)` with Benjamini–Hochberg
correction. **(A)** Dot heatmap of the 13 BM components across 22 cell types; fill is expression
z-scored **within gene** (so cell types are comparable across genes of very different absolute
abundance), dot size is the detection fraction. Genes are ordered by structural class (collagen
IV, laminins, linkers/proteoglycans); pericytes are bolded on the axis. **(B)** Selectivity per
gene: log₂ ratio of the pericyte marginal mean to the next-highest cell type, with the pericyte
rank among 22 cell types printed at each bar. Pericytes rank **first** for *COL4A1* (log₂ 1.09),
*COL4A2* (0.99), *COL18A1* (0.85), *LAMB1* (0.54), *NID1* (1.22) and *NID2* (0.34), but not for
the laminin α-chains — *LAMA4* ranks 3rd (highest in alveolar fibroblasts), *LAMA5* 9th and
*LAMA3* 17th (highest in AT1). *LAMA3* and *LAMA5* were **pre-specified** as negative controls
and behave as predicted, confirming the metric is not simply tracking abundance. **(C)** The
pre-specified primary endpoint, `BM − fibrillar` (mean of per-gene z across the BM panel minus
the same across a 10-gene fibrillar panel). Because any multiplicative cell-size or
capture-efficiency constant applies to both panels, it cancels in the difference; this
distinguishes *selective BM deposition* from *generalized matrix richness*. Points are marginal
means with 95% CI. Pericytes are strongly BM-shifted relative to every fibroblast population
(peribronchial +2.16, *P* = 5.2 × 10⁻³⁰¹; adventitial +1.91, *P* = 1.3 × 10⁻²⁹¹; myofibroblast
+1.78, *P* = 7.0 × 10⁻¹⁸²; alveolar +1.47, *P* = 8.6 × 10⁻²⁰¹) yet statistically
**indistinguishable from capillary endothelium** (aerocyte +0.004, *P* = 0.93; general capillary
−0.07, *P* = 0.098), with AT1 even further BM-shifted (−0.53, *P* = 2.8 × 10⁻³⁴) — i.e. BM
deposition is a shared property of the barrier-forming cells that build the vascular and alveolar
basement membranes. For this endpoint the study random-effect SD (0.130) did not exceed the
residual SD (0.245) and no gene was flagged study-dominated. **(D)** Consequence for the pericyte
state model: cluster × program relative enrichment, with the BM panel included as a sixth
program; stars mark the winning program per cluster. Leiden clustering runs on the
study-integrated `X_pca_harmony` embedding and is independent of the marker panels, so adding a
panel cannot move the clusters — only their annotation. All three clusters previously labelled
*fibroblast-like* (clusters 1, 3, 5; 4,200 cells, 36.0% of pericytes) become BM-dominant. In
cluster 1 every one of the five original programs was **negatively** enriched (fibroblast-like
−0.33, activated/migratory −0.57, inflammatory −0.62, synthetic/contractile −0.66,
vascular-stabilizing −0.77) while BM scored +0.45, so the former label was a least-negative
default rather than a positive identity. The reassignment survives removing *COL4A1* from the
fibroblast-like panel, and the laminin and linker sub-panels each win independently.
Source data: `basement_membrane/_m/stats_data/`, `basement_membrane/_m/state_gate_*.tsv`.

### Figure S10 — `figureS_ras_landscape` — The lung renin–angiotensin axis is distributed across cell types

**Figure S.** No lung cell type carries an autonomous angiotensin circuit. Donor × cell-type
pseudobulk, *n* = 4,376 units, 22 cell types, 417 donors (≥5 cells per unit). **(A)** Detection
fraction of the local renin–angiotensin machinery — substrate (*AGT*), proteases (*REN*, *ACE*,
*ACE2*, *CMA1*, *CTSG*, *CTSD*, *ENPEP*, *MME*), receptors (*AGTR1*, *AGTR2*, *LRP2*, *MAS1*) and
TGF-β comparators — across cell types, ordered by *AGT*. Colour is the detection fraction on a
square-root scale. The steps segregate by compartment: *AGT* is perivascular-stromal (vascular
smooth muscle 0.086, adventitial fibroblasts 0.036), *AGTR1* is overwhelmingly pericyte (0.342,
an order of magnitude above the next cell type), *ACE* is capillary-endothelial and myeloid
(aerocyte capillary 0.351, general capillary 0.144, alveolar macrophages 0.196), and the
renin-independent chymase route is mast-cell restricted (*CTSG* 0.109, *CMA1* 0.025). Renin is
effectively absent everywhere (maximum 0.023). **(B)** Number of circuit steps present per cell
type at a 0.05 donor-level detection threshold; fill marks whether a cell type carries a complete
autonomous circuit (substrate + an angiotensin II–generating protease + AT1R). **Zero of 22 cell
types** qualify, and no cell type carries more than **one** of the three requirements — the
segregation is complete. The axis therefore runs perivascular stroma → endothelial ACE or
mast-cell chymase → pericyte AT1R and requires a minimum of three distinct cell types. Pericytes additionally rank first for
*ACE2* (0.030) and second for *ENPEP* (0.094), so they are both the principal receiver and a site
of signal termination. Caveat: *AGT* transcript is a proxy for **local** angiotensinogen only —
the dominant physiological source is hepatic and is invisible to lung single-cell data — so a low
local *AGT* signal does not imply low local angiotensin II.
Source data: `agt_axis/_m/stats_data/ras_*.tsv`.

### Figure S15 — `figureS_bm_copd` — Basement-membrane remodeling in IPF but not detectably in COPD

**Figure S.** Disease contrast for basement-membrane components in GSE136831 (Adams/Kaminski),
the only cohort in this project with a COPD arm. Donor × compartment pseudobulk built from raw
counts (summed within unit → CP10K → log1p), ≥5 cells per donor, adjusted for depth, sex, age and
ever-smoker status (available in this cohort but not in HLCA). **(A)** COPD-vs-Control effects for
the two **pre-specified** primary genes, *LAMB1* and *LAMA4*, across the seven compartments with
≥5 donors on both arms (ATI, ATII, endothelial, fibroblast, mural, myofibroblast, SMC); 14 tests,
Benjamini–Hochberg corrected within the primary family. Points are estimates with 95% CI. **No
contrast survives correction** (all BH ≥ 0.81; nearest myofibroblast *LAMA4*, nominal *P* = 0.058),
and no exploratory BM gene reaches BH < 0.05. The null is interpretable rather than a sensitivity
failure because the internal positive control fires: IPF-vs-Control shows substantial BM
remodeling in the same compartments and models (nominal *P*: endothelial *HSPG2* 1.5 × 10⁻⁵,
endothelial *LAMA3* 1.0 × 10⁻⁴, fibroblast *COL4A1* 1.9 × 10⁻³, fibroblast *COL18A1* 2.1 × 10⁻³,
mural *COL18A1* 1.8 × 10⁻³, mural *LAMA4* 1.2 × 10⁻²; IPF p-values are uncorrected, as IPF was
designated a positive control rather than a hypothesis family). **(B)** Minimum detectable effect
for the pericyte compartment at 80% power. **Pericytes are reported descriptively with no
p-value**: at a 5-cell floor this cohort has 6 COPD but only **one** Control donor with ≥5
pericytes, so a pericyte-specific COPD-vs-Control contrast is not estimable and the MDE is
undefined. No claim is made about pericyte-specific BM dysregulation in COPD in either direction.
Source data: `basement_membrane/_m/stats_data/bm_copd_*.tsv`, `bm_pericyte_power.tsv`.

### Figure S16 — `figureS_disease_robustness` — The Fibrotic/ILD injury-program effect survives dropping any cohort, and is reproduced *within* studies

**Figure S16.** Study-level robustness of the Fibrotic/ILD − Healthy difference in the
donor-level pericyte injury-program score reported in `figure_disease_main` **A**. The two
panels answer two **different** objections and neither substitutes for the other, which is why
they travel together: **(A)** asks *does one cohort create the effect*, **(B)** asks *is the
effect a between-study batch artifact at all*. LOSO cannot answer the second (dropping a study
still leaves an estimate pooled across the rest), and the forest cannot answer the first (it
covers only the three studies that sampled both arms).

**(A) Leave-one-study-out.** The panel-A three-group model refit 23 times, dropping one dataset
each time (blue points, 95 % CIs; orange line = the full-data estimate, +0.79 SD). The estimate
stays between **+0.65 and +1.08 SD and remains significant in 23 of 23 refits**, so no single
cohort creates the association. Source:
`disease_association/_m/mixed_model_forest/leave_one_study_out_3group.tsv`.

**(B) Within-study random-effects meta-analysis.** Each blue point is the contrast estimated
**within a single study** that sampled both groups (≥ 2 donors per group after a
≥ 10-pericyte-per-donor filter); point area is proportional to inverse-variance weight and
whiskers are 95 % CIs. The orange diamond is the DerSimonian–Laird random-effects pooled
estimate: **+0.90 SD (95 % CI 0.47–1.33), *P* = 4.5×10⁻⁵, I² = 25 %** (low heterogeneity),
carried by the two balanced cohorts (Banovich_Kropski_2020, Sheppard_2020). Because every
contrast is estimated inside one study, the effect cannot be a between-study batch artifact.
This panel was **panel A of `figure_disease_main`** until the continuous-injury rebuild
(2026-07-28). Source:
`disease_association/_m/mixed_model_forest/forest_per_study.tsv`, `forest_pooled_RE.tsv`.

**Scale note — A and B are not on a common scale.** A refits the **three-group** model, whose
program scores are *z*-standardized over Healthy + Fibrotic/ILD + Other; B is the **two-group**
primary model, standardized over Healthy + Fibrotic/ILD only. So **+0.90 SD in B and +0.79 SD
in A (or in `figure_disease_main` A) are not the same estimate on the same scale** and must not
be quoted interchangeably; each axis title names its own standardization set, and no cross-panel
annotation is drawn. The agreement that matters here is **directional and in significance**, not
numeric.
