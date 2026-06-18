# Mechanism figures — index and figure legends

Publication-ready figures for the mechanistic revision (target: *Circulation Research*).
All panels are exported as `.pdf` (vector, cairo), `.svg`, and `.png` (350 dpi) by
`figures/_h/manuscript_mechanism_figure.R` and `figures/_h/sensitivity_robustness_figure.R`
(run via `figures/_h/step_figures.sh` from `figures/_m`). Shared visual language:
Okabe–Ito palette, `theme_ms` (8–9 pt base), no in-panel titles (interpretation lives in
the legends below), bold uppercase panel tags.

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

**Supplements:** `figureS_pericyte_layer`, `figureS_acta2_control`, `figureS_alluvial`,
`figureS_program_category`, `figureS_balance_by_state`, `figureS_sensitivity`,
`figureS_crossspecies_mouse`.

**Supporting files (not figures):**
- `tableS_acta2_control.tsv` — source-data values underlying `figureS_acta2_control`.
- `figure_panel_manifest.tsv` — maps the per-module PDFs (liana dotplots, NicheNet heatmap,
  state/DPT/PAGA UMAPs, etc.) that are assembled into the main figures in a vector editor;
  `exists` flags whether each source panel has been generated.
- `figureB_states_continuum_niche.{pdf,png}` — **working draft only** (legacy
  `assemble_mechanism_figures.R`, `theme_bw`, in-panel titles). Superseded for submission by
  `figure_mechanism_main` + the manifest image panels; kept for quick QC, not for the manuscript.

**State-model key (applies to every panel that mentions programs).** Six stable pericyte
Leiden subclusters collapse onto three discrete dominant programs
(`pericyte_states/_m/annotations/state_program_map.tsv`): vascular-stabilizing (clusters 0, 2),
fibroblast-like (clusters 1, 3, 5), activated/migratory (cluster 4). "Fibroblast-like" is a
pericyte program resembling fibroblasts, **not** bona-fide fibroblasts. Two further program
*scores* (inflammatory, synthetic/contractile) are continuous axes that are never a cell's
dominant discrete program, so they appear in score/continuum panels but not in discrete
program-assignment panels.

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
dominant state-program (vascular-stabilizing, fibroblast-like, activated/migratory). **(D)**
*AGTR1* across programs under three measurement lenses — raw expression, binary detection, and
scVI-denoised — as donor-aware marginal means centered within lens (±SE). The raw and detection
lenses share a vascular-stabilizing–high pattern that **reverses** under denoising (denoised
*AGTR1* is flat-to-injury), identifying the raw enrichment as a transcript-capture/dropout effect
and establishing that *AGTR1* is not a vascular-stabilizing state marker. **(E)** Diffusion
pseudotime on the same embedding, ordering cells along a vascular-stabilizing ↔ injury/fibrotic
continuum (ordering reflects transcriptional similarity, not time). **(F)** Donor-level Spearman
correlation between each program score / *AGTR1* and pseudotime; orange = *P* < 0.05. Program
scores rise monotonically along the continuum (donor ρ ≈ 0.46–0.52); *AGTR1* shows a weak positive
trend (ρ ≈ 0.25). Together: where *AGTR1* sits (A–C), why it is a compartment label not a state
marker (D), and how the compartment is organized as an injury continuum (E–F).

### `figure_ccc_nichenet` — Niche signaling drives pericyte target programs

**Figure.** Cell–cell communication and ligand–target inference identify the niche signals
received by lung pericytes. **(A)** NicheNet ligand-activity ranking (top 15 ligands) for
predicting the pericyte target-gene program; bars show AUPR (corrected). TGF-β–family ligands
(*TGFB1/2/3*, *CCN1*, *CCN2*) are highlighted (orange). **(B)** Ligand→target regulatory-potential
heatmap for the prioritized ligands (rows) and their top predicted pericyte target genes
(columns); fill encodes NicheNet regulatory potential. **(C)** Dot plot of the fraction of cells
expressing each prioritized ligand across the niche sender cell types (size and color = fraction
expressing), identifying which neighboring populations supply each signal. **(D)** Cross-program
specificity by projection: pericyte-learned CoGAPS expression programs projected (projectR/OLS)
onto a cell-type × donor pseudobulk; fill = within-program *z* across cell types. The
fibroblast-like injury program is mirrored most strongly in bona-fide fibroblasts, whereas the
vascular-stabilizing and synthetic/contractile programs remain mural-compartment specific.
Abbreviations: CCC, cell–cell communication; AUPR, area under the precision–recall curve.

### `figure_mechanism_main` — Donor-level disease phenotype and the stabilizing↔injury continuum

**Figure.** Pericyte niche state across health and disease at the donor level. **(A)**
Niche-stability index and **(B)** injury-stromal score per donor, grouped by disease
(Healthy, COPD, Fibrotic/ILD); box = median and IQR, whiskers = 1.5×IQR, white diamond = mean,
points = individual donors; brackets show Wilcoxon rank-sum comparisons versus Healthy. **(C)**
AT1R–AT2R pathway balance per donor (shown as a **corollary** of injury intensity: the disease
effect is redundant with the injury-stromal score and collapses after adjustment — see
`MECHANISM_ANALYSES.md`). **(D)** Mean donor-level pericyte state-program composition by disease
(stacked fractions across the three discrete programs plus the continuous-program assignments).
**(E)** Continuum trends: donor-level Spearman correlation between diffusion-pseudotime ordering
and each program score / *AGTR1*; points colored by significance (orange, *P* < 0.05). The
ordering reflects transcriptional similarity along a vascular-stabilizing ↔ injury/fibrotic
**continuum**, not a temporal axis. Vascular-stabilizing, inflammatory, synthetic/contractile,
and activated/migratory scores rise monotonically along the continuum (donor ρ ≈ 0.46–0.52,
*P* < 1×10⁻¹¹); *AGTR1* shows a weak positive trend (ρ ≈ 0.25). Assemble with the state UMAP,
DPT-pseudotime UMAP, and PAGA panels in `figure_panel_manifest.tsv`. *n* = 32 donors (donor-level
mixed-model marginal means; df = 31).

### `figureS_alluvial` — Stable subclusters → dominant program → effector class

**Figure S.** Alluvial flow linking the six stable pericyte Leiden subclusters (left axis) to
their three dominant programs (middle axis: vascular-stabilizing, fibroblast-like,
activated/migratory) to six functional effector molecule classes (right axis: signaling ligands;
chemokines/cytokines; adhesion molecules; ECM structural; matrix-remodeling enzymes; fibrotic
mediators). Left→middle ribbon width is proportional to subcluster cell count; middle→right width
is each program's effector composition (donor-averaged mean fraction of cells expressing the
class, normalized within program so flow is conserved). Ribbons are colored by **stable
subcluster**, with hue families grouped by program (blues = vascular-stabilizing [P0, P2];
oranges = fibroblast-like [P1, P3, P5]; pink = activated/migratory [P4]) so both the
subcluster→program merge and the program grouping are legible. Cluster→program assignments are
the canonical `state_program_map.tsv` used throughout the analysis.

### `figureS_program_category` — Program × protein-category enrichment

**Figure S.** Dot heatmap relating the five pericyte programs (rows) to eight curated marker
protein categories (columns: signaling ligands; chemokines; cytokines; adhesion molecules; ECM
structural; matrix-remodeling; fibrotic mediators; mural identity). Dot size = prevalence (mean
fraction of the program's cells expressing the category); fill = relative enrichment (*z* of
mean detection across programs, computed within category to expose program-specific signal
otherwise swamped by category-level baseline detectability). Cells are assigned to their dominant
program by *z*-scored argmax of the five program scores
(`cell_communication/_h/03b.program_category_enrichment.py`), retaining all five programs.

### `figureS_pericyte_layer` — Pericyte-layer supporting detail

**Figure S.** Companion detail for `figure_pericyte_layer`, all on the same shared UMAP
(*n* = 11,680 pericytes). **(A)** Per-program module-score overlays (vascular-stabilizing,
synthetic/contractile, activated/migratory, inflammatory, fibroblast-like), showing the spatial
organization of each program across the embedding. **(B)** *ACTA2* expression (log), the canonical
contractile-mural benchmark, which mirrors the raw *AGTR1* vascular-stabilizing–high pattern.
**(C)** Binary *AGTR1* detection, visualizing the dropout structure that drives the raw
enrichment reversed in `figure_pericyte_layer` panel D. **(D)** Donor-mean *AGTR1* versus *ACTA2*
expression by stable subcluster (donor × subcluster means, subclusters with ≥5 cells; box =
median/IQR, whiskers = 1.5×IQR), comparing the two markers' subcluster profiles.

### `figureS_acta2_control` — AGTR1 is not reducible to ACTA2⁺ contractile identity

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

### `figureS_balance_by_state` — AT1R–AT2R balance by pericyte program

**Figure S.** Donor-level AT1R–AT2R pathway balance across pericyte programs (box = median/IQR,
points = donors). Differences among programs are **not significant** (smallest pairwise
*P* = 0.067); this panel documents that the AT1R/AT2R skew tracks continuous injury intensity and
disease rather than discrete program identity, justifying its demotion from the main figure.

### `figureS_sensitivity` — Smoking and cohort robustness

**Figure S.** Robustness of the donor-level disease associations. **(A)** Smoking status is
recorded essentially only for Healthy donors, making a smoking-stratified disease contrast
inestimable (the underlying confound). **(B)** Among donors carrying a smoking label, the
injury/*AGTR1* readouts show no smoking gradient (marginal means ±95% CI), so the phenotype is
not a smoking artifact. **(C)** The Healthy-vs-Fibrotic/ILD disease effect is essentially
unchanged when smoking is added as a covariate (base vs +smoking marginal means ±95% CI).
**(D)** Leave-one-study-out analysis: the Fibrotic/ILD effect on the headline injury-state
fraction is stable across cohorts. HLCA carries no medication metadata; medication sensitivity is
noted as a limitation/future-cohort question (see Methods).

### `figureS_crossspecies_mouse` — Cross-species conservation of the AGTR1 mural signal

**Figure S.** Conservation of the *AGTR1/Agtr1a* mural-compartment signal in mouse lung. **(A)**
Per-cell scatter of *Agtr1a* (scVI batch-corrected) versus the vascular-stabilizing program
score in mouse mural cells, with per-dataset linear fits (color = dataset D1–D4). **(B)**
Per-dataset Spearman ρ between *Agtr1a* and the vascular-stabilizing score (*n* labeled per bar);
the positive relationship is consistent across all four mouse datasets (ρ ≈ 0.53–0.89). Because
the mouse mural set is sparse and ~91% homeostatic-contractile, conservation is reported at the
**compartment level** (*Agtr1a* labels the mural compartment), not as a discrete
stabilizing-vs-injury state axis, consistent with the human compartment-level framing and the
wet-lab pericyte-loss / losartan-rescue phenotype.
