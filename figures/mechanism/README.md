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
expressing), identifying which neighboring populations supply each signal. **(D)** Cross-program
specificity by projection: pericyte-learned CoGAPS expression patterns (de-novo-selected nP = 8,
one row per pattern labelled by its assigned program) projected (projectR/OLS) onto a cell-type ×
donor pseudobulk; fill = within-pattern *z* across cell types. The patterns that rank pericytes
first are the **vascular-stabilizing and the two basement-membrane patterns** (the structural/
support programs), whereas the fibrillar-matrix and inflammatory patterns are mirrored most
strongly in bona-fide fibroblasts — the injury axes a subset of pericytes adopt are the same ones
fibroblasts constitutively run. Abbreviations: CCC, cell–cell communication; AUPR, area under the
precision–recall curve.

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

### `figureS_alluvial` — Stable subclusters → dominant program → effector class

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

### `figureS_program_category` — Program × protein-category enrichment

**Figure S.** Dot heatmap relating the six pericyte programs (rows) to eight curated marker
protein categories (columns: signaling ligands; chemokines; cytokines; adhesion molecules; ECM
structural; matrix-remodeling; fibrotic mediators; mural identity). Dot size = prevalence (mean
fraction of the program's cells expressing the category); fill = relative enrichment (*z* of
mean detection across programs, computed within category to expose program-specific signal
otherwise swamped by category-level baseline detectability). Cells are assigned to their dominant
program by *z*-scored argmax of the six program scores — now including basement-membrane
(`cell_communication/_h/03b.program_category_enrichment.py`), retaining all six programs.

### `figureS_pericyte_layer` — Pericyte-layer supporting detail

**Figure S.** Companion detail for `figure_pericyte_layer`, all on the same shared UMAP
(*n* = 11,680 pericytes). **(A)** Per-program module-score overlays (vascular-stabilizing,
synthetic/contractile, activated/migratory, inflammatory, fibroblast-like, basement-membrane),
showing the spatial organization of each of the six programs across the embedding. **(B)** *ACTA2* expression (log), the canonical
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

### `figureS_crossspecies_mouse` — *Agtr1a* marks the mouse pericyte compartment

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

### `figureS_ras_landscape` — The lung renin–angiotensin axis is distributed across cell types

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

### `figureS_bm_copd` — Basement-membrane remodeling in IPF but not detectably in COPD

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
