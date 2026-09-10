# Cell type-specific expression of angiotensin receptors in the human lung

## Abstract
The renin-angiotensin system is a highly characterized integrative pathway in mammalian homeostasis whose clinical spectrum has been expanded to lung disorders such as chronic obstructive pulmonary disease (COPD)-emphysema, idiopathic pulmonary fibrosis (IPF), and COVID pathogenesis. Despite this widespread interest, specific localization of this receptor family in the mammalian lung is limited, partially due to the imprecision of available antibody reagents. In this study, we establish the expression pattern of the two predominant angiotensin receptors in the human lung, *AGTR1* and *AGTR2*, using complementary and comprehensive bulk and single-cell RNA-sequence datasets that are publicly available. We show these two receptors have distinct localization patterns and developmental trajectories in the human lung, pericytes for *AGTR1* and a subtype of alveolar epithelial type 2 cells for *AGTR2*. In the context of disease, we further pinpoint *AGTR2* localization to the COPD-associated subpopulation of alveolar epithelial type 2 (AT2B) and *AGTR1* localization to fibroblasts, where their expression is upregulated in individuals with COPD, but not in individuals with IPF. Finally, we examine the genetic variation of the angiotensin receptors, finding *AGTR2* associated with lung phenotype (i.e., cystic fibrosis) via rs1403543. Together, our findings provide a critical foundation for delineating this pathway’s role in lung homeostasis and constructing rational approaches for targeting specific lung disorders.


---

## Current state of this repository (as of 2026-09-07)

**The abstract above describes the 2024 preprint and is left unchanged.** The
repository has since grown a second, larger body of work — a mechanistic
single-cell analysis of lung pericytes — and **that work revises one of the
preprint's conclusions.** This section exists so a reader lands on the current
position rather than inferring it.

### What changed relative to the preprint

| Preprint claim | Current position |
| -------------- | ---------------- |
| *AGTR1* localizes to pericytes | **Held.** Confirmed and extended: *AGTR1* labels the pericyte/mural **compartment**. |
| *AGTR1* marks a distinct fibroblast population upregulated in COPD | ⚠️ **Not supported.** Six of seven stromal cell types show *negative* Fibrotic−Healthy estimates in the HLCA, nothing survives BH correction, and the independent GSE136831 evaluation disagrees in sign against an HLCA estimate that is itself null. **The analyses behind this now live in [`heart-gen/lung-pericyte-analysis`](https://github.com/heart-gen/lung-pericyte-analysis)** — see `writings/AGTR1_DISEASE_DIRECTION.md` there, which is the document that reconciles the two. |
| — | **New:** *AGTR1* is a **compartment label, not a state marker**. The raw enrichment reverses under denoising and under a no-imputation negative-binomial count model. |
| — | **New:** disease acts on pericytes through **continuous injury-program intensity**, not through discrete state composition. Established in [`heart-gen/lung-pericyte-analysis`](https://github.com/heart-gen/lung-pericyte-analysis), not here. |

### Repository map

**Mechanistic pipeline** (2026; each module is `<name>/_h` scripts → `<name>/_m` outputs)

| Module | What it does |
| ------ | ------------ |
| `pericyte_states/` | Six stable pericyte states, three programs, and the diffusion-pseudotime axis. **Root of the pipeline** — seven modules consume its labels. |
| `basement_membrane/` | Basement-membrane vs fibrillar-collagen matrix programs; NicheNet ligand→matrix inference. |
| `niche_index/` | Donor-level composite: stability minus injury. |
| `cell_communication/` | LIANA + NicheNet signalling into pericytes, with a permutation specificity control. |
| `pathway_balance/` | AT1R/AT2R transcriptional balance. |
| `pericyte_cogaps/` | Unsupervised NMF (CoGAPS) validation of the state model. |
| `agt_axis/` | Whether the AGT→AGTR1 edge is coherent (it is a three-cell relay). |
| `cross_species/` | Mouse lung mural comparison. |
| `localization/airspace_analysis/` | Airspace proximity and the *AGTR1* dropout audit. ⚠️ Current work, not legacy. |

**Support**: `inputs/` (data acquisition and the shared HLCA/GEO object builders),
`figures/` (assembly), `tables/` (supplementary tables), `writings/`
(manuscript-voice summaries and the defect ledger — **not tracked**, see below).

**Not here (moved 2026-09-10)**: the disease layer — `disease_association/`,
`sensitivity/`, the niche-index disease test, Figure 5, supplements S12/S16 and
Tables S13A/S13E/S13F/S14 — is in
[`heart-gen/lung-pericyte-analysis`](https://github.com/heart-gen/lung-pericyte-analysis).
**This repository makes no disease claims.** Two shared builders that lived under
`disease_association/` despite not being disease analyses stayed and were
relocated: `inputs/hlca/_h/02.preprocess_reference.py` (read by
`cell_communication/`) and `inputs/ipf/_h/01.generate_ipf_data.R` (read by the
`basement_membrane/` COPD arm). See `disease_association/MOVED.md`.

**Legacy (2024 preprint)**: `localization/{pericyte_analysis,lungmap_replication,enrichment_analysis}/`.
The retired 2024 disease cohorts (`copd/`, `ipf_analysis/`) went to the disease
repository with the rest of `disease_association/`; the one script that was
exempt from their retirement, the GSE136831 h5ad builder, is now
`inputs/ipf/_h/01.generate_ipf_data.R`.

### Running it

```bash
bash submit_pipeline.sh          # submits the module DAG to SLURM (Bridges-2)
```

Each module is also independently runnable as `cd <module>/_m && sbatch ../_h/step_N.sh`.

⚠️ `inputs/hlca/_h/step_2.sh` builds the all-cell HLCA object
`cell_communication/` reads. It is a 12-hour EM job and is **not** in
`submit_pipeline.sh`'s dependency chain — run it once, first, if
`inputs/hlca/_m/hlca_full.dataset.h5ad` is absent.

### A note on `writings/`

`writings/` is **gitignored** and contains the working documents: manuscript-voice
`*_SUMMARY.md` files, the per-analysis PI briefings in `writings/pi_briefings/`,
and `writings/TODO.md` — a defect ledger recording every known open issue with its
evidence. If you are picking this repository up, **`writings/TODO.md` is the file
to read first**; it is the authority on what is and is not currently trustworthy.


## Citation

If you the anything in this repository please cite the following pre-print:

Benjamin, K J M, Sauler, M, Poonyagariyagorn, H, and Enid R Neptune. "Cell type-specific expression of angiotensin receptors in the human lung with implications for health, aging, and chronic disease". *BioRxiv*. 2024. PMID: [38948835](https://biorxiv.org/cgi/content/short/2024.06.17.599425v1).

## License

<img src="https://licensebuttons.net/l/by-nc/3.0/88x31.png" alt width="88" height="31" scale="0">
Attribution-NonCommercial: CC BY-NC

This license lets others remix, tweak, and build upon our work non-commercially as long as they acknowledge our work.

[View License Deed](https://creativecommons.org/licenses/by-nc/4.0) | [View Legal Code](https://creativecommons.org/licenses/by-nc/4.0/legalcode)

