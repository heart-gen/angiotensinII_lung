# gwas_enrichment — human genetic support for the AGTR1⁺ injury-niche pericyte state

Tests whether **genetic risk for lung disease concentrates in the cell states**
defined by the rest of the study, using two orthogonal methods plus locus-level
colocalization. This converts "AGTR1 marks an injury-niche pericyte state" into a
human-genetics statement: *heritability for COPD/IPF/lung-function is enriched in the
genes that define that state, and the AGTR1 expression signal colocalizes with COPD risk.*

## Why this addresses the "descriptive" critique
Independent (germline genetic) evidence that the transcriptomic state is disease-relevant —
not derivable from expression alone. Two methods (MAGMA gene-property + S-LDSC) that must
agree, plus a focal AGTR1 colocalization, is the standard top-tier package.

## The AGTR2 / chromosome-X constraint (important)
The LD reference (1000G Phase3 EUR), baselineLD model, and HM3 weights are **autosomal**
(MHC excluded). So:
- **AGTR1 (chr3)** — fully testable; single-gene coloc is AGTR1-only.
- **AGTR2 (Xq23)** — *excluded* from every SNP-heritability/gene-level test. Its evidence
  comes from the expression / `pathway_balance` modules, not genetics.
- The headline test is **cell-state gene-program** enrichment (not a single gene), so it is
  robust to AGTR2 being unavailable. State specificity is built on **autosomal protein-coding
  genes only** (`00.build_specificity.py` drops chrX/Y).

## Pipeline (all hg19 / GRCh37: g1000_eur + 1000G Phase3 panels; rsID-keyed)

| Step | Script | What | Where to run |
|---|---|---|---|
| stage | `inputs/gwas/_h/00.stage_gwas.sh` | copy sumstats → `inputs/gwas/_m/` + `SOURCES.md`; optional Shrine download | **login node** |
| 0 | `step_0.sh` → `00.build_specificity.py` | autosomal cell-type + pericyte-state specificity; LDSC BEDs | RM-shared |
| 1 | `step_1.sh` → `01.munge_gwas.py` | harmonize GWAS → MAGMA (`SNP,P,N`) + LDSC (`SNP,A1,A2,N,Z`) | RM-shared |
| 2 | `02.run_magma.sh` | MAGMA gene-property (marginal + joint-states) | RM-shared |
| 3 | `03.make_annot.sh` | S-LDSC binary annotations per chr | RM-shared |
| 4 | `04.compute_ldscores.sh` | LD scores (array 1–22) | RM-shared array |
| 5 | `05.partitioned_h2.sh` | munge + partitioned h² per trait (array 0–12) | RM-shared array |
| 6 | `step_6.sh` → `06.coloc_agtr1.R` | AGTR1 coloc: COPD GWAS × lung cis-eQTL | RM-shared |
| 7 | `step_7.sh` → `07.summarize_and_figures.R` | collate + heatmap/forest/concordance | RM-shared |

```bash
bash inputs/gwas/_h/00.stage_gwas.sh           # login node
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/step_0.sh
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/step_1.sh
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/02.run_magma.sh
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/03.make_annot.sh
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/04.compute_ldscores.sh   # after 03
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/05.partitioned_h2.sh     # after 04
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/step_6.sh
sbatch -D gwas_enrichment/_m gwas_enrichment/_h/step_7.sh                # after 02 & 05 (& 06)
```

One-time R dependency (login node): `coloc` into the project `.Rlib`
```r
install.packages("coloc", lib=".Rlib", repos="https://cloud.r-project.org")
```

## Traits
- **Lung disease:** COPD (COPDGene spirometry + UKB), IPF (UKB), asthma (TAGC), lung function
  (Shrine FEV1/FVC — external download).
- **Smoking** (GSCAN SmkInit/CigDay): confound check *and* biology.
- **RAS/cardiovascular** (ICBP SBP/DBP, CARDIoGRAM CAD): AT1R-axis relevance.
- **Negative controls** (educational attainment, body-fat %): specificity guards.

COPDGene reports only |β| (unsigned) → MAGMA + coloc only (no S-LDSC; coloc.abf is
sign-agnostic). All other traits run all applicable methods.

## Tools / resources (verified on PSC)
- MAGMA v1.10 + `g1000_eur` + `NCBI37.3.gene.loc`: `/ocean/projects/bio250020p/shared/opt/magma-v1.10/`
- LDSC (Py3-ported) + `ldsc_wrapper.py`: `/ocean/projects/bio250020p/shared/opt/ldsc/`, env `…/opt/env/genomics`
- LDSC reference panels: `/ocean/projects/bio250020p/shared/resources/ldsc/` (1000G Phase3, baselineLD v2.2, HM3 weights)
- GTEx v11 lung cis-eQTL: `/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11/GTEx_Analysis_v11_eQTL/`

## Caveats to state in Methods
European-ancestry discovery + LD reference; AGTR2/X excluded from genetics; cell-state
enrichment reflects the gene program (not AGTR1 alone); negative-control traits included to
demonstrate specificity.
