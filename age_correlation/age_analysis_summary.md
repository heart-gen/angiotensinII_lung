# Age-Associated Analysis Summary: Angiotensin Receptors in Human Lung

## Overview

This document summarizes the age-associated changes in angiotensin receptor expression (AGTR1, AGTR2) and cell type proportions in human lung tissue, combining evidence from snRNA-seq data, bulk RNA-seq, and deconvolution analyses.

---

## Data Sources

| Analysis Type | Data Source | Location |
|---------------|-------------|----------|
| snRNA-seq | HLCA v2 (Human Lung Cell Atlas) | `age_correlation/_h/` |
| Bulk RNA-seq | GTEx v8 lung tissue | `age_correlation/bulk_correlation/` |
| Deconvolution | BayesPrism on GTEx | `deconvolution/` |

---

## 1. Gene Expression Changes with Age

### AGTR1 (Angiotensin Receptor Type 1)

| Method | Direction | Estimate | P-value | Significance |
|--------|-----------|----------|---------|--------------|
| Bulk GTEx | ↑ Increases | 0.0019 | 0.005 | **Significant** |
| Deconvolution-adjusted | ↑ Increases | 0.0065 | 0.004 | **Significant** |
| snRNA-seq (smooth muscle) | ↑ Trending | rho=0.36 | 0.057 | Marginal |

### AGTR2 (Angiotensin Receptor Type 2)

| Method | Direction | Estimate | P-value | Significance |
|--------|-----------|----------|---------|--------------|
| Bulk GTEx | ↓ Decreases | -0.0096 | 1.3e-8 | **Highly significant** |
| Deconvolution-adjusted | ↓ Decreases | -0.071 | 1.1e-9 | **Highly significant** |

**Summary**: AGTR1 modestly increases with age while AGTR2 strongly decreases with age in bulk lung tissue.

---

## 2. Cell Types with Highest AGTR1 Expression (snRNA-seq)

| Rank | Cell Type | Mean Expression | % Positive Cells | N Donors |
|------|-----------|-----------------|------------------|----------|
| 1 | Pericytes | 0.725 | 44.8% | 11 |
| 2 | Alveolar fibroblasts | 0.149 | 12.5% | 22 |
| 3 | Adventitial fibroblasts | 0.100 | 13.4% | 28 |
| 4 | Smooth muscle | 0.065 | 5.4% | 29 |
| 5 | Peribronchial fibroblasts | 0.049 | 5.1% | 17 |

---

## 3. Cell Type Proportion Changes with Age (Deconvolution)

### Significant Changes

| Cell Type | R (correlation) | P-value | Direction | Interpretation |
|-----------|-----------------|---------|-----------|----------------|
| **Alveolar Fibroblasts** | 0.22 | 1.1e-07 | ↑ Increase | Strongest expansion |
| **Pericytes** | 0.15 | 0.00029 | ↑ Increase | Moderate expansion |
| **Smooth Muscle** | 0.11 | 0.0074 | ↑ Increase | Modest expansion |

### Non-Significant Changes

| Cell Type | R (correlation) | P-value | Direction |
|-----------|-----------------|---------|-----------|
| Adventitial Fibroblasts | 0.067 | 0.11 | NS |
| Peribronchial Fibroblasts | 0.018 | 0.66 | NS |

**Key Finding**: No pericyte depletion observed. Both pericytes and alveolar fibroblasts expand with age, but alveolar fibroblasts show the strongest effect (~1.5x stronger correlation than pericytes).

---

## 4. Cell-Type-Specific AGTR1 Age Correlations (snRNA-seq)

| Cell Type | Spearman rho | P-value | FDR | N Donors |
|-----------|--------------|---------|-----|----------|
| Smooth muscle | 0.357 | 0.057 | 0.285 | 29 |
| Alveolar fibroblasts | 0.092 | 0.685 | 0.947 | 22 |
| Pericytes | -0.082 | 0.818 | 0.947 | 11 |
| Peribronchial fibroblasts | -0.058 | 0.825 | 0.947 | 17 |
| Adventitial fibroblasts | -0.013 | 0.947 | 0.947 | 28 |

**Note**: No cell type reaches FDR significance for within-cell-type AGTR1 expression changes with age. Smooth muscle shows a trend (p=0.057) but limited statistical power.

---

## 5. Summary: Where is the Age Association Seen?

| Finding | snRNA-seq | Bulk RNA-seq | Deconvolution |
|---------|-----------|--------------|---------------|
| AGTR1 increases with age | Marginal (smooth muscle only, p=0.057) | **Yes** (p=0.005) | **Yes** (p=0.004) |
| AGTR2 decreases with age | Not tested | **Yes** (p<1e-8) | **Yes** (p<1e-9) |
| Alveolar fibroblasts expand | Not tested | N/A | **Yes** (R=0.22, p<1e-7) |
| Pericytes expand | Not tested | N/A | **Yes** (R=0.15, p=0.0003) |
| Pericyte depletion | Not tested | N/A | **No** (opposite: expansion) |

---

## 6. Interpretation

1. **Bulk tissue AGTR1 increase with age** is supported by both direct bulk analysis and deconvolution-adjusted analysis.

2. **The mechanism appears to be compositional**: Alveolar fibroblasts and pericytes (the top AGTR1-expressing cell types) both expand with age, driving increased bulk AGTR1 levels.

3. **Alveolar fibroblast expansion is the dominant effect**: With R=0.22 vs R=0.15 for pericytes, fibroblast expansion contributes more to the age-associated AGTR1 increase.

4. **Within-cell-type expression changes are weak**: snRNA-seq data shows no significant cell-type-specific AGTR1 expression changes with age after FDR correction, suggesting the bulk effect is primarily driven by cell composition shifts rather than per-cell expression changes.

5. **AGTR2 shows opposite pattern**: Strong decrease with age in bulk tissue, potentially reflecting different cellular sources or regulatory mechanisms.

---

## Methods Summary

- **snRNA-seq**: Spearman correlation of donor-level mean AGTR1 expression vs age; GAM models for non-linear effects
- **Bulk RNA-seq**: Linear regression with TMM + voom normalization, adjusting for sex and race
- **Deconvolution**: BayesPrism reference-based deconvolution using HLCA v2 as reference; Pearson correlation of estimated cell proportions vs age
