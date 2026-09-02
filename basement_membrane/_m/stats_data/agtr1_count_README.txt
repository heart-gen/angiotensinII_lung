AGTR1 count-model lens -- no imputation
============================================================

AGTR1 integer counts (raw/X, not the soupX float matrix) are the RESPONSE of
a negative-binomial GLMM with offset(log(library size)) and (1|study) +
(1|donor_id). Dropout is therefore inside the likelihood, and nothing is
borrowed across genes, so this lens cannot manufacture or erase an
association the way a denoiser can.

ROLE REVERSAL: the matrix score is the predictor here and the outcome in
bm_vs_agtr1_models.tsv. The null is the same; the sign means 'one SD higher
matrix score goes with higher AGTR1 per transcript sampled'. This is an
association test, not a direction-of-effect test.

Depth appears twice by design: the offset fixes the AGTR1 exposure, and
mean_log10_counts is free because the matrix score is itself depth-dependent.

Pseudobulk units: 214 (>= 5 cells) from 95 donors.
Cells: 11680 from 194 donors; AGTR1 detected in 37.4%.

Rows marked model = 'Poisson+OLRE' are fallbacks where glmer.nb did not
converge; they are not NB fits and should be reported as such.
BH is applied within (level, model).

GROUP COMPARISONS -- two groupings, two files, different admissibility:
  agtr1_count_by_cluster{,_posthoc}.tsv  grouping = pericyte_state (Leiden).
    AGTR1 is not among the 2,000 HVGs that define these clusters, so the
    grouping is independent of both predictor and response. This is the
    arbiter for any AGTR1-across-clusters claim.
  agtr1_count_by_program{,_posthoc}.tsv  grouping = state_program.
    state_program is a marker-panel argmax that INCLUDES the BM panel, so it
    is not independent of the matrix score and must never arbitrate a
    BM-vs-AGTR1 test. AGTR1 is in no panel, so it is admissible for the
    AGTR1-across-programs question only (figureS_acta2_control panel A).
Both carry `converged` and `max_grad`; an unconverged row is not a result.
