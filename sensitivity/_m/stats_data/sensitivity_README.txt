Sensitivity summary:
- Disease effects re-estimated with +smoking and +BMI covariates (see covariate_robustness_emmeans.tsv).
- Smoking x disease availability in smoking_availability_by_disease.tsv.
- Smoking MAIN effect among donors with a smoking label in smoking_main_effect_healthy.tsv.
- Smoking-stratified disease effects in smoking_stratified_injury.tsv.
- LIMITATION (smoking confound): HLCA records smoking_status ONLY for Healthy donors; it is
  missing for every diseased donor. A smoking-STRATIFIED disease contrast is therefore inestimable
  (estimable strata: none). The smoking signal is instead summarized as a main effect among
  donors that carry a smoking label (smoking_main_effect_healthy.tsv); smoking_stratified_injury.tsv
  is expected to be empty under the current metadata.
- Leave-one-STUDY-out stability of the Fibrotic_ILD effect in leave_one_study_out.tsv
  (column `dropped_study`). This is the cohort-robustness arm and is PRIMARY.
- Leave-one-DATASET-out in leave_one_dataset_out.tsv (column `dropped_dataset`).
  Batch robustness only: 5 studies span >1 dataset, so a dataset drop leaves the
  rest of that study in the fit and cannot remove a cohort. Changed 2026-09-07 (P2-8);
  before that date the dataset-level arm was written under the study-level filename.
- LIMITATION: HLCA lacks medication metadata; ARB/ACEi use cannot be adjusted for here.
