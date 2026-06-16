## Confounder sensitivity analyses for the headline donor-level associations.
##
## Re-fits disease effects on injury-state fraction, niche index, injury-stromal
## score, and AGTR1+ fraction with progressively added covariates:
##   base:        resp ~ disease_group + age + sex
##   + smoking:   ... + smoking_status
##   + comorbid:  ... + smoking_status + BMI
## plus smoking-stratified emmeans and leave-one-study-out (LOSO) robustness.
##
## LIMITATION (documented): the HLCA has NO medication metadata, so medication
## (e.g., ACE inhibitor / ARB use) sensitivity cannot be tested here; it is a
## limitation to state in the manuscript and a question for future cohorts.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(emmeans) })

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc) ~ "Healthy",
        lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
        TRUE ~ "Other")
}
collapse_smoking <- function(s) {
    s <- as.character(s)
    dplyr::case_when(s %in% c("never") ~ "never", s %in% c("former") ~ "former",
                     s %in% c("active") ~ "active", TRUE ~ "other/unknown")
}
write_tsv_safe <- function(x, f, rn = FALSE) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = rn)
}

NICHE  <- "../../niche_index/_m/niche_index_per_donor.tsv.gz"
STATES <- "../../pericyte_states/_m/pericytes_states_metadata.tsv.gz"
outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)

donor <- data.table::fread(NICHE) |>
    mutate(disease_group = relevel(factor(map_disease_group(lung_condition)), "Healthy"),
           smoking = factor(collapse_smoking(smoking_status)),
           age = suppressWarnings(as.numeric(age)), sex = factor(sex),
           BMI = suppressWarnings(as.numeric(BMI)))

## donor -> dataset map for LOSO
dataset_map <- data.table::fread(STATES) |>
    group_by(donor_id) |> summarise(dataset = dplyr::first(dataset), .groups = "drop")
donor <- left_join(donor, dataset_map, by = "donor_id")

RESPONSES <- c("injury_frac", "niche_index", "injury_stromal_score", "AGTR1_pos_frac")

## ---- (1) covariate robustness ------------------------------------------
covariate_robustness <- function(resp) {
    d <- donor |> tidyr::drop_na(all_of(resp), age, sex) |>
        mutate(disease_group = droplevels(disease_group))
    models <- list(
        base      = reformulate(c("disease_group", "age", "sex"), resp),
        smoking   = reformulate(c("disease_group", "age", "sex", "smoking"), resp),
        comorbid  = reformulate(c("disease_group", "age", "sex", "smoking", "BMI"), resp))
    rows <- lapply(names(models), function(m) {
        dm <- if (m == "comorbid") tidyr::drop_na(d, BMI) else d
        fit <- lm(models[[m]], data = dm)
        e <- as.data.frame(emmeans(fit, ~ disease_group))
        e$model <- m; e$response <- resp; e$n <- nrow(dm); e
    })
    bind_rows(rows)
}
cov_res <- bind_rows(lapply(RESPONSES, covariate_robustness))
write_tsv_safe(cov_res, file.path(outdir, "covariate_robustness_emmeans.tsv"))

## ---- (2) smoking analyses -----------------------------------------------
## NOTE: in the HLCA donor metadata, smoking_status is recorded ONLY for Healthy
## (incl. tumor-adjacent) donors and is missing for every diseased donor (COPD,
## IPF, ILD, etc.). A smoking-STRATIFIED disease effect is therefore inestimable
## -- each smoking stratum contains a single disease group (Healthy) -- so that
## model is dropped by the nlevels(disease_group) < 2 guard. We (a) record that
## confound explicitly, and (b) fit the smoking MAIN effect among the donors that
## do carry a smoking label (effectively the Healthy donors), which is what the
## data can actually support.
strat <- donor |> filter(smoking %in% c("never", "former", "active")) |>
    tidyr::drop_na(injury_stromal_score, age, sex)

## (2a) smoking x disease availability table (drives the confound note)
avail <- donor |> mutate(has_smk = smoking %in% c("never", "former", "active")) |>
    group_by(disease_group) |>
    summarise(n_donors = dplyr::n(), n_with_smoking = sum(has_smk),
              smoking_levels = paste(sort(unique(smoking[has_smk])), collapse = ","),
              .groups = "drop")
write_tsv_safe(avail, file.path(outdir, "smoking_availability_by_disease.tsv"))

## (2b) attempt the smoking-stratified disease effect (kept for completeness;
##      will be empty under the current metadata, documented in the note below)
strat_rows <- lapply(c("never", "former", "active"), function(sm) {
    sub <- strat |> filter(smoking == sm) |> mutate(disease_group = droplevels(disease_group))
    if (nlevels(sub$disease_group) < 2 || nrow(sub) < 8) return(NULL)
    fit <- lm(injury_stromal_score ~ disease_group + age, data = sub)
    e <- as.data.frame(emmeans(fit, ~ disease_group)); e$smoking <- sm; e$n <- nrow(sub); e
})
strat_res <- bind_rows(strat_rows)
estimable_strata <- c("never", "former", "active")[sapply(c("never","former","active"), function(sm) {
    sub <- strat |> filter(smoking == sm) |> mutate(disease_group = droplevels(disease_group))
    nlevels(sub$disease_group) >= 2 && nrow(sub) >= 8 })]
write_tsv_safe(strat_res, file.path(outdir, "smoking_stratified_injury.tsv"))

## (2c) smoking MAIN effect among donors with a smoking label (Healthy-dominated):
##      resp ~ smoking (+ age, + sex when >1 level). Answers "does smoking alone
##      shift the injury/niche phenotype?" using the data that exist.
smk_main <- function(resp) {
    d <- donor |> filter(smoking %in% c("never", "former", "active")) |>
        mutate(smoking = factor(smoking, levels = c("never", "former", "active"))) |>
        tidyr::drop_na(all_of(resp), age) |> droplevels()
    if (nlevels(d$smoking) < 2 || nrow(d) < 8) return(NULL)
    terms <- c("smoking", "age"); if (nlevels(factor(d$sex)) > 1) terms <- c(terms, "sex")
    fit <- lm(reformulate(terms, resp), data = d)
    e <- as.data.frame(emmeans(fit, ~ smoking)); e$response <- resp; e$n <- nrow(d)
    e
}
smk_main_res <- bind_rows(lapply(RESPONSES, smk_main))
write_tsv_safe(smk_main_res, file.path(outdir, "smoking_main_effect_healthy.tsv"))

## ---- (3) leave-one-study-out (LOSO) -------------------------------------
loso <- function(resp = "injury_stromal_score") {
    d <- donor |> tidyr::drop_na(all_of(resp), age, sex, dataset) |>
        mutate(disease_group = droplevels(disease_group))
    if (!"Fibrotic_ILD" %in% levels(d$disease_group)) return(NULL)
    out <- lapply(unique(d$dataset), function(ds) {
        sub <- d |> filter(dataset != ds) |> mutate(disease_group = droplevels(disease_group))
        if (!all(c("Healthy", "Fibrotic_ILD") %in% levels(sub$disease_group))) return(NULL)
        fit <- lm(reformulate(c("disease_group", "age", "sex"), resp), data = sub)
        ct <- summary(fit)$coefficients
        term <- grep("Fibrotic_ILD", rownames(ct), value = TRUE)[1]
        data.frame(dropped_dataset = ds, response = resp, n = nrow(sub),
                   estimate = ct[term, 1], se = ct[term, 2], p = ct[term, 4])
    })
    bind_rows(out)
}
write_tsv_safe(bind_rows(lapply(RESPONSES, loso)), file.path(outdir, "leave_one_study_out.tsv"))

writeLines(c(
    "Sensitivity summary:",
    "- Disease effects re-estimated with +smoking and +BMI covariates (see covariate_robustness_emmeans.tsv).",
    "- Smoking x disease availability in smoking_availability_by_disease.tsv.",
    "- Smoking MAIN effect among donors with a smoking label in smoking_main_effect_healthy.tsv.",
    "- Smoking-stratified disease effects in smoking_stratified_injury.tsv.",
    sprintf("- LIMITATION (smoking confound): HLCA records smoking_status ONLY for Healthy donors; it is"),
    "  missing for every diseased donor. A smoking-STRATIFIED disease contrast is therefore inestimable",
    sprintf("  (estimable strata: %s). The smoking signal is instead summarized as a main effect among",
            if (length(estimable_strata)) paste(estimable_strata, collapse = ",") else "none"),
    "  donors that carry a smoking label (smoking_main_effect_healthy.tsv); smoking_stratified_injury.tsv",
    "  is expected to be empty under the current metadata.",
    "- Leave-one-study-out stability of the Fibrotic_ILD effect in leave_one_study_out.tsv.",
    "- LIMITATION: HLCA lacks medication metadata; ARB/ACEi use cannot be adjusted for here."),
    file.path(outdir, "sensitivity_README.txt"))

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
