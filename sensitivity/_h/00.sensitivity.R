## Confounder sensitivity analyses for the headline donor-level associations.
##
## Re-fits disease effects on injury-state fraction, niche index, injury-stromal
## score, and AGTR1+ fraction with progressively added covariates:
##   base:        resp ~ disease_group + sex + (1 | study)      <- the PRIMARY model
##   + age:       ... + age                                     (age-complete donors only)
##   + smoking:   ... + smoking                                 (all donors; label is a level)
##   + comorbid:  ... + smoking + BMI
## plus smoking-stratified emmeans and leave-one-dataset-out (LOSO) robustness.
##
## ** REVISED 2026-09-07 (defect P1-10). ** The base model here used to be
## `lm(~ disease_group + age + sex)` -- no study term, and `drop_na(age)` applied
## before every fit. That mirrored `niche_index/_h/01`, which has now moved to
## `lmer(~ disease_group + sex + (1 | study))` for the reasons written into that
## script's header (age missingness is study-structured and deletes 42 of 89
## donors, three-quarters of them fibrotic). A covariate-robustness module whose
## OWN base model differs from the primary it is testing cannot do its job, so the
## base is realigned and `+ age` becomes one of the arms being tested -- which is
## what a covariate-sensitivity analysis should have been doing with it anyway.
## Each row carries `n` and `restriction` so a shrinking sample is visible rather
## than inferred.
##
## LOSO is over `dataset`, not `study`, and the column and axis say so (P1-11).
##
## LIMITATION (documented): the HLCA has NO medication metadata, so medication
## (e.g., ACE inhibitor / ARB use) sensitivity cannot be tested here; it is a
## limitation to state in the manuscript and a question for future cohorts.

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(emmeans); library(lmerTest)
})
emm_options(lmerTest.limit = 30000, pbkrtest.limit = 30000)

## Fit an lmer when the formula carries a random term, an lm otherwise.
fit_model <- function(form, data) {
    if (any(grepl("\\|", labels(terms(form)))))
        suppressMessages(lmerTest::lmer(form, data = data)) else lm(form, data = data)
}

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

## `dataset` and `study` now arrive on the niche table itself (added to
## 00.niche_index.py 2026-09-07). Joining a second copy here produced
## `dataset.x`/`dataset.y` and a hard failure, so the map is only built when the
## upstream table predates that change.
if (!all(c("dataset", "study") %in% names(donor))) {
    dataset_map <- data.table::fread(STATES) |>
        group_by(donor_id) |>
        summarise(dataset = dplyr::first(dataset), study = dplyr::first(study),
                  .groups = "drop")
    donor <- left_join(donor, dataset_map, by = "donor_id")
}
stopifnot(all(c("dataset", "study") %in% names(donor)))
donor <- donor |> mutate(study = factor(study), dataset = factor(dataset))
cat("donors:", nrow(donor), "| studies:", nlevels(droplevels(donor$study)),
    "| datasets:", nlevels(droplevels(donor$dataset)), "\n")
print(table(donor$disease_group))

RESPONSES <- c("injury_frac", "niche_index", "injury_stromal_score", "AGTR1_pos_frac")

## ---- (1) covariate robustness ------------------------------------------
covariate_robustness <- function(resp) {
    d <- donor |> tidyr::drop_na(all_of(resp), sex) |>
        mutate(disease_group = droplevels(disease_group))
    n_base <- nrow(d)
    ## Each arm names the rows it needs, so a covariate that is really a cohort
    ## filter shows up as a drop in `n` rather than as a quiet change in estimate.
    arms <- list(
        base     = list(terms = c("disease_group", "sex", "(1 | study)"),
                        need = character(0), note = "PRIMARY -- all donors"),
        age      = list(terms = c("disease_group", "sex", "age", "(1 | study)"),
                        need = "age",
                        note = "RESTRICTED to age-reporting cohorts -- not an age adjustment"),
        smoking  = list(terms = c("disease_group", "sex", "smoking", "(1 | study)"),
                        need = character(0),
                        note = "smoking label is a factor level; diseased donors are all other/unknown"),
        comorbid = list(terms = c("disease_group", "sex", "smoking", "BMI", "(1 | study)"),
                        need = "BMI", note = "RESTRICTED to BMI-reporting donors"))
    rows <- lapply(names(arms), function(m) {
        a  <- arms[[m]]
        dm <- if (length(a$need)) tidyr::drop_na(d, all_of(a$need)) else d
        dm <- mutate(dm, disease_group = droplevels(disease_group))
        if (nlevels(dm$disease_group) < 2 || dplyr::n_distinct(dm$study) < 2) {
            cat("SKIP", resp, m, "-- groups:", nlevels(dm$disease_group),
                "studies:", dplyr::n_distinct(dm$study), "\n")
            return(NULL)
        }
        fit <- try(fit_model(reformulate(a$terms, resp), dm), silent = TRUE)
        if (inherits(fit, "try-error")) { cat("FAIL", resp, m, "\n"); return(NULL) }
        e <- as.data.frame(emmeans(fit, ~ disease_group))
        e$model <- m; e$response <- resp; e$n <- nrow(dm)
        e$n_dropped_vs_base <- n_base - nrow(dm); e$restriction <- a$note
        e$n_donors_group <- paste(names(table(dm$disease_group)),
                                  as.integer(table(dm$disease_group)),
                                  sep = "=", collapse = ";")
        e
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
    tidyr::drop_na(injury_stromal_score, sex)

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
    fit <- lm(injury_stromal_score ~ disease_group, data = sub)
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
        tidyr::drop_na(all_of(resp)) |> droplevels()
    if (nlevels(d$smoking) < 2 || nrow(d) < 8) return(NULL)
    ## Age is not required here either: this arm is already restricted to the
    ## smoking-labelled donors, and stacking a second missingness filter on top
    ## makes the reported n uninterpretable.
    terms <- "smoking"; if (nlevels(factor(d$sex)) > 1) terms <- c(terms, "sex")
    if (dplyr::n_distinct(d$study) > 1) terms <- c(terms, "(1 | study)")
    fit <- try(fit_model(reformulate(terms, resp), d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- as.data.frame(emmeans(fit, ~ smoking)); e$response <- resp; e$n <- nrow(d)
    e$model <- paste(terms, collapse = " + ")
    e
}
smk_main_res <- bind_rows(lapply(RESPONSES, smk_main))
write_tsv_safe(smk_main_res, file.path(outdir, "smoking_main_effect_healthy.tsv"))

## ---- (3) leave-one-study-out (LOSO) -------------------------------------
## Refits the PRIMARY model with one dataset removed. `age` is deliberately NOT a
## drop_na target here: requiring it would silently make every LOSO refit a
## different, age-restricted analysis than the effect being tested.
loso <- function(resp = "injury_stromal_score") {
    d <- donor |> tidyr::drop_na(all_of(resp), sex, dataset, study) |>
        mutate(disease_group = droplevels(disease_group))
    if (!"Fibrotic_ILD" %in% levels(d$disease_group)) return(NULL)
    out <- lapply(levels(droplevels(d$dataset)), function(ds) {
        sub <- d |> filter(dataset != ds) |>
            mutate(disease_group = droplevels(disease_group), study = droplevels(study))
        if (!all(c("Healthy", "Fibrotic_ILD") %in% levels(sub$disease_group))) return(NULL)
        if (dplyr::n_distinct(sub$study) < 2) return(NULL)
        fit <- try(suppressMessages(lmerTest::lmer(
            reformulate(c("disease_group", "sex", "(1 | study)"), resp), data = sub)),
            silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        ct   <- summary(fit)$coefficients
        term <- grep("Fibrotic_ILD", rownames(ct), value = TRUE)[1]
        if (is.na(term)) return(NULL)
        data.frame(dropped_dataset = ds, response = resp, n = nrow(sub),
                   n_studies = dplyr::n_distinct(sub$study),
                   estimate = ct[term, "Estimate"], se = ct[term, "Std. Error"],
                   p = ct[term, ncol(ct)],
                   singular = lme4::isSingular(fit),
                   model = "lmer(~ disease_group + sex + (1 | study))")
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
