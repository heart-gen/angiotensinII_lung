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
        grp_n <- table(dm$disease_group)
        e$n_donors_group <- paste(names(grp_n), as.integer(grp_n),
                                  sep = "=", collapse = ";")

        ## ---- NON-ESTIMABLE GUARD (P2-9, added 2026-09-07) -----------------
        ## Announcing a restriction in a `restriction` column was not enough: the
        ## `+BMI` arm still printed an "Other" marginal mean of 2.023 from ONE
        ## donor, and a number on the page is what gets quoted. Values that
        ## cannot support quotation are now blanked to NA and the reason is
        ## carried in the row, so a reader cannot reach the estimate at all.
        e$n_group <- as.integer(grp_n[as.character(e$disease_group)])
        reason <- rep(NA_character_, nrow(e))

        ## (a) row-level: a marginal mean from < MIN_GROUP donors is not a result.
        reason[e$n_group < MIN_GROUP_DONORS] <-
            sprintf("group has %d donor(s), floor is %d",
                    e$n_group[e$n_group < MIN_GROUP_DONORS], MIN_GROUP_DONORS)

        ## (b) arm-level: this block exists to re-estimate the Healthy-vs-Fibrotic
        ## contrast under extra covariates. An arm that cannot see both arms of
        ## that contrast is not a robustness check of it, whatever its other
        ## groups look like, so the whole arm is withheld.
        ok_contrast <- all(c("Healthy", "Fibrotic_ILD") %in% names(grp_n)) &&
            all(grp_n[c("Healthy", "Fibrotic_ILD")] >= MIN_GROUP_DONORS)
        if (!ok_contrast)
            reason <- sprintf("arm cannot estimate Healthy vs Fibrotic_ILD (%s)",
                              e$n_donors_group)

        ## (c) identification: if only ONE level of an added covariate contains
        ## more than one disease group, the disease contrast is identified inside
        ## that level alone rather than adjusted for the covariate. Computed, not
        ## assumed, so it keeps reporting the truth if the metadata changes.
        cov_extra <- setdiff(a$terms, c("disease_group", "sex", "(1 | study)"))
        ## dm may be a data.table, where dm[chr] is a join rather than a column
        ## subset -- select one column at a time.
        cov_extra <- cov_extra[cov_extra %in% names(dm)]
        cov_extra <- cov_extra[!vapply(cov_extra, function(cv) is.numeric(dm[[cv]]), logical(1))]
        ident <- vapply(cov_extra, function(cv) {
            tb <- table(dm[[cv]], dm$disease_group) > 0
            sum(rowSums(tb) > 1) <= 1
        }, logical(1))
        e$identified_within_one_level <- if (length(ident)) paste(
            names(ident)[ident], collapse = ",") else NA_character_
        e$identified_within_one_level[!nzchar(e$identified_within_one_level)] <- NA

        e$non_estimable <- !is.na(reason)
        e$non_estimable_reason <- reason
        for (cl in intersect(c("emmean", "SE", "df", "lower.CL", "upper.CL",
                               "asymp.LCL", "asymp.UCL"), names(e)))
            e[[cl]][e$non_estimable] <- NA_real_
        if (any(e$non_estimable))
            cat(sprintf("WITHHELD %s/%s: %d of %d marginal means (%s)\n", resp, m,
                        sum(e$non_estimable), nrow(e), reason[which(e$non_estimable)[1]]))
        e
    })
    bind_rows(rows)
}
## Floor for quoting a marginal mean. Three donors is not a defensible estimate
## either, but it is the point below which the mean IS the donor.
MIN_GROUP_DONORS <- 3L

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
    emm <- emmeans(fit, ~ smoking)
    e <- as.data.frame(emm); e$response <- resp; e$n <- nrow(d)
    e$model <- paste(terms, collapse = " + ")
    ## PAIRWISE CONTRASTS (P2-10, added 2026-09-07). The block used to write the
    ## marginal means and stop, which left the module's ONLY estimable smoking
    ## question untested -- and the means are not flat: active smokers separate
    ## from never/former on all four endpoints. Reporting means whose separation
    ## a reader can see but cannot test is the defect. `n_smoking` records the
    ## donors per level, because these strata are small and the contrast is only
    ## as good as the smaller arm.
    ct <- as.data.frame(pairs(emm, adjust = "none"))
    ct$response <- resp; ct$n <- nrow(d)
    ct$model <- paste(terms, collapse = " + ")
    ct$n_smoking <- paste(sprintf("%s=%d", levels(d$smoking),
                                  as.integer(table(d$smoking))), collapse = ";")
    list(emmeans = e, contrasts = ct)
}
smk_main_all  <- lapply(RESPONSES, smk_main)
smk_main_res  <- bind_rows(lapply(smk_main_all, `[[`, "emmeans"))
smk_main_ct   <- bind_rows(lapply(smk_main_all, `[[`, "contrasts"))
write_tsv_safe(smk_main_res, file.path(outdir, "smoking_main_effect_healthy.tsv"))
## BH across every contrast this block writes (3 pairs x 4 responses), declared
## as one family so the adjustment is not silently per-response.
if (nrow(smk_main_ct)) {
    smk_main_ct$p_BH <- p.adjust(smk_main_ct$p.value, method = "BH")
    smk_main_ct$bh_family <- sprintf("smoking main effect: %d contrasts x %d responses",
                                     nrow(smk_main_ct) / max(1L, dplyr::n_distinct(smk_main_ct$response)),
                                     dplyr::n_distinct(smk_main_ct$response))
}
write_tsv_safe(smk_main_ct, file.path(outdir, "smoking_main_effect_contrasts.tsv"))
cat("\n== smoking main effect: pairwise contrasts (P2-10) ==\n")
print(smk_main_ct)

## ---- (3) leave-one-out refits ------------------------------------------
## GROUPING FIXED 2026-09-07 (P2-8). This block looped over `dataset` while
## writing `leave_one_study_out.tsv` -- a filename that contradicted the column
## inside it. The distinction is not cosmetic: 18 studies are split across 23
## datasets here, and five studies (Sun_2020, Thienpont_2018, Lafyatis_Rojas_2019,
## Meyer_2021, Regev_2021) contribute more than one. Dropping one dataset leaves
## the rest of its study in the fit, so for those five the analysis could not
## test the objection it exists to answer -- "is this effect carried by one
## cohort?" -- no matter how many refits it ran.
##
## Both groupings are now emitted, because they answer different questions:
##   study-level   -> leave_one_study_out.tsv   (cohort robustness; the PRIMARY)
##   dataset-level -> leave_one_dataset_out.tsv (batch robustness; the old arm)
##
## `age` is deliberately NOT a drop_na target here: requiring it would silently
## make every refit a different, age-restricted analysis than the effect tested.
loso <- function(resp = "injury_stromal_score", by = "study") {
    d <- donor |> tidyr::drop_na(all_of(resp), sex, dataset, study) |>
        mutate(disease_group = droplevels(disease_group))
    if (!"Fibrotic_ILD" %in% levels(d$disease_group)) return(NULL)
    out <- lapply(levels(droplevels(d[[by]])), function(g) {
        sub <- d |> filter(.data[[by]] != g) |>
            mutate(disease_group = droplevels(disease_group), study = droplevels(study))
        if (!all(c("Healthy", "Fibrotic_ILD") %in% levels(sub$disease_group))) return(NULL)
        ## A study random effect needs >= 2 remaining studies. Dropping a whole
        ## study can take the fit below that where dropping a dataset would not,
        ## so this guard bites harder on the study-level arm -- by design.
        if (dplyr::n_distinct(sub$study) < 2) return(NULL)
        fit <- try(suppressMessages(lmerTest::lmer(
            reformulate(c("disease_group", "sex", "(1 | study)"), resp), data = sub)),
            silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        ct   <- summary(fit)$coefficients
        term <- grep("Fibrotic_ILD", rownames(ct), value = TRUE)[1]
        if (is.na(term)) return(NULL)
        r <- data.frame(dropped = g, dropped_level = by, response = resp,
                        n = nrow(sub), n_dropped = sum(d[[by]] == g, na.rm = TRUE),
                        n_studies = dplyr::n_distinct(sub$study),
                        n_datasets = dplyr::n_distinct(sub$dataset),
                        estimate = ct[term, "Estimate"], se = ct[term, "Std. Error"],
                        p = ct[term, ncol(ct)],
                        singular = lme4::isSingular(fit),
                        model = "lmer(~ disease_group + sex + (1 | study))")
        ## Keep the historical column name on each arm so downstream readers bind
        ## to something that says what was actually dropped.
        names(r)[names(r) == "dropped"] <- paste0("dropped_", by)
        r
    })
    bind_rows(out)
}
loso_study   <- bind_rows(lapply(RESPONSES, loso, by = "study"))
loso_dataset <- bind_rows(lapply(RESPONSES, loso, by = "dataset"))
write_tsv_safe(loso_study,   file.path(outdir, "leave_one_study_out.tsv"))
write_tsv_safe(loso_dataset, file.path(outdir, "leave_one_dataset_out.tsv"))
cat(sprintf("\nLOSO refits (P2-8): %d study-level over %d studies, %d dataset-level over %d datasets\n",
            nrow(loso_study), dplyr::n_distinct(donor$study),
            nrow(loso_dataset), dplyr::n_distinct(donor$dataset)))

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
    "- Leave-one-STUDY-out stability of the Fibrotic_ILD effect in leave_one_study_out.tsv",
    "  (column `dropped_study`). This is the cohort-robustness arm and is PRIMARY.",
    "- Leave-one-DATASET-out in leave_one_dataset_out.tsv (column `dropped_dataset`).",
    "  Batch robustness only: 5 studies span >1 dataset, so a dataset drop leaves the",
    "  rest of that study in the fit and cannot remove a cohort. Changed 2026-09-07 (P2-8);",
    "  before that date the dataset-level arm was written under the study-level filename.",
    "- LIMITATION: HLCA lacks medication metadata; ARB/ACEi use cannot be adjusted for here."),
    file.path(outdir, "sensitivity_README.txt"))

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
