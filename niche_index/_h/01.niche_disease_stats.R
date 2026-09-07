## Donor-level disease association of the pericyte-endothelial niche index.
##
## Hypothesis (mirrors the in-vivo mouse pericyte loss + losartan rescue):
## fibrotic/ILD lungs have LOWER niche-stability and HIGHER injury-stromal
## scores than healthy lungs.
##
## ** REVISED 2026-09-07 (defect P1-10). The primary model changed twice over. **
##
## (1) `+ age` is dropped from the primary. It was never a covariate adjustment
##     here: age is missing in a *study*-structured way, so `drop_na(age)` deleted
##     42 of 89 donors and fitted 47 -- the same mechanism established for
##     `pericyte_states/_h/01.state_stats.R` (P1-2) and already handled correctly
##     by `disease_association/_h/03.disease_forest.R`. This module was the sixth
##     instance of that trap in the repository and the last one still unguarded.
##     The age-adjusted fit is still produced, on the age-complete subset, with the
##     `_ageadj` suffix and an `arm` column. It is a RESTRICTION TO AGE-REPORTING
##     COHORTS, not an age adjustment, and must be labelled that way.
##
## (2) `(1 | study)` is added. The previous model was a plain `lm` with no study
##     term of any kind, in the one compartment this repository has repeatedly
##     documented as carrying between-study confounding (see `sensitivity/`, and
##     the leave-one-study-out result that the injury effect loses significance
##     whenever a fibrosis-heavy cohort is dropped). Dropping age WITHOUT adding a
##     study term is strictly worse than doing neither -- restoring the deleted
##     donors restores the study imbalance they were hiding.
##
## HC3-robust SEs are retained only on the age-free fixed-effect `lm`, as a
## labelled comparison; they are not defined for the mixed fit and are not the
## primary readout. See P2-12 on their reliability at these group sizes.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(emmeans)
    library(lmerTest)
})
emm_options(lmerTest.limit = 30000, pbkrtest.limit = 30000)

map_disease_group <- function(lung_condition) {
    lc <- as.character(lung_condition)
    dplyr::case_when(
        grepl("^Healthy", lc)                                      ~ "Healthy",
        lc %in% c("COPD")                                          ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis",
              lc, ignore.case = TRUE)                              ~ "Fibrotic_ILD",
        TRUE                                                       ~ "Other"
    )
}

save_ggplots <- function(fn, p, w, h)
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), plot = p, width = w, height = h)

write_tsv_safe <- function(x, file, row_names = FALSE) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = row_names, col.names = TRUE)
}

## `--suffix` selects which donor table (and hence which donor cell-count
## threshold) to model: "" is the primary >=10 run, "_mincells20" the sensitivity
## run. Outputs inherit the same suffix so the two never overwrite each other.
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
SFX <- parse_arg("--suffix", "")

infile <- paste0("niche_index_per_donor", SFX, ".tsv.gz")
if (!file.exists(infile)) stop("missing ", infile, " -- run 00.niche_index.py first")
df <- data.table::fread(infile) |>
    mutate(disease_group = relevel(factor(map_disease_group(lung_condition)), "Healthy"),
           sex = factor(sex), age = suppressWarnings(as.numeric(age))) |>
    filter(!is.na(disease_group))

## `study` arrives from 00.niche_index.py. A donor maps to exactly one study, so
## this is a clean grouping factor. Guard rather than assume: an older per-donor
## table has no such column, and silently falling back to a plain `lm` is exactly
## the failure this revision exists to remove.
HAS_STUDY <- "study" %in% names(df) && dplyr::n_distinct(df$study) > 1
if (!HAS_STUDY)
    stop("no usable `study` column in ", infile, " -- re-run 00.niche_index.py. ",
         "Fitting these models without a study random intercept is defect P1-10; ",
         "do not proceed by dropping the guard.")
cat("study levels:", dplyr::n_distinct(df$study), "\n")

MIN_CELLS <- if ("min_cells" %in% names(df)) unique(df$min_cells)[1] else NA_integer_
cat("input:", infile, "| min_cells:", MIN_CELLS, "| donors:", nrow(df), "\n")
print(table(df$disease_group))

outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)

## The two arms. `_ageadj` is a restriction to age-reporting cohorts, and the
## `arm` column on every output says so, so that no downstream reader has to infer
## it from an `n_donors` that happens to be smaller.
COVARS_PRIMARY  <- c("disease_group", "sex", "(1 | study)")
COVARS_AGE_SENS <- c("disease_group", "sex", "age", "(1 | study)")

## SMALL-GROUP GUARD. Dropping `+ age` restored COPD to these models -- and the
## COPD arm is ONE donor. A singleton group is a fitted value with no residual
## information, and it lands on the omnibus F and on three of the six contrasts.
## The pattern here is the one `pericyte_states/_h/01.state_stats.R` settled on:
## fit the full model, then refit with the small groups removed and ship BOTH
## p-values side by side, so a reader never has to reconstruct which groups were
## load-bearing. `estimable` is FALSE when the guard could not be evaluated.
MIN_GROUP_N <- 3L

refit_excl_small <- function(sub, response, covars) {
    tab   <- table(sub$disease_group)
    small <- names(tab)[tab < MIN_GROUP_N]
    keep  <- sub[!sub$disease_group %in% small, ]
    keep$disease_group <- droplevels(keep$disease_group)
    out <- list(small_groups = if (length(small)) paste(small, collapse = ";") else "",
                n_excl = nrow(sub) - nrow(keep), p_omnibus = NA_real_,
                posthoc = NULL, estimable = FALSE)
    if (!length(small)) {                       # nothing to exclude
        out$estimable <- TRUE; out$small_groups <- "none"; return(out)
    }
    if (nlevels(keep$disease_group) < 2 || dplyr::n_distinct(keep$study) < 2) return(out)
    f2 <- try(suppressMessages(lmerTest::lmer(reformulate(covars, response), data = keep)),
              silent = TRUE)
    if (inherits(f2, "try-error")) return(out)
    a2 <- as.data.frame(anova(f2))
    out$p_omnibus <- a2[["Pr(>F)"]][rownames(a2) == "disease_group"]
    out$posthoc   <- as.data.frame(pairs(emmeans(f2, ~ disease_group), adjust = "BH"))
    out$estimable <- TRUE
    out
}

## One run per (response, arm). Returns the fitted model invisibly.
run_one <- function(response, covars = COVARS_PRIMARY, arm = "", drop_age_rows = FALSE) {
    sub <- df |> tidyr::drop_na(all_of(response), sex)
    if (drop_age_rows) sub <- tidyr::drop_na(sub, age)
    sub <- mutate(sub, disease_group = droplevels(disease_group))

    ## Guard against a degenerate refit: if the age restriction has left a single
    ## study or a single disease group, an lmer either fails or returns a fit whose
    ## marginal means invite quotation (P2-9). Say so and skip, rather than write it.
    if (dplyr::n_distinct(sub$study) < 2 || nlevels(sub$disease_group) < 2) {
        cat("SKIP", response, "arm", if (nzchar(arm)) arm else "(primary)",
            "-- studies:", dplyr::n_distinct(sub$study),
            "groups:", nlevels(sub$disease_group), "\n")
        return(invisible(NULL))
    }

    form  <- reformulate(covars, response)
    mtag  <- paste0("lmer(", deparse(form[[3]]), ")")
    fit   <- suppressMessages(lmerTest::lmer(form, data = sub))
    emm   <- emmeans(fit, ~ disease_group)
    gd    <- refit_excl_small(sub, response, covars)
    vc    <- as.data.frame(lme4::VarCorr(fit))
    sd_st <- vc$sdcor[vc$grp == "study"][1]
    sing  <- lme4::isSingular(fit)

    ## n_donors, min_cells, arm and the study-variance diagnostics travel with
    ## every output: a table that cannot say which arm and which threshold produced
    ## an estimate is a table that will be quoted as the primary one.
    tag <- function(x) as.data.frame(x) |>
        mutate(n_donors = nrow(sub), min_cells = MIN_CELLS, arm = if (nzchar(arm)) arm else "primary",
               n_studies = dplyr::n_distinct(sub$study), study_sd = sd_st,
               singular = sing,
               n_donors_group = paste(names(table(sub$disease_group)),
                                      as.integer(table(sub$disease_group)),
                                      sep = "=", collapse = ";"),
               small_groups = gd$small_groups, estimable = gd$estimable,
               model = mtag)

    av <- as.data.frame(anova(fit))
    av$n_donors <- nrow(sub); av$min_cells <- MIN_CELLS
    av$arm <- if (nzchar(arm)) arm else "primary"
    av$n_studies <- dplyr::n_distinct(sub$study); av$study_sd <- sd_st
    av$singular <- sing
    ## The omnibus p with every group under MIN_GROUP_N donors removed. NA when the
    ## refit is not identified; `estimable` says which.
    av$p_excl_small_groups <- ifelse(rownames(av) == "disease_group", gd$p_omnibus, NA_real_)
    av$small_groups <- gd$small_groups; av$estimable <- gd$estimable
    av$n_donors_group <- paste(names(table(sub$disease_group)),
                               as.integer(table(sub$disease_group)),
                               sep = "=", collapse = ";")
    av$model <- mtag
    write_tsv_safe(av, file.path(outdir, paste0(response, "_anova", SFX, arm, ".tsv")), TRUE)
    write_tsv_safe(tag(emm), file.path(outdir, paste0(response, "_emmeans", SFX, arm, ".tsv")))
    ph <- tag(pairs(emm, adjust = "BH"))
    ## Carry the guarded p per contrast where the guarded refit still contains it.
    ph$p_excl_small_groups <- NA_real_
    if (!is.null(gd$posthoc)) {
        m <- match(as.character(ph$contrast), as.character(gd$posthoc$contrast))
        ph$p_excl_small_groups <- gd$posthoc$p.value[m]
    }
    write_tsv_safe(ph, file.path(outdir, paste0(response, "_posthoc", SFX, arm, ".tsv")))

    ## HC3 on the age-free FIXED-EFFECT fit only, kept as a labelled comparison so
    ## the change of estimator is visible. It is NOT the primary readout: it has no
    ## study term, which is the whole point of the revision above.
    lfit <- lm(reformulate(setdiff(covars, "(1 | study)"), response), data = sub)
    rb <- as.data.frame(unclass(lmtest::coeftest(lfit, vcov = sandwich::vcovHC(lfit, type = "HC3"))))
    rb$term <- rownames(rb)
    rb <- tag(rb); rb$model <- "lm + HC3 -- UNGUARDED, comparison only"
    write_tsv_safe(rb, file.path(outdir, paste0(response, "_robust_coefs", SFX, arm, ".tsv")))

    if (!nzchar(arm)) {
        p <- ggboxplot(sub, x = "disease_group", y = response, add = "jitter",
                       fill = "disease_group", palette = "jco",
                       add.params = list(alpha = 0.5, size = 1.2),
                       xlab = "", ylab = response, legend = "none",
                       ggtheme = theme_pubr(base_size = 13)) +
            rotate_x_text(30) +
            stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                         fill = "white", color = "black")
        save_ggplots(file.path(outdir, paste0("box_", response, SFX)), p, 5, 5)
    }
    invisible(fit)
}

## PRIMARY responses (main narrative + main figures).
PRIMARY_RESP <- c("niche_stability_score", "injury_stromal_score", "niche_index")
## SENSITIVITY responses: injury / index composites that additionally fold in the
## AGTR1+ pericyte fraction. Reported in the supplement only -- the AGTR1+ fraction
## is kept out of the primary composite to avoid circularity with the focal receptor.
SENS_RESP    <- c("injury_stromal_score_sens_agtr1", "niche_index_sens_agtr1")

for (resp in c(PRIMARY_RESP, SENS_RESP)) {
    if (!resp %in% names(df)) next
    run_one(resp, COVARS_PRIMARY,  "",         drop_age_rows = FALSE)
    run_one(resp, COVARS_AGE_SENS, "_ageadj",  drop_age_rows = TRUE)
}

cat("\nNOTE: A lower niche-index / higher injury-stromal score in fibrotic/ILD\n",
    "lungs is the human transcriptomic counterpart of the in-vivo mouse\n",
    "pericyte-loss phenotype rescued by losartan (AT1 blockade).\n")

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
