## Donor-aware statistics on the de novo pericyte states (NVU pattern).
##
## States are STABLE Leiden clusters on the study-integrated embedding
## (`pericyte_state`), annotated to a dominant curated program (`state_program`)
## by 00.state_discovery.py. The unit of replication is the donor throughout.
##
## ** REVISED 2026-09-02: the donor-level disease models DO now carry `(1 | study)`. **
## This file previously argued that study was "handled once, by the integration
## that the clustering runs on". That is true of the EMBEDDING but not of the
## donor-level composition fractions: integration harmonises where a cell lands,
## not how many pericytes of each state a given cohort's donors contribute, which
## still varies with tissue sampling, dissociation and disease definition.
## Measured, once the age filter stopped hiding it: without a study term
## basement_membrane composition gives p = 0.0015 and vascular_stabilizing
## p = 0.0018 against disease; with `(1 | study)` the same fits give p = 0.797 and
## p = 0.912. The apparent effect sat in the "Other" group, 14 of whose 22 donors
## are Regev_2021. The integration argument did not survive contact with the
## restored donors.
##
##   (A) AGTR1 across states / programs (donor x group mixed model) -- this is the
##       RAW-EXPRESSION lens ONLY. Its apparent vascular-stabilizing enrichment is a
##       dropout/transcript-capture artifact that REVERSES under scVI denoising, so it
##       is superseded by the three-lens analysis (03.agtr1_lenses.R). Kept as a
##       diagnostic; do NOT cite AGTR1_by_* as a state-marker claim (AGTR1 labels the
##       pericyte/mural compartment, not a discrete state).
##   (B) Composition vs disease: do specific stable clusters / programs expand
##       in fibrosis/ILD? (donor-level ANCOVA, BH across clusters).
##   (C) Injury-program fraction vs disease (programs grouped; headline contrast).
##
## NOTE: COPD is very sparse among pericytes (~41 cells) and few fibrotic donors
## clear the donor cell-count filter (>=10 primary, >=20 sensitivity); the powered
## contrast is Healthy vs Fibrotic/ILD and even that is donor-limited. The COPD AGTR1 signal is carried by the whole-stroma
## disease_association analysis.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(rlang)
    library(ggpubr)
    library(lme4)
    library(lmerTest)
    library(emmeans)
})
emm_options(lmerTest.limit = 30000, pbkrtest.limit = 30000)

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
## Donor cell-count thresholds. The FIRST is primary and writes the canonical
## unsuffixed filenames; the rest are sensitivity analyses and are suffixed
## `_mincells<N>`. >=10 matches disease_association/_h/03.disease_forest.R, which
## has always used 10 -- the modules previously disagreed (10 vs 20) with nothing
## in the outputs revealing it, so two incompatible donor denominators were being
## reported side by side.
MIN_CELLS <- as.integer(strsplit(parse_arg("--min-cells", "10,20"), ",")[[1]])

## Programs counted as "injury-associated" for the grouped endpoint (C).
##
## HISTORY -- read before editing. This vector used to be
##   c("inflammatory", "fibroblast_like", "activated_migratory")
## but after the 2026-07-21 relabel the only `state_program` levels that exist are
## vascular_stabilizing, basement_membrane and activated_migratory:
## `fibroblast_like` became `basement_membrane` (a STRUCTURAL program that is
## deliberately not injury -- see basement_membrane/), and `inflammatory` is not
## dominant for any stable cluster. The stale names matched nothing, so the
## endpoint silently narrowed to activated/migratory alone while the summaries
## kept quoting the old three-program means (0.488/0.381/0.317). Naming only
## levels that exist, plus the hard check in injury_fraction_by_disease(), is what
## stops that from recurring.
INJURY_PROGRAMS <- c("activated_migratory")

map_disease_group <- function(lung_condition) {
    lc <- as.character(lung_condition)
    dplyr::case_when(
        grepl("^Healthy", lc)                                   ~ "Healthy",
        lc %in% c("COPD")                                       ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|Systemic sclerosis",
              lc, ignore.case = TRUE)                           ~ "Fibrotic_ILD",
        TRUE                                                     ~ "Other"
    )
}

load_meta <- function(path = "pericytes_states_metadata.tsv.gz") {
    data.table::fread(path) |>
        mutate(
            disease_group  = factor(map_disease_group(lung_condition)),
            age            = suppressWarnings(as.numeric(age_or_mean_of_age_range)),
            sex            = factor(sex),
            ethnicity      = factor(self_reported_ethnicity),
            donor_id       = factor(donor_id),
            pericyte_state = factor(pericyte_state),
            state_program  = factor(state_program)
        )
}

save_ggplots <- function(fn, p, w, h)
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), plot = p, width = w, height = h)

write_tsv_safe <- function(x, file, row_names = FALSE) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = row_names, col.names = TRUE)
}

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

## Post-hoc contrasts carrying BOTH a BH-adjusted p-value and an interval, so the
## composition results can be drawn as a forest plot. BH has no interval analogue,
## so the CI is the NOMINAL 95% interval (adjust = "none") while the p-value stays
## BH-adjusted across contrasts; interval and p-value therefore need not agree at
## the 0.05 boundary. Column order matches the pre-existing *_posthoc.tsv schema
## (contrast, estimate, SE, df, t.ratio, p.value) with lower/upper appended.
posthoc_with_ci <- function(emm) {
    pr  <- pairs(emm, adjust = "BH")
    est <- as.data.frame(pr)
    ci  <- as.data.frame(confint(pr, adjust = "none"))
    m   <- match(est$contrast, ci$contrast)
    est$lower.CL <- ci$lower.CL[m]
    est$upper.CL <- ci$upper.CL[m]
    est
}

## ----- (A) AGTR1 across states / programs (donor x group mixed model) -----
agtr1_by_group <- function(df, group, outdir, tag, min_cells = 5) {
    agg <- df |>
        group_by(donor_id, .data[[group]]) |>
        summarise(AGTR1_mean = mean(AGTR1_expr, na.rm = TRUE), n_cells = n(),
                  disease_group = first(disease_group), sex = first(sex),
                  age = mean(age, na.rm = TRUE), .groups = "drop") |>
        filter(n_cells >= min_cells) |>
        ## `age` deliberately not required here either -- see the rationale block
        ## below composition_by_disease(). Requiring it would fit this model on a
        ## different (47-donor, 5-studies-deleted) cohort than the composition and
        ## injury models, which is how the module came to report disease results
        ## on incompatible donor sets in the first place.
        tidyr::drop_na(AGTR1_mean, sex)
    agg[[group]]       <- droplevels(factor(agg[[group]]))
    agg$disease_group  <- droplevels(agg$disease_group)
    agg$sex            <- droplevels(agg$sex)
    if (nlevels(agg[[group]]) < 2) return(invisible(NULL))

    # Donor random intercept accounts for within-donor correlation across groups.
    form <- reformulate(c(group, "disease_group", "sex", "(1 | donor_id)"),
                        "AGTR1_mean")
    fit <- suppressMessages(lmerTest::lmer(form, data = agg))
    emm <- emmeans(fit, specs = group)
    write_tsv_safe(as.data.frame(anova(fit)),
                   file.path(outdir, paste0("AGTR1_by_", tag, "_anova.tsv")), TRUE)
    write_tsv_safe(as.data.frame(emm),
                   file.path(outdir, paste0("AGTR1_by_", tag, "_emmeans.tsv")))
    write_tsv_safe(as.data.frame(pairs(emm, adjust = "BH")),
                   file.path(outdir, paste0("AGTR1_by_", tag, "_posthoc.tsv")))

    p <- ggboxplot(agg, x = group, y = "AGTR1_mean", add = "jitter",
                   fill = group, palette = "npg",
                   add.params = list(alpha = 0.5, size = 1),
                   xlab = tag, ylab = "Mean AGTR1 (norm. expr.)",
                   legend = "none", ggtheme = theme_pubr(base_size = 13)) +
        rotate_x_text(35) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black")
    save_ggplots(file.path(outdir, paste0("boxplot_AGTR1_by_", tag)), p, 6, 5)
    invisible(agg)
}

## ----- (B) Composition vs disease (per stable cluster / program) ----------
## ---------------------------------------------------------------------------
## WHY THE PRIMARY MODELS NO LONGER CARRY `+ age`  (changed 2026-09-02)
##
## `+ age` was never a covariate adjustment here. Age is missing for 46 of the 93
## donors passing the cell filter, and the missingness is a STUDY property, not a
## donor property: 17 of 18 studies are all-or-nothing (only Banovich_Kropski_2020
## is partial, 8/15). So dropping incomplete rows deletes five whole studies.
##
## That deletion is not random with respect to the exposure. Retention by group:
##
##     Healthy       38 / 42   (90%)
##     Fibrotic/ILD   6 / 19   (32%)
##     Other          3 / 28   (11%)
##     COPD           0 / 1
##
## It removes entire dedicated fibrosis cohorts -- Kaminski_2020 (6 IPF + 1 COPD)
## and Sheppard_2020 (5 fibrotic) -- leaving the disease contrast resting on one
## or two remaining fibrosis studies. Adding age to control confounding therefore
## makes disease and study MORE collinear than they were: the cure introduces a
## worse confound than the one it treats.
##
## And age does not predict these outcomes anyway. Within the 47 age-complete
## donors, with study in the model, NO outcome shows an age effect at BH < 0.05
## (min BH = 0.225; only activated_migratory_score is even nominal, p = 0.025).
## `age` also spans 0-75 years here, crossing developmental stages, so a linear
## term is questionable regardless.
##
## Multiple imputation is deliberately NOT used: five studies have zero age
## observations, so there is no within-study information to borrow and imputation
## would extrapolate across studies -- and study is confounded with disease, so it
## would fabricate the very structure under test.
##
## `disease_association/_h/03.disease_forest.R` already dropped age for exactly
## this reason. Matching it here also ends the state of two modules reporting
## disease results on incompatible donor sets.
##
## The age-adjusted fit is still produced, as an explicit SENSITIVITY on the
## age-complete subset, written with the `_ageadj` suffix. It is labelled for what
## it is -- a restriction to age-reporting cohorts -- not as "age-adjusted".
##
## THE STUDY RANDOM EFFECT IS NOT OPTIONAL ONCE AGE IS DROPPED.
## Restoring the 42 donors puts five whole studies back, and the composition
## models had no study term -- so the restored between-study variance flowed
## straight into the disease coefficient and manufactured an effect. Measured:
## basement_membrane p = 0.0015 and vascular_stabilizing p = 0.0018 without a
## study term, but p = 0.797 and p = 0.912 with `(1 | study)`. The apparent
## signal sat in the "Other" group, of which 14 of 22 donors are Regev_2021.
## Dropping age and omitting `(1 | study)` is strictly worse than doing neither.
COVARS_PRIMARY   <- c("disease_group", "sex", "(1 | study)")
COVARS_AGE_SENS  <- c("disease_group", "sex", "age", "(1 | study)")

## Fit an lmer when the formula carries a random term, an lm otherwise, and take
## the omnibus disease test from whichever was fitted.
fit_model <- function(covars, response, data) {
    f <- reformulate(covars, response)
    if (any(grepl("\\|", covars)))
        suppressMessages(lmerTest::lmer(f, data = data))
    else lm(f, data = data)
}
disease_omnibus <- function(fit) {
    if (inherits(fit, "merMod")) {
        a <- as.data.frame(anova(fit))
        r <- a["disease_group", , drop = FALSE]
        data.frame(Df = r[["NumDF"]], F.value = r[["F value"]],
                   `Pr..F.` = r[["Pr(>F)"]], check.names = FALSE)
    } else {
        a <- as.data.frame(car::Anova(fit, type = 2))["disease_group", ]
        data.frame(Df = a[["Df"]], F.value = a[["F value"]],
                   `Pr..F.` = a[["Pr(>F)"]], check.names = FALSE)
    }
}

## n_donors must come from the FITTED object. lm() performs NA deletion after
## nrow() is taken, so counting rows reported 93 donors for models that fitted 47
## -- the exported df of 42 (= 47 - 5 parameters) proved it against itself.
fit_n <- function(fit) tryCatch(stats::nobs(fit), error = function(e) NA_integer_)

composition_by_disease <- function(df, group, outdir, tag, min_cells_per_donor = 10,
                                   sfx = "") {
    donor_tot <- df |> count(donor_id, name = "n_total") |>
        filter(n_total >= min_cells_per_donor)
    comp <- df |>
        semi_join(donor_tot, by = "donor_id") |>
        count(donor_id, .data[[group]], name = "n") |>
        tidyr::complete(donor_id, !!rlang::sym(group), fill = list(n = 0)) |>
        left_join(donor_tot, by = "donor_id") |>
        mutate(frac = n / n_total)
    donor_meta <- df |> group_by(donor_id) |>
        summarise(disease_group = first(disease_group), sex = first(sex),
                  dataset = first(dataset), study = first(study),
                  age = mean(age, na.rm = TRUE),
                  .groups = "drop")
    comp <- comp |> left_join(donor_meta, by = "donor_id") |>
        mutate(disease_group = relevel(factor(disease_group), ref = "Healthy"))

    ## Donor-level source data for the supplementary composition figure. The models
    ## below only emit marginal means, so without this the per-donor points cannot
    ## be redrawn downstream.
    ## donor_id is a factor carrying every donor in the dataset, so the complete()
    ## above re-introduces donors that the >=20-cell filter dropped, as n = 0 with
    ## n_total = NA. lm() silently ignores them; an exported table must not.
    ## `min_cells` is carried on every exported table so the primary and
    ## sensitivity fits can be stacked without relying on the filename.
    write_tsv_safe(
        comp |> filter(!is.na(n_total)) |>
            rename(level = !!rlang::sym(group)) |>
            select(donor_id, level, disease_group, dataset, sex, age,
                   n, n_total, frac) |>
            mutate(min_cells = min_cells_per_donor) |>
            arrange(level, disease_group, donor_id),
        file.path(outdir, paste0("composition_", tag, "_by_donor", sfx, ".tsv")))

    ## Two arms per level: the primary fit on every donor, and the age-restricted
    ## sensitivity. Emitted through one loop so they can never drift apart.
    run_arm <- function(covars, arm_sfx, drop_age_rows) {
        results <- list()
        for (g in levels(factor(comp[[group]]))) {
            sub <- comp |> filter(.data[[group]] == g)
            sub <- if (drop_age_rows) tidyr::drop_na(sub, age, sex) else
                tidyr::drop_na(sub, sex)
            sub <- sub |> mutate(disease_group = droplevels(disease_group))
            if (nlevels(sub$disease_group) < 2) next
            ## Skip degenerate levels (e.g. a single-program grouping makes frac==1
            ## everywhere -> zero residual variance -> Anova.lm errors).
            if (sd(sub$frac, na.rm = TRUE) < 1e-9) next
            fit <- fit_model(covars, "frac", sub)
            emm <- emmeans(fit, ~ disease_group)
            key <- sanitize(g)
            ## Per-group donor counts and a small-group flag, for the same reason
            ## they are on the injury endpoint: COPD is n = 1 here, and on its own
            ## it drives cluster 5 to BH = 0.0005 (p = 0.51 once it is excluded).
            ## A significant omnibus that rests on one donor must not reach a
            ## reader as a bare p-value.
            ## Count on the rows the MODEL used. `tidyr::complete()` above
            ## re-introduces every donor that failed the cell filter as
            ## n_total = NA, so counting `sub` directly would tally donors the
            ## fit never saw -- the same error as P1-1, one level down. COPD
            ## then looked like a normal-sized group and its single-donor
            ## contrast went unflagged.
            fitted_rows <- sub |> filter(!is.na(frac), !is.na(n_total))
            grp_n <- fitted_rows |> distinct(donor_id, disease_group) |>
                count(disease_group, name = "n_donors_group")
            small <- grp_n$disease_group[grp_n$n_donors_group < 3]
            write_tsv_safe(as.data.frame(emm) |>
                               left_join(grp_n, by = "disease_group") |>
                               mutate(min_cells = min_cells_per_donor,
                                      n_donors = fit_n(fit), arm = arm_sfx,
                                      estimable = n_donors_group >= 3),
                           file.path(outdir, paste0("composition_", tag, "_", key,
                                                    "_emmeans", sfx, arm_sfx, ".tsv")))
            write_tsv_safe(posthoc_with_ci(emm) |>
                               mutate(min_cells = min_cells_per_donor, arm = arm_sfx,
                                      touches_small_group = Reduce(`|`,
                                          lapply(small, function(x)
                                              grepl(x, contrast, fixed = TRUE)), FALSE)),
                           file.path(outdir, paste0("composition_", tag, "_", key,
                                                    "_posthoc", sfx, arm_sfx, ".tsv")))
            ## Refit without any <3-donor group so the omnibus can be read against
            ## a version that no single donor can carry.
            p_nosmall <- NA_real_
            if (length(small)) {
                sub2 <- fitted_rows |> filter(!disease_group %in% small) |>
                    mutate(disease_group = droplevels(disease_group))
                if (nlevels(sub2$disease_group) >= 2)
                    p_nosmall <- tryCatch(
                        disease_omnibus(fit_model(covars, "frac", sub2))[["Pr..F."]],
                        error = function(e) NA_real_)
            }
            results[[g]] <- data.frame(
                level = g, n_donors = fit_n(fit), arm = arm_sfx,
                disease_omnibus(fit),
                p_excl_small_groups = p_nosmall,
                small_groups = paste(small, collapse = ","))
        }
        anova_all <- bind_rows(results)
        if (!nrow(anova_all)) return(invisible(NULL))
        pcol <- grep("^Pr", names(anova_all), value = TRUE)[1]
        if (!is.na(pcol)) anova_all$p_BH <- p.adjust(anova_all[[pcol]], method = "BH")
        anova_all$min_cells <- min_cells_per_donor
        write_tsv_safe(anova_all, file.path(outdir, paste0("composition_", tag,
                                            "_disease_anova_all", sfx, arm_sfx, ".tsv")))
        invisible(anova_all)
    }
    a_pri <- run_arm(COVARS_PRIMARY, "", drop_age_rows = FALSE)
    a_sen <- run_arm(COVARS_AGE_SENS, "_ageadj", drop_age_rows = TRUE)
    if (!is.null(a_pri) && !is.null(a_sen))
        cat(sprintf("  [%s] primary n=%s | age-restricted sensitivity n=%s\n", tag,
                    paste(unique(a_pri$n_donors), collapse = "/"),
                    paste(unique(a_sen$n_donors), collapse = "/")))

    p <- ggboxplot(comp, x = "disease_group", y = "frac", add = "jitter",
                   fill = "disease_group", palette = "jco",
                   add.params = list(alpha = 0.4, size = 0.8),
                   xlab = "", ylab = "Fraction per donor",
                   legend = "none", ggtheme = theme_pubr(base_size = 12)) +
        facet_wrap(vars(.data[[group]]), scales = "free_y") + rotate_x_text(35)
    save_ggplots(file.path(outdir, paste0("composition_", tag, "_by_disease", sfx)), p, 9, 7)
    invisible(comp)
}

## ----- (C) Injury-program fraction vs disease (headline) ------------------
injury_fraction_by_disease <- function(df, outdir, min_cells_per_donor = 10, sfx = "") {
    ## Hard fail if a named injury program is absent from the data. Silently
    ## matching nothing is exactly how this endpoint drifted from three programs to
    ## one without any output changing shape.
    have <- levels(factor(df$state_program))
    missing <- setdiff(INJURY_PROGRAMS, have)
    if (length(missing))
        stop("INJURY_PROGRAMS names levels absent from state_program: ",
             paste(missing, collapse = ", "), "\nPresent levels: ",
             paste(have, collapse = ", "),
             "\nUpdate INJURY_PROGRAMS deliberately -- do not leave stale names.")
    cat("injury programs:", paste(INJURY_PROGRAMS, collapse = " + "),
        "| min_cells =", min_cells_per_donor, "\n")

    donor_tot <- df |> count(donor_id, name = "n_total") |>
        filter(n_total >= min_cells_per_donor)
    donor_meta <- df |> group_by(donor_id) |>
        summarise(disease_group = first(disease_group), sex = first(sex),
                  dataset = first(dataset), study = first(study),
                  age = mean(age, na.rm = TRUE),
                  .groups = "drop")
    inj <- df |>
        semi_join(donor_tot, by = "donor_id") |>
        mutate(is_injury = state_program %in% INJURY_PROGRAMS) |>
        group_by(donor_id) |>
        summarise(injury_frac = mean(is_injury), n_injury = sum(is_injury),
                  n_total = n(), .groups = "drop") |>
        left_join(donor_meta, by = "donor_id") |>
        tidyr::drop_na(sex) |>
        mutate(disease_group = relevel(droplevels(factor(disease_group)), "Healthy"))

    write_tsv_safe(
        inj |> select(donor_id, disease_group, dataset, sex, age,
                      n_injury, n_total, injury_frac) |>
            mutate(min_cells = min_cells_per_donor) |>
            arrange(disease_group, donor_id),
        file.path(outdir, paste0("injury_fraction_by_donor", sfx, ".tsv")))

    ## Primary on every donor; age-restricted arm as an explicit sensitivity.
    ## This is the project's headline disease endpoint, so the two must be
    ## emitted together and the fitted N carried on both.
    ## Match the fibrotic level against the levels actually present rather than
    ## a hardcoded string: the level is `Fibrotic_ILD` here, and a literal
    ## "Fibrotic/ILD" silently counted zero.
    fibro_lvl <- grep("fibro", levels(inj$disease_group), value = TRUE,
                      ignore.case = TRUE)
    if (!length(fibro_lvl))
        warning("no fibrotic level found in disease_group; levels are: ",
                paste(levels(inj$disease_group), collapse = ", "), call. = FALSE)

    run_inj <- function(covars, arm_sfx, drop_age_rows) {
        dat <- if (drop_age_rows) tidyr::drop_na(inj, age) else inj
        dat <- dat |> mutate(disease_group = droplevels(disease_group))
        if (nlevels(dat$disease_group) < 2) {
            warning("injury_fraction[", arm_sfx, "]: <2 disease groups after ",
                    "filtering; arm skipped", call. = FALSE)
            return(invisible(NULL))
        }
        fit <- fit_model(covars, "injury_frac", dat)
        emm <- emmeans(fit, ~ disease_group)
        n_fit <- fit_n(fit)
        ## Per-group donor counts travel with the estimates. Without them a
        ## single-donor group (COPD is n = 1 here) produces a significant
        ## contrast -- Healthy vs COPD p = 0.036 -- with nothing on the row to
        ## warn the reader that it rests on one donor.
        grp_n <- dat |> distinct(donor_id, disease_group) |>
            count(disease_group, name = "n_donors_group")
        emm_df <- as.data.frame(emm) |>
            left_join(grp_n, by = "disease_group") |>
            mutate(min_cells = min_cells_per_donor, n_donors = n_fit,
                   arm = arm_sfx,
                   n_fibrotic = sum(dat$disease_group %in% fibro_lvl, na.rm = TRUE),
                   estimable = n_donors_group >= 3,
                   injury_programs = paste(INJURY_PROGRAMS, collapse = "+"))
        write_tsv_safe(emm_df, file.path(outdir, paste0("injury_fraction_emmeans",
                                                        sfx, arm_sfx, ".tsv")))
        ## Same for contrasts: flag any comparison touching a <3-donor group.
        small <- grp_n$disease_group[grp_n$n_donors_group < 3]
        ph <- posthoc_with_ci(emm) |>
            mutate(min_cells = min_cells_per_donor, arm = arm_sfx, n_donors = n_fit,
                   touches_small_group = Reduce(`|`, lapply(small, function(g)
                       grepl(g, contrast, fixed = TRUE)), FALSE))
        write_tsv_safe(ph, file.path(outdir, paste0("injury_fraction_posthoc", sfx,
                                                    arm_sfx, ".tsv")))
        if (length(small))
            cat(sprintf("    NOTE: %s has <3 donors; its contrasts are flagged\n",
                        paste(small, collapse = ", ")))
        cat(sprintf("  injury_fraction[%s]: n=%s donors, %d fibrotic\n",
                    if (nzchar(arm_sfx)) arm_sfx else "primary", n_fit,
                    sum(dat$disease_group %in% fibro_lvl, na.rm = TRUE)))
        invisible(fit)
    }
    fit <- run_inj(COVARS_PRIMARY, "", drop_age_rows = FALSE)
    run_inj(COVARS_AGE_SENS, "_ageadj", drop_age_rows = TRUE)
    emm <- emmeans(fit, ~ disease_group)
    ylab_txt <- paste0("Injury-program fraction\n(",
                       paste(INJURY_PROGRAMS, collapse = " + "), " states)")
    p <- ggboxplot(inj, x = "disease_group", y = "injury_frac", add = "jitter",
                   fill = "disease_group", palette = "jco",
                   xlab = "", ylab = ylab_txt,
                   legend = "none", ggtheme = theme_pubr(base_size = 13)) +
        rotate_x_text(35) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black")
    save_ggplots(file.path(outdir, paste0("injury_fraction_by_disease", sfx)), p, 5, 5)
    invisible(inj)
}

## ----- Main ---------------------------------------------------------------
df <- load_meta()
cat("cells:", nrow(df), " donors:", nlevels(df$donor_id), "\n")
cat("stable states:\n"); print(table(df$pericyte_state))
cat("programs:\n"); print(table(df$state_program))
print(table(df$disease_group))

outdir <- "stats_data"
if (!dir.exists(outdir)) dir.create(outdir)

agtr1_by_group(df, "state_program", outdir, "program")
agtr1_by_group(df, "pericyte_state", outdir, "state")
## Primary threshold first (canonical, unsuffixed filenames); the rest are
## sensitivity fits written alongside with a `_mincells<N>` suffix.
for (i in seq_along(MIN_CELLS)) {
    mc  <- MIN_CELLS[i]
    sfx <- if (i == 1L) "" else paste0("_mincells", mc)
    cat("\n===== donor filter: >=", mc, "pericytes",
        if (i == 1L) "(PRIMARY)" else "(sensitivity)", "=====\n")
    composition_by_disease(df, "pericyte_state", outdir, "state", mc, sfx)
    composition_by_disease(df, "state_program", outdir, "program", mc, sfx)
    injury_fraction_by_disease(df, outdir, mc, sfx)
}

cat("\nReproducibility information:\n")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
