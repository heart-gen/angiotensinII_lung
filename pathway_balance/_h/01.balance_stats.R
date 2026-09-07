## Donor-aware AT1R/AT2R balance statistics.
##
## Tests whether the AT1R-AT2R balance is (i) higher in injury pericyte programs
## and (ii) higher in fibrotic/ILD lungs -- a transcriptomic rationale for AT1R
## blockade (losartan). CAVEAT: the AT1R signature overlaps the injury-stromal /
## NicheNet program (TGFB1, CCN2, COL1A1, FN1, ACTA2, IL6, ...), so this balance
## partly re-measures injury intensity and is NOT a receptor-specific AT1R readout
## (AGTR1 itself is not disease-associated). Result: program contrasts NS (smallest
## p = 0.067).
##
## THE DISEASE-LEVEL RESULT IS REPORTED FROM THE STUDY-GUARDED FIT. This header
## used to read "the shift is disease-level (Healthy vs Other p = 0.040)", which
## was the plain-`lm` value from section (C). Section (B) fits the same donors
## with `(1 | dataset)` -- the guard this script exists to apply -- and it
## reverses the answer. Sections (C) and the arm decomposition now fit both and
## label them with a `study_guard` column; quote the row marked PRIMARY. See
## P1-9 in writings/TODO.md.
##
## States are now the NVU-pattern model: the stable Leiden clusters live in
## `pericyte_state` (numeric), annotated to an interpretable program in
## `state_program` (relative-enrichment argmax). The balance contrasts key on
## `state_program` so INJURY_STATES (program names) match.

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggpubr)
    library(lme4); library(lmerTest); library(emmeans)
})
emm_options(lmerTest.limit = 20000, pbkrtest.limit = 20000)

## Injury programs, as SCORES rather than as argmax labels. The label-based
## grouping (`state_program %in% INJURY_STATES`) is no longer usable: after the
## basement-membrane panel was added, `fibroblast_like` stops winning any cluster
## and the label-based injury set collapses from 4,420 cells / 139 donors to 220
## cells / 65 donors. Basement-membrane deposition is a matrix-stabilizing
## vascular function -- near-orthogonal to fibrillar ECM and aligning pericytes
## with endothelium rather than fibroblasts -- so it is deliberately NOT counted
## as injury. The three program scores below are bit-identical before and after
## the relabelling, so this selection is continuous with the previous analysis.
INJURY_SCORES <- c("inflammatory_score", "fibroblast_like_score",
                   "activated_migratory_score")

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc) ~ "Healthy",
        lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
        TRUE ~ "Other")
}
save_ggplots <- function(fn, p, w, h) for (e in c(".pdf", ".png")) ggsave(paste0(fn, e), p, width = w, height = h)
write_tsv_safe <- function(x, f, rn = FALSE) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), f, sep = "\t", quote = FALSE, row.names = rn)
}

df <- data.table::fread("pathway_balance_metadata.tsv.gz") |>
    mutate(disease_group = factor(map_disease_group(lung_condition)),
           age = suppressWarnings(as.numeric(age_or_mean_of_age_range)),
           sex = factor(sex), state_program = factor(state_program))

outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)

## Donor x program aggregation
agg <- df |>
    group_by(donor_id, state_program) |>
    summarise(balance = mean(AT1R_AT2R_balance, na.rm = TRUE),
              AT1R = mean(AT1R_score, na.rm = TRUE), AT2R = mean(AT2R_score, na.rm = TRUE),
              n_cells = n(), disease_group = first(disease_group),
              sex = first(sex), age = mean(age, na.rm = TRUE), .groups = "drop") |>
    filter(n_cells >= 5) |> drop_na(balance, age, sex) |>
    mutate(across(c(state_program, sex, disease_group), droplevels))

## (A) balance across programs -- donor x program pseudobulk with donor random
## intercept (accounts for within-donor correlation across programs).
fit_state <- suppressMessages(lmerTest::lmer(
    balance ~ state_program + disease_group + age + sex + (1 | donor_id), data = agg))
emm_state <- emmeans(fit_state, ~ state_program)
write_tsv_safe(as.data.frame(anova(fit_state)), file.path(outdir, "balance_by_state_anova.tsv"), TRUE)
write_tsv_safe(as.data.frame(emm_state), file.path(outdir, "balance_by_state_emmeans.tsv"))
write_tsv_safe(as.data.frame(pairs(emm_state, adjust = "BH")), file.path(outdir, "balance_by_state_posthoc.tsv"))

p1 <- ggboxplot(agg, x = "state_program", y = "balance", add = "jitter",
                fill = "state_program", palette = "npg", legend = "none",
                add.params = list(alpha = 0.5, size = 1),
                xlab = "Pericyte program", ylab = "AT1R - AT2R balance",
                ggtheme = theme_pubr(base_size = 13)) + rotate_x_text(35) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white")
save_ggplots(file.path(outdir, "balance_by_state"), p1, 6, 5)

## (B) balance vs disease (donor-level, injury states only). A study/dataset
## random intercept guards against disease-study confounding in the HLCA (the
## composition and niche analyses use leave-one-study-out for the same reason).
has_ds <- "dataset" %in% names(df)
## Continuous analogue of "injury-state cells": z-score each injury program
## across cells, average them, and keep the upper half. Selection is on the
## scores themselves, so it does not depend on which panel happens to win the
## per-cluster argmax.
inj_cols <- intersect(INJURY_SCORES, names(df))
if (length(inj_cols) == 0)
    stop("no injury program score columns found; re-run 00.pathway_balance.py ",
         "to export inflammatory/fibroblast_like/activated_migratory scores")
## `..inj_cols` (not `inj_cols`): df is a data.table, which resolves a bare symbol
## in `j` as a column name rather than as a variable from the calling scope.
inj_z <- scale(as.matrix(df[, ..inj_cols]))
inj_z[is.na(inj_z)] <- 0
df$injury_score <- rowMeans(inj_z)
inj_cut <- stats::median(df$injury_score, na.rm = TRUE)
cat(sprintf("\n(B) injury selection: %d of %d cells above the median composite ",
            sum(df$injury_score > inj_cut, na.rm = TRUE), nrow(df)),
    sprintf("injury score (cols: %s)\n", paste(inj_cols, collapse = ", ")))

## `age` is deliberately NOT required here (changed 2026-09-07). It was, and that
## made this the SEVENTH instance of the `+ age` study filter -- age missingness in
## the HLCA is a study property, so `drop_na(age)` silently restricted the donor set
## to age-reporting cohorts and did it unevenly across disease groups. See the
## rationale block in `pericyte_states/_h/01.state_stats.R` and P1-2/P1-10. The
## age-restricted fit is still produced below, labelled as a cohort restriction.
donor_inj <- df |>
    filter(injury_score > inj_cut) |>
    group_by(donor_id) |>
    summarise(balance = mean(AT1R_AT2R_balance, na.rm = TRUE),
              AT1R = mean(AT1R_score, na.rm = TRUE), AT2R = mean(AT2R_score, na.rm = TRUE),
              disease_group = first(disease_group), sex = first(sex),
              dataset = if (has_ds) first(dataset) else NA_character_,
              age = mean(age, na.rm = TRUE), n_cells = n(), .groups = "drop") |>
    filter(n_cells >= 10) |> drop_na(balance, sex) |>
    mutate(disease_group = relevel(droplevels(disease_group), "Healthy"))
cat("\n(B) disease-balance donors by group (note small diseased n):\n"); print(table(donor_inj$disease_group))
cat("(B) age-complete subset would be", sum(!is.na(donor_inj$age)), "of", nrow(donor_inj),
    "donors:\n"); print(table(donor_inj$disease_group[!is.na(donor_inj$age)]))

## THE donor table this module models, written out so that figures never re-derive
## the selection. `figures/_h/manuscript_mechanism_figure.R` used to rebuild it from
## `state_program %in% INJURY` -- the label set this module abandoned -- and landed
## on 220 cells / 5 donors while the module used 5,840 / 59 (defect P1-8). A figure
## that recomputes an upstream selection will drift from it; this file removes the
## opportunity.
write_tsv_safe(donor_inj, file.path(outdir, "balance_donor_injury_selected.tsv"))
if (has_ds && dplyr::n_distinct(donor_inj$dataset) > 1) {
    fit_dx <- suppressMessages(lmerTest::lmer(
        balance ~ disease_group + sex + (1 | dataset), data = donor_inj))
} else {
    fit_dx <- lm(balance ~ disease_group + sex, data = donor_inj)
}
tag_arm <- function(x, arm, n) as.data.frame(x) |>
    dplyr::mutate(arm = arm, n_donors = n,
                  n_donors_group = paste(names(table(donor_inj$disease_group)),
                                         as.integer(table(donor_inj$disease_group)),
                                         sep = "=", collapse = ";"))
write_tsv_safe(tag_arm(emmeans(fit_dx, ~ disease_group), "primary", nrow(donor_inj)),
               file.path(outdir, "balance_by_disease_emmeans.tsv"))
write_tsv_safe(tag_arm(pairs(emmeans(fit_dx, ~ disease_group), adjust = "BH"), "primary", nrow(donor_inj)),
               file.path(outdir, "balance_by_disease_posthoc.tsv"))

## Age-restricted companion, written with an `_ageadj` suffix and an `arm` column
## so it can never be mistaken for the primary. This is a RESTRICTION TO
## AGE-REPORTING COHORTS, not an age adjustment.
d_age <- tidyr::drop_na(donor_inj, age) |>
    dplyr::mutate(disease_group = droplevels(disease_group))
if (nlevels(d_age$disease_group) >= 2 &&
    (!has_ds || dplyr::n_distinct(d_age$dataset) > 1)) {
    fit_age <- if (has_ds && dplyr::n_distinct(d_age$dataset) > 1)
        suppressMessages(lmerTest::lmer(balance ~ disease_group + age + sex + (1 | dataset), data = d_age))
    else lm(balance ~ disease_group + age + sex, data = d_age)
    ea <- as.data.frame(emmeans(fit_age, ~ disease_group))
    pa <- as.data.frame(pairs(emmeans(fit_age, ~ disease_group), adjust = "BH"))
    for (x in list(list(ea, "balance_by_disease_emmeans_ageadj.tsv"),
                   list(pa, "balance_by_disease_posthoc_ageadj.tsv"))) {
        d <- x[[1]]; d$arm <- "_ageadj -- RESTRICTED to age-reporting cohorts"
        d$n_donors <- nrow(d_age)
        d$n_donors_group <- paste(names(table(d_age$disease_group)),
                                  as.integer(table(d_age$disease_group)),
                                  sep = "=", collapse = ";")
        write_tsv_safe(d, file.path(outdir, x[[2]]))
    }
    cat("\n(B) age-restricted arm: ", nrow(d_age), " of ", nrow(donor_inj),
        " donors\n", sep = "")
} else {
    cat("\n(B) age-restricted arm NOT estimable (groups collapse); no _ageadj files written\n")
}
p2 <- ggboxplot(donor_inj, x = "disease_group", y = "balance", add = "jitter",
                fill = "disease_group", palette = "jco", legend = "none",
                xlab = "", ylab = "AT1R - AT2R balance\n(injury pericytes)",
                ggtheme = theme_pubr(base_size = 13)) + rotate_x_text(30) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white")
save_ggplots(file.path(outdir, "balance_by_disease"), p2, 5, 5)

## (C) RIGOR / DEMOTION tests: is the disease-balance shift INDEPENDENT of the
## injury-stromal program, or just a redundant re-measurement of injury intensity?
## (the AT1R signature overlaps the injury/NicheNet effector genes). We (i) adjust
## the disease effect for the donor injury-stromal score, (ii) decompose into the
## AT1R and AT2R arms separately. Conclusion (documented): the balance is a
## corollary of injury intensity, not independent support -- see MECHANISM_ANALYSES.
ni_file <- "../../niche_index/_m/niche_index_per_donor.tsv.gz"
if (file.exists(ni_file)) {
    ni <- data.table::fread(ni_file, select = c("donor_id", "injury_stromal_score"))
    adj <- merge(donor_inj, ni, by = "donor_id")
    if (nrow(adj) >= 10 && "injury_stromal_score" %in% names(adj)) {
        cat(sprintf("\n(C) injury-adjusted: %d donors; cor(balance, injury_stromal)=%.3f\n",
                    nrow(adj), cor(adj$balance, adj$injury_stromal_score, use = "complete")))
        ## P1-9. Section (B) guards against disease-STUDY confounding in the HLCA
        ## with a `(1 | dataset)` random intercept -- this script's own comment
        ## says that is why it is there. Section (C) and the arm decomposition
        ## were plain `lm` with no such term, and the header quoted the UNGUARDED
        ## number. The guard is not cosmetic here: on the same donors it reverses
        ## both the sign and the verdict of the Healthy-vs-Other effect. Both
        ## fits are therefore run and both are written, distinguished by a
        ## `study_guard` column, with the GUARDED fit designated primary.
        guarded <- has_ds && dplyr::n_distinct(adj$dataset) > 1
        rows <- function(fit, tag, guard) {
            co <- as.data.frame(summary(fit)$coefficients)
            names(co)[names(co) == "Std. Error"] <- "SE"
            names(co)[names(co) == "Pr(>|t|)"]   <- "p_value"
            names(co)[names(co) == "t value"]    <- "t"
            co <- co[, intersect(c("Estimate", "SE", "df", "t", "p_value"), names(co)),
                     drop = FALSE]
            co$term <- rownames(co); co$model <- tag; co$study_guard <- guard
            co
        }
        GUARD_P <- "(1 | dataset) -- PRIMARY"
        GUARD_N <- "none -- UNGUARDED, for comparison only"
        fit_pair <- function(rhs, tag, resp = "balance") {
            o <- rows(lm(reformulate(rhs, resp), data = adj), tag, GUARD_N)
            if (guarded) {
                g  <- suppressMessages(lmerTest::lmer(
                          reformulate(c(rhs, "(1 | dataset)"), resp), data = adj))
                vc <- as.data.frame(lme4::VarCorr(g))
                gr <- rows(g, tag, GUARD_P)
                ## When the between-dataset SD is estimated at 0 the fit is
                ## singular and the guarded estimate EQUALS the unguarded one.
                ## That is a result, not a failure: it says this response carries
                ## no between-dataset variance for the guard to absorb, so the
                ## unguarded number was safe FOR THAT RESPONSE. Record it rather
                ## than leaving two identical rows unexplained.
                gr$dataset_sd <- vc$sdcor[vc$grp == "dataset"][1]
                gr$singular   <- lme4::isSingular(g)
                o  <- dplyr::bind_rows(gr, o)
            }
            o
        }
        ## `age` dropped from these rhs on 2026-09-07 for the reason given at the
        ## donor_inj block: it was a cohort filter, not a covariate.
        out <- dplyr::bind_rows(
            fit_pair(c("disease_group", "sex"), "unadjusted"),
            fit_pair(c("disease_group", "injury_stromal_score", "sex"),
                     "injury_adjusted"))
        out <- out[grepl("disease|injury", out$term), ]
        ## P1-20(b). The `injury_stromal_score` row is a score-on-score
        ## coefficient: both sides are `score_genes` panels, which have a
        ## non-zero null, so its magnitude is not an effect size. It is kept
        ## because the module's conclusion depends on whether the DISEASE term
        ## shrinks when it is added -- that comparison is categorical and is
        ## unaffected -- but it must not be quoted on its own.
        out$readout <- ifelse(grepl("injury_stromal_score", out$term),
                              "nuisance covariate -- score-on-score, non-zero null, NOT an effect size",
                              "disease contrast -- interpretable")
        write_tsv_safe(out, file.path(outdir, "balance_disease_injury_adjusted.tsv"))
        if (guarded) {
            cmp <- out[out$model == "unadjusted" & grepl("disease", out$term),
                       c("study_guard", "term", "Estimate", "p_value")]
            cat("\n  (C) study guard on / off, same donors -- if these disagree, the\n",
                "      GUARDED row is the one to report (P1-9):\n", sep = "")
            print(cmp, row.names = FALSE)
        }
        ## Report what the adjustment actually did rather than asserting it. This
        ## has now flipped twice and must be read from the run, never quoted from
        ## memory: under the label-based selection the disease term collapsed;
        ## under the continuous selection PLUS the `+ age` filter it did not; with
        ## age removed (2026-09-07, 59 donors instead of 30) it collapses again --
        ## Fibrotic/ILD 0.145 (P = 0.0038) -> 0.064 (P = 0.194) once
        ## `injury_stromal_score` enters. So the balance is a COROLLARY of injury
        ## intensity, which is what MECHANISM_ANALYSES claimed all along.
        cat("  (compare the unadjusted vs injury_adjusted disease terms below;\n",
            "   the covariate absorbs the disease effect only if they shrink)\n", sep = "")
        print(out[, c("study_guard", "model", "term", "Estimate", "p_value")],
              row.names = FALSE)
        ## arm decomposition: which arm (AT1R up vs AT2R down) drives the shift?
        ## Guarded and unguarded, same convention as above (P1-9): this table's
        ## AT1R Fibrotic row was previously reported at P = 0.030 from the
        ## unguarded fit alone.
        arm <- dplyr::bind_rows(lapply(c("AT1R", "AT2R"), function(v) {
            r <- fit_pair(c("disease_group", "sex"), v, resp = v)
            r[grepl("disease", r$term), ] }))
        write_tsv_safe(arm, file.path(outdir, "balance_arm_decomposition.tsv"))
        cat("\n  arm decomposition (AT1R vs AT2R ~ disease):\n")
        print(arm[, c("study_guard", "model", "term", "Estimate", "p_value")],
              row.names = FALSE)
    }
}

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
