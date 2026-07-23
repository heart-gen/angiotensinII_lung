## Donor-aware AT1R/AT2R balance statistics.
##
## Tests whether the AT1R-AT2R balance is (i) higher in injury pericyte programs
## and (ii) higher in fibrotic/ILD lungs -- a transcriptomic rationale for AT1R
## blockade (losartan). CAVEAT: the AT1R signature overlaps the injury-stromal /
## NicheNet program (TGFB1, CCN2, COL1A1, FN1, ACTA2, IL6, ...), so this balance
## partly re-measures injury intensity and is NOT a receptor-specific AT1R readout
## (AGTR1 itself is not disease-associated). Result: program contrasts NS (smallest
## p = 0.067); the shift is disease-level (Healthy vs Other p = 0.040).
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

donor_inj <- df |>
    filter(injury_score > inj_cut) |>
    group_by(donor_id) |>
    summarise(balance = mean(AT1R_AT2R_balance, na.rm = TRUE),
              AT1R = mean(AT1R_score, na.rm = TRUE), AT2R = mean(AT2R_score, na.rm = TRUE),
              disease_group = first(disease_group), sex = first(sex),
              dataset = if (has_ds) first(dataset) else NA_character_,
              age = mean(age, na.rm = TRUE), n_cells = n(), .groups = "drop") |>
    filter(n_cells >= 10) |> drop_na(balance, age, sex) |>
    mutate(disease_group = relevel(droplevels(disease_group), "Healthy"))
cat("\n(B) disease-balance donors by group (note small diseased n):\n"); print(table(donor_inj$disease_group))
if (has_ds && dplyr::n_distinct(donor_inj$dataset) > 1) {
    fit_dx <- suppressMessages(lmerTest::lmer(
        balance ~ disease_group + age + sex + (1 | dataset), data = donor_inj))
} else {
    fit_dx <- lm(balance ~ disease_group + age + sex, data = donor_inj)
}
write_tsv_safe(as.data.frame(emmeans(fit_dx, ~ disease_group)), file.path(outdir, "balance_by_disease_emmeans.tsv"))
write_tsv_safe(as.data.frame(pairs(emmeans(fit_dx, ~ disease_group), adjust = "BH")), file.path(outdir, "balance_by_disease_posthoc.tsv"))
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
        f_un  <- lm(balance ~ disease_group + age + sex, data = adj)
        f_adj <- lm(balance ~ disease_group + injury_stromal_score + age + sex, data = adj)
        rows <- function(fit, tag) { co <- as.data.frame(summary(fit)$coefficients)
            co$term <- rownames(co); co$model <- tag; co }
        out <- rbind(rows(f_un, "unadjusted"), rows(f_adj, "injury_adjusted"))
        out <- out[grepl("disease|injury", out$term), ]
        write_tsv_safe(out, file.path(outdir, "balance_disease_injury_adjusted.tsv"))
        ## Report what the adjustment actually did rather than asserting it: under the
        ## continuous injury selection the disease term does NOT collapse, which is the
        ## opposite of what the earlier label-based selection showed.
        cat("  (compare the unadjusted vs injury_adjusted disease terms below;\n",
            "   the covariate absorbs the disease effect only if they shrink)\n", sep = "")
        print(out[, c("model", "term", "Estimate", "Pr(>|t|)")])
        ## arm decomposition: which arm (AT1R up vs AT2R down) drives the shift?
        arm <- do.call(rbind, lapply(c("AT1R", "AT2R"), function(v) {
            fit <- lm(reformulate(c("disease_group", "age", "sex"), v), data = adj)
            r <- rows(fit, v); r[grepl("disease", r$term), ] }))
        write_tsv_safe(arm, file.path(outdir, "balance_arm_decomposition.tsv"))
        cat("\n  arm decomposition (AT1R vs AT2R ~ disease):\n")
        print(arm[, c("model", "term", "Estimate", "Pr(>|t|)")])
    }
}

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
