## Donor-aware statistics on the de novo pericyte states (NVU pattern).
##
## States are STABLE Leiden clusters on the study-integrated embedding
## (`pericyte_state`), annotated to a dominant curated program (`state_program`)
## by 00.state_discovery.py. Study is handled once, by the integration that the
## clustering runs on, so the donor-level models below do NOT add a study term.
## The unit of replication is the donor throughout.
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
## clear the >=20-cell filter; the powered contrast is Healthy vs Fibrotic/ILD and
## even that is donor-limited. The COPD AGTR1 signal is carried by the whole-stroma
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

INJURY_PROGRAMS <- c("inflammatory", "fibroblast_like", "activated_migratory")

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

## ----- (A) AGTR1 across states / programs (donor x group mixed model) -----
agtr1_by_group <- function(df, group, outdir, tag, min_cells = 5) {
    agg <- df |>
        group_by(donor_id, .data[[group]]) |>
        summarise(AGTR1_mean = mean(AGTR1_expr, na.rm = TRUE), n_cells = n(),
                  disease_group = first(disease_group), sex = first(sex),
                  age = mean(age, na.rm = TRUE), .groups = "drop") |>
        filter(n_cells >= min_cells) |>
        tidyr::drop_na(AGTR1_mean, age, sex)
    agg[[group]]       <- droplevels(factor(agg[[group]]))
    agg$disease_group  <- droplevels(agg$disease_group)
    agg$sex            <- droplevels(agg$sex)
    if (nlevels(agg[[group]]) < 2) return(invisible(NULL))

    # Donor random intercept accounts for within-donor correlation across groups.
    form <- reformulate(c(group, "disease_group", "sex", "age", "(1 | donor_id)"),
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
composition_by_disease <- function(df, group, outdir, tag, min_cells_per_donor = 20) {
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
                  age = mean(age, na.rm = TRUE), .groups = "drop")
    comp <- comp |> left_join(donor_meta, by = "donor_id") |>
        mutate(disease_group = relevel(factor(disease_group), ref = "Healthy"))

    results <- list()
    for (g in levels(factor(comp[[group]]))) {
        sub <- comp |> filter(.data[[group]] == g) |> tidyr::drop_na(age, sex) |>
            mutate(disease_group = droplevels(disease_group))
        if (nlevels(sub$disease_group) < 2) next
        ## Skip degenerate levels (e.g. a single-program grouping makes frac==1
        ## everywhere -> zero residual variance -> Anova.lm errors).
        if (sd(sub$frac, na.rm = TRUE) < 1e-9) next
        fit <- lm(frac ~ disease_group + age + sex, data = sub)
        emm <- emmeans(fit, ~ disease_group)
        key <- sanitize(g)
        write_tsv_safe(as.data.frame(emm),
                       file.path(outdir, paste0("composition_", tag, "_", key, "_emmeans.tsv")))
        write_tsv_safe(as.data.frame(pairs(emm, adjust = "BH")),
                       file.path(outdir, paste0("composition_", tag, "_", key, "_posthoc.tsv")))
        results[[g]] <- data.frame(level = g,
                                   as.data.frame(car::Anova(fit, type = 2))["disease_group", ])
    }
    anova_all <- bind_rows(results)
    pcol <- grep("^Pr", names(anova_all), value = TRUE)[1]
    if (!is.na(pcol)) anova_all$p_BH <- p.adjust(anova_all[[pcol]], method = "BH")
    write_tsv_safe(anova_all, file.path(outdir, paste0("composition_", tag, "_disease_anova_all.tsv")))

    p <- ggboxplot(comp, x = "disease_group", y = "frac", add = "jitter",
                   fill = "disease_group", palette = "jco",
                   add.params = list(alpha = 0.4, size = 0.8),
                   xlab = "", ylab = "Fraction per donor",
                   legend = "none", ggtheme = theme_pubr(base_size = 12)) +
        facet_wrap(vars(.data[[group]]), scales = "free_y") + rotate_x_text(35)
    save_ggplots(file.path(outdir, paste0("composition_", tag, "_by_disease")), p, 9, 7)
    invisible(comp)
}

## ----- (C) Injury-program fraction vs disease (headline) ------------------
injury_fraction_by_disease <- function(df, outdir, min_cells_per_donor = 20) {
    donor_tot <- df |> count(donor_id, name = "n_total") |>
        filter(n_total >= min_cells_per_donor)
    donor_meta <- df |> group_by(donor_id) |>
        summarise(disease_group = first(disease_group), sex = first(sex),
                  age = mean(age, na.rm = TRUE), .groups = "drop")
    inj <- df |>
        semi_join(donor_tot, by = "donor_id") |>
        mutate(is_injury = state_program %in% INJURY_PROGRAMS) |>
        group_by(donor_id) |>
        summarise(injury_frac = mean(is_injury), .groups = "drop") |>
        left_join(donor_meta, by = "donor_id") |>
        tidyr::drop_na(age, sex) |>
        mutate(disease_group = relevel(droplevels(factor(disease_group)), "Healthy"))

    fit <- lm(injury_frac ~ disease_group + age + sex, data = inj)
    emm <- emmeans(fit, ~ disease_group)
    write_tsv_safe(as.data.frame(emm), file.path(outdir, "injury_fraction_emmeans.tsv"))
    write_tsv_safe(as.data.frame(pairs(emm, adjust = "BH")),
                   file.path(outdir, "injury_fraction_posthoc.tsv"))
    p <- ggboxplot(inj, x = "disease_group", y = "injury_frac", add = "jitter",
                   fill = "disease_group", palette = "jco",
                   xlab = "", ylab = "Injury-program fraction\n(inflammatory+fibroblast+activated states)",
                   legend = "none", ggtheme = theme_pubr(base_size = 13)) +
        rotate_x_text(35) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black")
    save_ggplots(file.path(outdir, "injury_fraction_by_disease"), p, 5, 5)
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
composition_by_disease(df, "pericyte_state", outdir, "state")
composition_by_disease(df, "state_program", outdir, "program")
injury_fraction_by_disease(df, outdir)

cat("\nReproducibility information:\n")
Sys.time(); proc.time()
options(width = 120)
sessioninfo::session_info()
