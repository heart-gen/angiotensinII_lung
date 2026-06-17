## Donor-level disease association of the pericyte-endothelial niche index.
##
## Hypothesis (mirrors the in-vivo mouse pericyte loss + losartan rescue):
## fibrotic/ILD lungs have LOWER niche-stability and HIGHER injury-stromal
## scores than healthy lungs. Tested with ANCOVA + HC3-robust SEs + emmeans.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(emmeans)
})

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

df <- data.table::fread("niche_index_per_donor.tsv.gz") |>
    mutate(disease_group = relevel(factor(map_disease_group(lung_condition)), "Healthy"),
           sex = factor(sex), age = suppressWarnings(as.numeric(age))) |>
    filter(!is.na(disease_group))

cat("donors:", nrow(df), "\n"); print(table(df$disease_group))

outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)

run_one <- function(response) {
    sub <- df |> tidyr::drop_na(all_of(response), age, sex) |>
        mutate(disease_group = droplevels(disease_group))
    fit <- lm(reformulate(c("disease_group", "age", "sex"), response), data = sub)
    robust <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC3"))
    emm <- emmeans(fit, ~ disease_group)
    write_tsv_safe(as.data.frame(car::Anova(fit, type = 2)),
                   file.path(outdir, paste0(response, "_anova.tsv")), TRUE)
    write_tsv_safe(as.data.frame(emm), file.path(outdir, paste0(response, "_emmeans.tsv")))
    write_tsv_safe(as.data.frame(pairs(emm, adjust = "BH")), file.path(outdir, paste0(response, "_posthoc.tsv")))
    rb <- as.data.frame(unclass(robust)); rb$term <- rownames(rb)
    write_tsv_safe(rb, file.path(outdir, paste0(response, "_robust_coefs.tsv")))

    p <- ggboxplot(sub, x = "disease_group", y = response, add = "jitter",
                   fill = "disease_group", palette = "jco",
                   add.params = list(alpha = 0.5, size = 1.2),
                   xlab = "", ylab = response, legend = "none",
                   ggtheme = theme_pubr(base_size = 13)) +
        rotate_x_text(30) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black")
    save_ggplots(file.path(outdir, paste0("box_", response)), p, 5, 5)
    invisible(fit)
}

## PRIMARY responses (main narrative + main figures).
for (resp in c("niche_stability_score", "injury_stromal_score", "niche_index"))
    run_one(resp)

## SENSITIVITY responses: injury / index composites that additionally fold in the
## AGTR1+ pericyte fraction. Reported in the supplement only -- the AGTR1+ fraction
## is kept out of the primary composite to avoid circularity with the focal receptor.
for (resp in c("injury_stromal_score_sens_agtr1", "niche_index_sens_agtr1"))
    if (resp %in% names(df)) run_one(resp)

cat("\nNOTE: A lower niche-index / higher injury-stromal score in fibrotic/ILD\n",
    "lungs is the human transcriptomic counterpart of the in-vivo mouse\n",
    "pericyte-loss phenotype rescued by losartan (AT1 blockade).\n")

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
