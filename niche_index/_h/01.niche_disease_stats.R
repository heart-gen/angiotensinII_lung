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

MIN_CELLS <- if ("min_cells" %in% names(df)) unique(df$min_cells)[1] else NA_integer_
cat("input:", infile, "| min_cells:", MIN_CELLS, "| donors:", nrow(df), "\n")
print(table(df$disease_group))

outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)

run_one <- function(response) {
    sub <- df |> tidyr::drop_na(all_of(response), age, sex) |>
        mutate(disease_group = droplevels(disease_group))
    fit <- lm(reformulate(c("disease_group", "age", "sex"), response), data = sub)
    robust <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC3"))
    emm <- emmeans(fit, ~ disease_group)
    ## n_donors and min_cells travel with every output: n was previously only
    ## recoverable by adding the residual df back to the parameter count.
    tag <- function(x) as.data.frame(x) |>
        mutate(n_donors = nrow(sub), min_cells = MIN_CELLS,
               model = "lm(~ disease_group + age + sex)")
    ## The ANOVA is written with row names (term labels), so it cannot go through
    ## tag(); add the same provenance columns by hand. Without them a downstream
    ## table cannot tell which donor threshold an F statistic came from.
    av <- as.data.frame(car::Anova(fit, type = 2))
    av$n_donors <- nrow(sub); av$min_cells <- MIN_CELLS
    av$model <- "lm(~ disease_group + age + sex)"
    write_tsv_safe(av, file.path(outdir, paste0(response, "_anova", SFX, ".tsv")), TRUE)
    write_tsv_safe(tag(emm), file.path(outdir, paste0(response, "_emmeans", SFX, ".tsv")))
    write_tsv_safe(tag(pairs(emm, adjust = "BH")),
                   file.path(outdir, paste0(response, "_posthoc", SFX, ".tsv")))
    rb <- as.data.frame(unclass(robust)); rb$term <- rownames(rb)
    write_tsv_safe(tag(rb), file.path(outdir, paste0(response, "_robust_coefs", SFX, ".tsv")))

    p <- ggboxplot(sub, x = "disease_group", y = response, add = "jitter",
                   fill = "disease_group", palette = "jco",
                   add.params = list(alpha = 0.5, size = 1.2),
                   xlab = "", ylab = response, legend = "none",
                   ggtheme = theme_pubr(base_size = 13)) +
        rotate_x_text(30) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                     fill = "white", color = "black")
    save_ggplots(file.path(outdir, paste0("box_", response, SFX)), p, 5, 5)
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
