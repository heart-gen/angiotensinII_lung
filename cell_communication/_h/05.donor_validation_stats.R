## Donor-level validation of the NicheNet prediction.
##
## Tests, across donors, whether fibroblast/mural TGF-beta (and the top
## prioritized ligand set) predicts the WHOLE-pericyte injury/activation
## target-program expression (receiver = all pericytes, NOT a binary AGTR1+ split,
## which would be dropout-driven) -- a donor-replicated check of the network-based
## NicheNet ranking. Unit = donor (one row each); dataset random effect where
## available. This converts the CCC result from a ranking into an inferential
## donor-level association.

suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(ggpubr)
    library(lme4); library(lmerTest); library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
indir <- if (length(args) >= 1) args[1] else "."
outdir <- file.path(indir, "donor_validation"); dir.create(outdir, showWarnings = FALSE)

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(grepl("^Healthy", lc) ~ "Healthy", lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
        TRUE ~ "Other")
}
write_tsv <- function(x, f) write.table(as.data.frame(x, check.names = FALSE), f,
                                        sep = "\t", quote = FALSE, row.names = FALSE)

d <- data.table::fread(file.path(indir, "donor_validation_table.tsv.gz"))
# disease_group may already be harmonized; normalize defensively
if (!"disease_group" %in% names(d) || all(is.na(d$disease_group)))
    d$disease_group <- map_disease_group(d$lung_condition)
d <- d |> mutate(disease_group = relevel(factor(disease_group), "Healthy"))
cat("donors:", nrow(d), "\n")

## Two donor-level specifications are fit for every predictor, and BOTH are
## published -- they are not competing versions of one analysis:
##
##   lmer(+1|dataset)  dataset as a RANDOM effect. Figure S8C/D.
##   lm(+dataset)      dataset as a FIXED effect. Figure 4D, where the panel is a
##                     Frisch-Waugh-Lovell partial-regression plot and the drawn
##                     residual slope IS the fixed-effect coefficient. A random
##                     effect would shrink toward the grand mean and the annotated
##                     estimate would no longer equal the line in the figure.
##
## They give materially different estimates (composite: 0.649 vs 0.50), so the
## `model` column below is load-bearing: never quote one for the other's figure.
## Figure 4D reads its annotation from this file rather than refitting inline.
fit_one <- function(sub, pred, spec) {
    if (spec == "lmer") {
        fit <- suppressMessages(lmerTest::lmer(
            reformulate(c(pred, "disease_group", "(1 | dataset)"), "receiver_target_expr"),
            data = sub))
        list(co = summary(fit)$coefficients[pred, ], model = "lmer(+1|dataset)",
             figure = "Figure S8C")
    } else {
        fit <- lm(reformulate(c(pred, "dataset", "disease_group"), "receiver_target_expr"),
                  data = sub)
        list(co = summary(fit)$coefficients[pred, ], model = "lm(+dataset fixed)",
             figure = "Figure 4D")
    }
}

results <- list()
for (pred in c("sender_TGFB1", "sender_TGFB2", "sender_ligand_mean")) {
    if (!pred %in% names(d)) { cat("skipping absent predictor:", pred, "\n"); next }
    sub <- d |> tidyr::drop_na(receiver_target_expr, all_of(pred))
    if (nrow(sub) < 10 || all(!is.finite(sub[[pred]]))) {
        cat("skipping predictor with too few usable donors:", pred, "\n"); next
    }
    ct <- cor.test(sub[[pred]], sub$receiver_target_expr, method = "spearman")
    has_ds <- "dataset" %in% names(sub) && length(unique(sub$dataset)) > 1
    ## Without >1 dataset neither adjustment is identifiable; fall back to the
    ## disease-only lm and label it so, rather than silently mislabelling it.
    specs <- if (has_ds) c("lmer", "lm_fixed") else c("lm_fixed")
    for (spec in specs) {
        f <- if (has_ds) fit_one(sub, pred, spec) else
            list(co = summary(lm(reformulate(c(pred, "disease_group"),
                                             "receiver_target_expr"), data = sub))$coefficients[pred, ],
                 model = "lm(no dataset adj)", figure = NA_character_)
        est <- unname(f$co["Estimate"]); se <- unname(f$co["Std. Error"])
        results[[paste(pred, spec)]] <- data.frame(
            predictor = pred, n_donors = nrow(sub), model = f$model, figure = f$figure,
            spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
            adj_estimate = est, adj_se = se,
            adj_ci_lo = est - 1.96 * se, adj_ci_hi = est + 1.96 * se,
            adj_p = unname(f$co[grep("^Pr\\(", names(f$co))]))
    }
}
res <- bind_rows(results)
## BH WITHIN each model family: the specifications are alternative analyses of the
## same donors, so pooling them would adjust each comparison twice over.
res <- res |> group_by(model) |> mutate(adj_p_BH = p.adjust(adj_p, method = "BH")) |>
    ungroup() |> as.data.frame()
write_tsv(res, file.path(outdir, "donor_validation_results.tsv"))
print(res)

## scatter: fibroblast TGFB1 vs whole-pericyte target program (per donor)
ds <- d |> tidyr::drop_na(sender_TGFB1, receiver_target_expr)
p <- ggscatter(ds, x = "sender_TGFB1", y = "receiver_target_expr",
               add = "reg.line", conf.int = TRUE, color = "disease_group",
               palette = c(Healthy = "#0072B2", COPD = "#E69F00",
                           Fibrotic_ILD = "#D55E00", Other = "#999999"),
               add.params = list(color = "black"),
               xlab = "Fibroblast/mural TGFB1 (per donor)",
               ylab = "Pericyte target program (per donor)",
               ggtheme = theme_pubr(base_size = 12)) +
    stat_cor(method = "spearman", size = 3)
for (e in c(".pdf", ".png")) ggsave(file.path(outdir, paste0("donor_validation_scatter", e)),
                                    p, width = 5, height = 4)

cat("\nReproducibility information:\n"); options(width = 120); sessioninfo::session_info()
