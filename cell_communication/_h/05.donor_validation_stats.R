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

results <- list()
for (pred in c("sender_TGFB1", "sender_ligand_mean")) {
    sub <- d |> tidyr::drop_na(receiver_target_expr, all_of(pred))
    ct <- cor.test(sub[[pred]], sub$receiver_target_expr, method = "spearman")
    ## donor-level model, disease-adjusted; dataset random effect if available
    if ("dataset" %in% names(sub) && length(unique(sub$dataset)) > 1) {
        fit <- suppressMessages(lmerTest::lmer(
            reformulate(c(pred, "disease_group", "(1 | dataset)"), "receiver_target_expr"),
            data = sub))
        co <- summary(fit)$coefficients[pred, ]
        est <- co["Estimate"]; se <- co["Std. Error"]; p <- co["Pr(>|t|)"]; model <- "lmer(+1|dataset)"
    } else {
        fit <- lm(reformulate(c(pred, "disease_group"), "receiver_target_expr"), data = sub)
        co <- summary(fit)$coefficients[pred, ]
        est <- co["Estimate"]; se <- co["Std. Error"]; p <- co["Pr(>|t|)"]; model <- "lm"
    }
    results[[pred]] <- data.frame(
        predictor = pred, n_donors = nrow(sub), model = model,
        spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
        adj_estimate = unname(est), adj_se = unname(se), adj_p = unname(p))
}
res <- bind_rows(results)
res$adj_p_BH <- p.adjust(res$adj_p, method = "BH")
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
