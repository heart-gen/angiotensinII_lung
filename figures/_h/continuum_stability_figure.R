## Supplementary figure S7: is the pericyte transcriptional continuum an artifact of
## the analyst's choices?
##
## The headline result is that along a diffusion pseudotime rooted at the
## vascular-stabilizing pole, the injury / activation / ECM program scores rise and
## the basement-membrane score falls. A reviewer can object that this depends on the
## root cell, the neighborhood size, the number of diffusion components, or the
## particular cells included. pericyte_states/_h/02b.continuum_sensitivity.py re-runs
## DPT across all of those choices; this figure shows the resulting distribution.
##
##   A  donor-level rho for each feature across the 18 parameter settings
##   B  the spread of those 18 estimates per feature, with sign consistency
##   C  donor-level rho after re-rooting at each alternative cluster / latent root
##   D  |rho| under the canonical root vs |rho| under each alternative root
##
## The distinction A/B vs C/D matters: parameters should leave sign AND magnitude
## alone, whereas re-rooting at the opposite pole legitimately REVERSES pseudotime,
## so only |rho| is expected to be preserved. D is the panel that makes that explicit.
##
## Correlations are donor-level (Spearman of per-donor mean pseudotime vs per-donor
## mean score) -- the conservative, donor-aware read, aggregated exactly as the
## headline 02.continuum_dpt.py does.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

PS <- P("pericyte_states", "_m")
LEVEL <- "donor"

FEAT_ORDER <- c("vascular_stabilizing_score", "basement_membrane_score",
                "inflammatory_score", "synthetic_contractile_score",
                "activated_migratory_score", "fibroblast_like_score",
                "AGTR1_expr")
FEAT_LABS <- c(vascular_stabilizing_score = "Vascular-stabilizing",
               basement_membrane_score    = "Basement-membrane",
               inflammatory_score         = "Inflammatory",
               synthetic_contractile_score = "Synthetic/contractile",
               activated_migratory_score  = "Activated/migratory",
               fibroblast_like_score      = "Fibroblast-like",
               AGTR1_expr                 = "AGTR1")
feat_factor <- function(x) factor(x, levels = rev(FEAT_ORDER))

ROOT_LABS <- c(pc1min = "PC1 min", pc1max = "PC1 max")
root_lab <- function(x) ifelse(x %in% names(ROOT_LABS), ROOT_LABS[x],
                               ifelse(grepl("^[0-9]+$", x), paste0("Cluster ", x), x))

runs <- fread(file.path(PS, "continuum_sensitivity_runs.tsv"))
if (!"level" %in% names(runs))
    stop("continuum_sensitivity_runs.tsv has no `level` column -- rerun ",
         "pericyte_states/_h/step_2b.sh with the donor-aware 02b script.")
runs <- runs[level == LEVEL & feature %in% FEAT_ORDER]
summ <- fread(file.path(PS, "continuum_sensitivity_summary.tsv"))[
    level == LEVEL & feature %in% FEAT_ORDER]

## Diverging fill shared by the two heatmaps so A and C are directly comparable.
rho_lim <- max(abs(runs$spearman_rho), na.rm = TRUE)
fill_rho <- list(
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-rho_lim, rho_lim),
                         name = expression(Spearman~rho)),
    theme(legend.position = "right", legend.key.width = unit(3, "mm"),
          legend.key.height = unit(8, "mm"), legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6)))

## ===== Panel A -- parameter sweep at the canonical root ===================
prm <- runs[block == "param"]
## Compact column codes: 18 settings will not fit as spelled-out labels, so the
## axis title carries the key and the tick is just the three varying numbers. The
## two 80% draws differ only by seed and are distinguished a/b.
prm[, subsample := fifelse(frac >= 1, "100",
                           paste0(round(frac * 100), fifelse(seed == 13, "a", "b")))]
prm[, setting := sprintf("%d/%d/%s", neighbors, n_dcs, subsample)]
setting_order <- unique(prm[order(neighbors, n_dcs, -frac, seed), setting])
prm[, setting := factor(setting, levels = setting_order)]

pA <- ggplot(prm, aes(setting, feat_factor(feature), fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_y_discrete(labels = FEAT_LABS) +
    labs(x = sprintf("n_neighbors / n_DCs / subsample %% (%d settings)",
                     length(setting_order)), y = NULL) +
    theme_ms() + fill_rho +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.4),
          axis.text.y = element_text(size = 6.5),
          panel.grid = element_blank(), panel.border = element_blank())

## ===== Panel B -- spread of the 18 estimates ==============================
sgn <- summ[, .(feature, sign_consistency, mean_rho,
                lab = sprintf("%.0f%%", 100 * sign_consistency))]
pB <- ggplot(prm, aes(spearman_rho, feat_factor(feature))) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    stat_summary(fun.min = min, fun.max = max, geom = "linerange",
                 linewidth = 0.5, colour = "grey60") +
    geom_point(size = 1.1, alpha = 0.55, colour = "#0072B2",
               position = position_jitter(height = 0.12, width = 0)) +
    geom_text(data = sgn, aes(x = mean_rho, y = feat_factor(feature), label = lab),
              vjust = -1.3, size = 2.0, colour = "grey30") +
    scale_y_discrete(labels = FEAT_LABS) +
    labs(x = expression(Spearman~rho~"(pseudotime vs score)"), y = NULL) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6.5))

## ===== Panel C -- re-rooting ==============================================
rt <- runs[block == "root"]
rt[, root_label := root_lab(root)]
root_order <- unique(rt[order(!grepl("^Cluster", root_label), root_label), root_label])
rt[, root_label := factor(root_label, levels = root_order)]

pC <- ggplot(rt, aes(root_label, feat_factor(feature), fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_y_discrete(labels = FEAT_LABS) +
    labs(x = "Alternative root", y = NULL) +
    theme_ms() + fill_rho +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
          axis.text.y = element_text(size = 6.5),
          panel.grid = element_blank(), panel.border = element_blank())

## ===== Panel D -- |rho| canonical vs |rho| alternative root ===============
## Canonical = the vascular-stabilizing root at the default parameters that the
## re-rooting block also uses (k30 / dc10 / all cells), so the two differ only in root.
canon <- prm[neighbors == 30 & n_dcs == 10 & frac >= 1,
             .(feature, rho_canon = spearman_rho)]
cmp <- merge(rt[, .(feature, root_label, rho_alt = spearman_rho)], canon, by = "feature")
cmp[, `:=`(abs_canon = abs(rho_canon), abs_alt = abs(rho_alt),
           reversed = sign(rho_alt) != sign(rho_canon))]
d_rho <- cor(cmp$abs_canon, cmp$abs_alt, method = "spearman")
ax_max <- max(c(cmp$abs_canon, cmp$abs_alt)) * 1.05

pD <- ggplot(cmp, aes(abs_canon, abs_alt)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = 2,
                linewidth = 0.35) +
    geom_point(aes(colour = feature, shape = reversed), size = 1.5, alpha = 0.85) +
    annotate("text", x = 0, y = ax_max, hjust = 0, vjust = 1, size = 2.3,
             colour = "grey20",
             label = sprintf("Spearman rho(|rho|) = %.3f\nacross %d root x feature pairs",
                             d_rho, nrow(cmp))) +
    scale_colour_manual(values = setNames(OKABE[seq_along(FEAT_ORDER)], FEAT_ORDER),
                        labels = FEAT_LABS, name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), name = NULL,
                       labels = c(`FALSE` = "same sign", `TRUE` = "sign reversed")) +
    coord_equal(xlim = c(0, ax_max), ylim = c(0, ax_max)) +
    labs(x = expression("|"*rho*"|, canonical vascular-stabilizing root"),
         y = expression("|"*rho*"|, alternative root")) +
    theme_ms() +
    theme(legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 5.8),
          legend.spacing.y = unit(1, "mm"))

## ---- assemble ------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_continuum_stability", fig, 9.5, 7.0)

cat("Wrote figureS_continuum_stability to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
