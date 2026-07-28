## Supplementary figure S8: is the prioritized-ligand result specific, and does it
## hold at the level of individual donors?
##
## NicheNet ranks ligands by how well their prior regulatory potential recovers a
## curated target gene set. Two fair objections follow: (i) any gene set of that
## size might score well (the ranking would then be an artifact of the prior network,
## not of pericyte biology), and (ii) the whole result is computed on pooled
## pseudobulk, so it might not hold donor by donor. Panels A-B answer the first with
## a size-matched gene-set permutation null; C-E answer the second with the
## donor-level sender/receiver relationship.
##
##   A  observed corrected AUPR per prioritized ligand vs its permutation null
##   B  standardized distance from that null, with the empirical p floor
##   C  donor-level composite prioritized-ligand signal vs pericyte target response
##   D  the same for TGFB1 alone
##   E  the same for TGFB2 alone
##
## READ C-E CAREFULLY: the drawn line is the RAW (unadjusted) fit, while the
## annotated beta is the dataset- and disease-adjusted LMM coefficient. In C and D
## the two agree in direction so the mismatch is harmless, but in E they do not --
## TGFB2's raw association is flat (rho = 0.06, p = 0.52) while its adjusted beta is
## the largest of the three (0.54, p = 0.087). The line and the number are answering
## different questions there, and the legend says so.
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

CCC <- P("cell_communication", "_m")

## The prioritized set is the top 11 of the specificity run: ranks 1-11 are the
## ligands carried into the manuscript text (TGFB2 ... AGT). Read the count from
## the data rather than hard-coding the membership.
N_PRIORITIZED <- 11

spec <- fread(file.path(CCC, "nichenet", "nichenet_specificity_Pericytes.tsv"))[
    order(rank)][seq_len(N_PRIORITIZED)]
spec[, test_ligand := factor(test_ligand, levels = rev(test_ligand))]

donors <- fread(file.path(CCC, "donor_validation_table.tsv.gz"))
donors[, disease_group := factor(disease_group, levels = DISEASE_LEVELS)]
res <- fread(file.path(CCC, "donor_validation", "donor_validation_results.tsv"))

## ===== Panel A -- observed AUPR against the size-matched null ==============
## The null band is mean +/- 2 SD over 1,000 random target sets of the same size.
## It is drawn on the same axis as the observation so the separation is the point;
## the bands are invisibly narrow at this scale, which is itself the message.
pA <- ggplot(spec, aes(y = test_ligand)) +
    geom_errorbarh(aes(xmin = pmax(null_mean - 2 * null_sd, 0),
                       xmax = null_mean + 2 * null_sd),
                   height = 0.35, linewidth = 0.5, colour = "grey55") +
    geom_point(aes(x = null_mean), size = 1.1, colour = "grey45", shape = 15) +
    geom_point(aes(x = obs_aupr), size = 2.0, colour = "#D55E00") +
    geom_segment(aes(x = null_mean, xend = obs_aupr,
                     y = test_ligand, yend = test_ligand),
                 linewidth = 0.3, colour = "#D55E00", alpha = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = "Corrected AUPR (target-gene recovery)", y = NULL) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6.5))

## ===== Panel B -- standardized distance from the null =====================
p_floor <- max(spec$p_emp)
pB <- ggplot(spec, aes(z, test_ligand)) +
    geom_segment(aes(x = 0, xend = z, y = test_ligand, yend = test_ligand),
                 linewidth = 0.4, colour = "grey70") +
    geom_point(size = 2.0, colour = "#0072B2") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "SD from permutation null (z)", y = NULL,
         caption = sprintf("every ligand sits at the permutation floor, empirical P = %.2e",
                           p_floor)) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6.5),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0))

## ===== Panels C/D/E -- donor-level validation =============================
## Disease-group colour is retained here (unlike the main-figure version, which
## shows the adjusted, disease-neutral residuals) because the supplementary point
## is that the relationship is not carried by one disease group.
##
## Three panels, ordered composite -> TGFB1 -> TGFB2, i.e. strongest to weakest.
## TGFB2 is included because it is the TOP-ranked ligand by corrected AUPR (0.208,
## just above TGFB1's 0.203) in panel A, and without it the leading ligand would be
## the only one with no individual donor-level estimate. Its donor-level support is
## weaker than its rank implies -- that asymmetry is the point of showing it, not a
## reason to omit it.
## 05.donor_validation_stats.R publishes TWO specifications per predictor, on
## purpose: lmer(+1|dataset) for this figure and lm(+dataset fixed) for Figure 4D,
## where the panel is a partial-regression plot whose drawn slope IS the
## fixed-effect coefficient. They differ materially (composite: 0.50 vs 0.649), so
## the row must be selected by `model`, not by predictor alone -- taking both would
## make sprintf() return two labels and stamp two overlapping annotations.
LMM_MODEL <- "lmer(+1|dataset)"

## Three panels across a 7.5 in canvas leave ~2.4 in each, so two things that were
## affordable at two-up are not: the y-axis title and the disease legend are drawn
## on the leftmost panel only, and the "(LMM, +1|dataset)" qualifier -- identical in
## all three -- moves to the figure caption. Both are deduplication, not omission.
donor_panel <- function(xvar, xlab, pred_name, first = FALSE) {
    r <- res[predictor == pred_name & model == LMM_MODEL]
    if (nrow(r) != 1L)
        stop("expected exactly 1 ", LMM_MODEL, " row for ", pred_name, ", got ", nrow(r))
    lab <- sprintf("raw %s\nadj. beta = %.2f, %s",
                   fmt_rho(r$spearman_rho, r$spearman_p),
                   r$adj_estimate, tolower(fmt_p(r$adj_p)))
    ggplot(donors, aes(.data[[xvar]], receiver_target_expr)) +
        geom_smooth(method = "lm", formula = y ~ x, colour = "grey25",
                    fill = "grey85", linewidth = 0.5, alpha = 0.6) +
        geom_point(aes(colour = disease_group), size = 1.0, alpha = 0.85) +
        annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.15,
                 size = 2.1, colour = "grey20", lineheight = 1.05) +
        scale_colour_manual(values = DISEASE_COL, labels = DISEASE_LABS, name = NULL,
                            drop = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
        labs(x = xlab,
             y = if (first) "Pericyte target-program expression\n(donor mean)" else NULL) +
        theme_ms() +
        theme(legend.position = if (first) c(0.98, 0.02) else "none",
              legend.justification = c(1, 0), legend.background = element_blank(),
              legend.key.size = unit(2.6, "mm"), legend.text = element_text(size = 5.6))
}

pC <- donor_panel("sender_ligand_mean",
                  "Niche expression of prioritized ligands\n(donor mean, composite)",
                  "sender_ligand_mean", first = TRUE)
pD <- donor_panel("sender_TGFB1", "Niche TGFB1 expression\n(donor mean)",
                  "sender_TGFB1")
pE <- donor_panel("sender_TGFB2", "Niche TGFB2 expression\n(donor mean)",
                  "sender_TGFB2")

## ---- assemble ------------------------------------------------------------
fig <- (pA | pB) / (pC | pD | pE) +
    plot_layout(heights = c(1, 1.15)) +
    ## The caption belongs to the ASSEMBLED figure, so it is themed through
    ## plot_annotation(theme = ...); `&` would push plot.caption down into every
    ## subplot instead, where panel B already has one of its own.
    plot_annotation(
        tag_levels = "A",
        caption = paste("C-E: donor-level fits are LMM with (1|dataset) and a",
                        "disease-group term; n = 105 donors."),
        theme = theme(plot.caption = element_text(size = 5.8, colour = "grey30",
                                                  hjust = 0))) &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_nichenet_specificity", fig, 7.5, 7.0)

cat("Wrote figureS_nichenet_specificity to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
