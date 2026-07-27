## Supplementary figure S8: is the prioritized-ligand result specific, and does it
## hold at the level of individual donors?
##
## NicheNet ranks ligands by how well their prior regulatory potential recovers a
## curated target gene set. Two fair objections follow: (i) any gene set of that
## size might score well (the ranking would then be an artifact of the prior network,
## not of pericyte biology), and (ii) the whole result is computed on pooled
## pseudobulk, so it might not hold donor by donor. Panels A-B answer the first with
## a size-matched gene-set permutation null; C-D answer the second with the
## donor-level sender/receiver relationship.
##
##   A  observed corrected AUPR per prioritized ligand vs its permutation null
##   B  standardized distance from that null, with the empirical p floor
##   C  donor-level composite prioritized-ligand signal vs pericyte target response
##   D  the same for TGFB1 alone, raw and dataset-adjusted
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

## ===== Panels C/D -- donor-level validation ===============================
## Disease-group colour is retained here (unlike the main-figure version, which
## shows the adjusted, disease-neutral residuals) because the supplementary point
## is that the relationship is not carried by one disease group.
donor_panel <- function(xvar, xlab, pred_name) {
    r <- res[predictor == pred_name]
    lab <- sprintf("raw %s\nadj. beta = %.2f, %s (LMM, +1|dataset)",
                   fmt_rho(r$spearman_rho, r$spearman_p),
                   r$adj_estimate, tolower(fmt_p(r$adj_p)))
    ggplot(donors, aes(.data[[xvar]], receiver_target_expr)) +
        geom_smooth(method = "lm", formula = y ~ x, colour = "grey25",
                    fill = "grey85", linewidth = 0.5, alpha = 0.6) +
        geom_point(aes(colour = disease_group), size = 1.1, alpha = 0.85) +
        annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.15,
                 size = 2.2, colour = "grey20", lineheight = 1.05) +
        scale_colour_manual(values = DISEASE_COL, labels = DISEASE_LABS, name = NULL,
                            drop = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
        labs(x = xlab, y = "Pericyte target-program expression\n(donor mean)") +
        theme_ms() +
        theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
              legend.background = element_blank(),
              legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6))
}

pC <- donor_panel("sender_ligand_mean",
                  "Niche expression of prioritized ligands\n(donor mean, composite)",
                  "sender_ligand_mean")
pD <- donor_panel("sender_TGFB1", "Niche TGFB1 expression\n(donor mean)",
                  "sender_TGFB1")

## ---- assemble ------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
    plot_layout(heights = c(1, 1.15)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_nichenet_specificity", fig, 7.5, 7.0)

cat("Wrote figureS_nichenet_specificity to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
