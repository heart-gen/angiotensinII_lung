## Supplementary figure S4: the AGTR1-undetected pericytes are a detection
## phenomenon, not a second population.
##
## AGTR1 is detected in a minority of pericytes, which invites the reading that the
## rest are a distinct AGTR1-negative subpopulation living somewhere else in the
## lung. This figure closes that off from two directions:
##
##   A  the zero count itself is unremarkable -- AGTR1 sits inside the distribution
##      of zero rates for genes of the same pooled abundance in the same cells
##   B  the two groups occupy the same airspace-affinity distribution
##   C  what survives adjustment is a ~0.11 SD shift that vanishes at donor level
##   D  the zero-count cells are not AGTR1-low (their denoised AGTR1 spans two
##      orders of magnitude and its median EXCEEDS that of the detected cells), and
##      within them denoised AGTR1 carries no airspace signal
##   E  the zeros concentrate in shallow libraries, tracking the matched-gene
##      expectation across depth almost exactly
##
## Panels A and E are complementary, not redundant: A asks whether the NUMBER of
## zeros is unusual, E asks whether the zeros land where dropout predicts they
## should. A gene with a normal zero count but depth-independent zeros would pass A
## and fail E.
##
## WHAT THIS FIGURE DOES NOT CLAIM. The per-cell effect in C is significant, and it
## is significant on the scVI-denoised continuous lens too (-0.087 SD per SD,
## P = 6e-09), so there IS a real graded association between AGTR1 level and
## distance from the airspace -- as a perivascular receptor should have. The claim
## is narrower: detectability is not the boundary of a population. The zero count
## is what dropout predicts (A), the two groups' distributions coincide (B), the
## contrast does not survive donor aggregation (C, P = 0.79), and the zeros are
## internally graded with no spatial structure of their own (D).
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

DROP <- P("localization", "airspace_analysis", "_m", "agtr1_dropout")

## Grey for the undetected cells throughout, so the eye reads "the grey group is
## not somewhere else" in B and D without consulting a legend twice.
DET_LEVELS <- c("AGTR1-undetected", "AGTR1-detected")
DET_COL    <- c(`AGTR1-undetected` = "#999999", `AGTR1-detected` = "#D55E00")

summ  <- fread(file.path(DROP, "dropout_model_summary.tsv"))
null  <- fread(file.path(DROP, "matched_gene_null.tsv"))
depth <- fread(file.path(DROP, "detection_by_depth.tsv"))
ladd  <- fread(file.path(DROP, "airspace_effect_ladder.tsv"))
cells <- fread(file.path(DROP, "agtr1_cells.tsv.gz"))

cells[, detect_lab := factor(fifelse(AGTR1_detect == 1L, "AGTR1-detected",
                                     "AGTR1-undetected"), levels = DET_LEVELS)]

## ===== Panel A -- observed zero count against the matched-gene null ========
## The null is the zero rate of the 200 genes matched to AGTR1 on pooled abundance.
## Drawn as the distribution rather than as a single expected value because the
## claim is "AGTR1 is inside this", and a bar chart of two numbers cannot show that.
##
## Observed (0.626) and matched median (0.639) land within one bin of each other,
## so they are distinguished by a mapped legend rather than by in-plot rotated
## labels, which would collide.
refs <- data.table(
    kind = factor(c("AGTR1 (observed)", "matched-gene median"),
                  levels = c("AGTR1 (observed)", "matched-gene median")),
    x = c(summ$observed_undetected_frac, summ$matched_expected_undetected_frac))
REF_COL <- c(`AGTR1 (observed)` = DET_COL[["AGTR1-detected"]],
             `matched-gene median` = "grey30")

pA <- ggplot(null, aes(frac_undetected)) +
    geom_histogram(bins = 40, fill = "grey85", colour = "grey60", linewidth = 0.2) +
    geom_vline(data = refs, aes(xintercept = x, colour = kind, linetype = kind),
               linewidth = 0.7) +
    ## Right-aligned: the reference rules sit at ~63%, near the left of the axis,
    ## and a left-aligned block is struck through by both of them.
    annotate("text", x = Inf, y = Inf, hjust = 1.04, vjust = 1.25, size = 2.1,
             colour = "grey20", lineheight = 1.15,
             label = sprintf(
                 "observed %s / %s (%.1f%%)\nexpected %s [%s-%s] (%.1f%%)\nz = %.2f, empirical P = %.2f",
                 format(summ$observed_undetected_n, big.mark = ","),
                 format(summ$n_cells, big.mark = ","),
                 100 * summ$observed_undetected_frac,
                 format(round(summ$matched_expected_undetected_n), big.mark = ","),
                 format(round(summ$matched_lo_n), big.mark = ","),
                 format(round(summ$matched_hi_n), big.mark = ","),
                 100 * summ$matched_expected_undetected_frac,
                 summ$z, summ$p_empirical)) +
    scale_colour_manual(values = REF_COL, name = NULL) +
    scale_linetype_manual(values = c(`AGTR1 (observed)` = "solid",
                                     `matched-gene median` = "dashed"), name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
    labs(x = sprintf("Pericytes with zero counts\n(%d genes matched to AGTR1 on abundance)",
                     summ$n_matched_genes),
         y = "Matched genes") +
    theme_ms() +
    ## The annotation block already occupies the top of the panel, so the key goes
    ## underneath rather than inside; at this panel width the two collide.
    theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0),
          legend.key = element_blank(), legend.key.size = unit(3.2, "mm"),
          legend.text = element_text(size = 5.8))

## ===== Panel B -- airspace affinity by detectability ======================
## Violin + boxplot rather than two densities: the distributions are near-identical,
## and overlaid densities at that similarity read as a single curve.
nlab <- cells[, .(n = .N, y = max(airspace_score)), by = detect_lab]
pB <- ggplot(cells, aes(detect_lab, airspace_score, fill = detect_lab)) +
    geom_violin(scale = "width", width = 0.85, colour = NA, alpha = 0.55) +
    geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3,
                 fill = "white", colour = "grey25") +
    geom_text(data = nlab, inherit.aes = FALSE,
              aes(x = detect_lab, y = y,
                  label = paste0("n = ", format(n, big.mark = ","))),
              vjust = -0.6, size = 2.1, colour = "grey20") +
    scale_fill_manual(values = DET_COL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = NULL, y = "Airspace affinity\n(mean cosine similarity to alveolar centroids)") +
    theme_ms() +
    theme(axis.text.x = element_text(size = 7))

## ===== Panel C -- the specification ladder ================================
## Every estimate is in SD of airspace affinity, so the rows are comparable. The
## donor-level row has a 0-1 predictor and therefore the largest possible contrast;
## its interval straddling zero is the load-bearing result of the panel.
ladd[, label := factor(label, levels = rev(label))]
ladd[, sig := pval < 0.05]
## fmt_p() branches on a scalar, so it cannot be called inside aes() over 5 rows.
ladd[, p_lab := vapply(pval, fmt_p, character(1))]
pC <- ggplot(ladd, aes(estimate, label)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_linerange(aes(xmin = ci_lower, xmax = ci_upper), linewidth = 0.5,
                   colour = "grey40") +
    geom_point(aes(fill = sig), shape = 21, size = 2.4, colour = "grey20",
               stroke = 0.3) +
    geom_text(aes(x = ci_upper, label = paste0("  ", p_lab)),
              hjust = 0, size = 2.0, colour = "grey30") +
    scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "white"),
                      guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.30))) +
    labs(x = "Effect on airspace affinity (SD)", y = NULL,
         caption = sprintf(
             "per-cell rows: LMM with (1|donor) + dataset variance component, n = %s cells / %d donors; donor row: OLS, n = %d donors",
             format(ladd$n_cells[1], big.mark = ","), ladd$n_donors[1],
             ladd[spec == "donor-level OLS", n_donors])) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6.8),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0))

## ===== Panel D -- inside the zeros ========================================
## Restricted to cells with NO observed AGTR1 transcript. Two things are visible at
## once. Horizontally: denoised AGTR1 spans two orders of magnitude among these
## cells and its median sits ABOVE the median of the detected cells (dashed rule),
## so "undetected" is not a low-AGTR1 label. Vertically: the fit is flat, so within
## the undetected set the denoised level carries no airspace signal -- there is no
## hidden subpopulation being masked by the zeros.
##
## x is log10 because the denoised values are strongly right-skewed (median 0.77,
## max 10.8); on a linear axis 95% of the cells pile into the left tenth.
zero    <- cells[AGTR1_detect == 0L]
det_med <- cells[AGTR1_detect == 1L, median(AGTR1_scvi)]
rho     <- suppressWarnings(cor.test(zero$AGTR1_scvi, zero$airspace_score,
                                     method = "spearman"))
pD <- ggplot(zero, aes(AGTR1_scvi, airspace_score)) +
    geom_point(size = 0.2, alpha = 0.12, colour = "grey40") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "#0072B2",
                fill = "#0072B2", alpha = 0.25, linewidth = 0.6) +
    geom_vline(xintercept = det_med, linetype = "dashed",
               colour = DET_COL[["AGTR1-detected"]], linewidth = 0.5) +
    annotate("text", x = det_med, y = -Inf, hjust = 0.5, vjust = -0.8,
             size = 1.9, colour = DET_COL[["AGTR1-detected"]],
             label = "median of detected cells") +
    ## Left-aligned: the detected-cell median rule lands mid-axis on the log scale
    ## and cuts through a right-aligned block. Anchored to the data minimum rather
    ## than -Inf, which scale_x_log10() turns into NaN and silently drops.
    annotate("text", x = min(zero$AGTR1_scvi), y = Inf, hjust = 0, vjust = 1.3,
             size = 2.1,
             colour = "grey20", lineheight = 1.15,
             label = sprintf("%s\nn = %s zero-count cells",
                             fmt_rho(unname(rho$estimate), rho$p.value),
                             format(nrow(zero), big.mark = ","))) +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = "scVI-denoised AGTR1, log scale\n(pericyte-trained model)",
         y = "Airspace affinity") +
    theme_ms()

## ===== Panel E -- detection against sequencing depth ======================
## Ribbon = the 2.5-97.5 percentile of the matched genes' detection rate in the
## SAME cells, so the comparison is within-bin and needs no model.
## The two bands mean different things and must not be read as comparable error
## bars: the AGTR1 interval is sampling error on one gene (Wilson, binomial), the
## grey band is BETWEEN-GENE spread across the 200 matched genes. That is why the
## grey is far wider, and why the legend keys name what each one is.
SER_COL <- c(`AGTR1 (95% Wilson CI)` = unname(DET_COL[["AGTR1-detected"]]),
             `matched genes (median)` = "grey35")
pE <- ggplot(depth, aes(depth_median)) +
    geom_ribbon(aes(ymin = matched_lo, ymax = matched_hi,
                    fill = "matched genes (2.5-97.5 pct)"), alpha = 0.55) +
    geom_line(aes(y = matched_rate, colour = "matched genes (median)"),
              linetype = "dashed", linewidth = 0.45) +
    geom_linerange(aes(ymin = obs_lo, ymax = obs_hi,
                       colour = "AGTR1 (95% Wilson CI)"), linewidth = 0.4) +
    geom_line(aes(y = obs_rate, colour = "AGTR1 (95% Wilson CI)"),
              linewidth = 0.6) +
    geom_point(aes(y = obs_rate, colour = "AGTR1 (95% Wilson CI)"), size = 1.5) +
    scale_colour_manual(values = SER_COL, name = NULL) +
    scale_fill_manual(values = c(`matched genes (2.5-97.5 pct)` = "grey80"),
                      name = NULL) +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.05, 0.28))) +
    labs(x = "Library size (UMIs per cell, depth-decile median)",
         y = "Cells with a detected transcript") +
    theme_ms() +
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.background = element_blank(), legend.key = element_blank(),
          legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 5.8),
          legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 0, 0, 0))

## ---- assemble ------------------------------------------------------------
## C is a five-row forest: short, and it needs the full width for its labels and
## p-value annotations, so it gets its own band rather than sharing a row.
##
## wrap_elements() is load-bearing, not decoration. patchwork aligns panel regions
## down the grid, so C's long y-axis labels ("binary detectability, depth-adjusted")
## otherwise reserve that width in every row, pushing A and D into the right third
## of their cells and leaving a blank left margin. Wrapping C makes it opaque to the
## alignment machinery; A/B and D/E still align with each other.
fig <- (pA | pB) / wrap_elements(full = pC) / (pD | pE) +
    plot_layout(heights = c(1.05, 0.80, 1.00)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_agtr1_dropout", fig, 7.5, 8.6)

cat("Wrote figureS_agtr1_dropout to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
