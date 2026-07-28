## Supplementary figure S3: stability and annotation of human lung pericyte states.
##
## The six pericyte states used throughout this study are data-driven Leiden
## clusters on the study-integrated embedding, selected by a bootstrap stability
## sweep -- they are NOT defined by marker-score argmax. This figure is the audit
## trail for that claim, in the order a sceptical reader would ask for it:
##
##   A  the sweep itself -- median bootstrap Jaccard across neighbours x
##      resolution, with the selected solution ringed
##   B  the selected solution broken out per cluster, so the median in A is not
##      hiding one unstable cluster
##   C  what the clusters express: top row = data-driven Wilcoxon markers per
##      cluster, bottom row = the curated program panels. Markers CHARACTERISE the
##      clusters here; they played no part in defining them
##   D  relative program enrichment per cluster, and the resulting assignment
##
## Panels A-B are about the clustering; C-D are about what the clusters turned out
## to be. The direction of that dependency is the point of the figure.
##
## There is deliberately no program-score UMAP panel: figureS_pericyte_layer (S1)
## panel A already shows exactly those six overlays on this same embedding, and a
## panel duplicated across two supplements is worse than a cross-reference.
## Panel C's curated block carries the cluster-to-program link instead.
##
## Clusters are shown in NUMERIC order P0-P5 (descending size) rather than the
## program-grouped order used in S1/S12, because this figure is about the
## clustering solution itself and the program assignment is an OUTPUT (panel E),
## not an ordering key. Cluster hues are the shared ones, so identity still
## carries across figures.
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

PS <- P("pericyte_states", "_m")

CLUST_ORDER <- paste0("P", 0:5)
## Shared cluster hues (figureS_alluvial, figureS_pericyte_layer): blues =
## vascular-stabilizing, oranges = basement-membrane, pink = activated/migratory.
CLUST_COL <- c(P0 = "#08519C", P2 = "#6BAED6", P1 = "#D94801",
               P3 = "#FD8D3C", P5 = "#FDD0A2", P4 = "#CC79A7")

PROG_LEVELS <- c("vascular_stabilizing", "synthetic_contractile",
                 "activated_migratory", "inflammatory", "fibroblast_like",
                 "basement_membrane")
PROG_LABS <- c(vascular_stabilizing = "Vascular-stabilizing",
               synthetic_contractile = "Synthetic/contractile",
               activated_migratory = "Activated/migratory",
               inflammatory = "Inflammatory",
               fibroblast_like = "Fibroblast-like",
               basement_membrane = "Basement-membrane")

pclust <- function(x, levels = CLUST_ORDER)
    factor(paste0("P", x), levels = levels)

## ===== Panel A -- the stability sweep =====================================
grid   <- fread(file.path(PS, "stability", "cluster_stability_grid.tsv"))
chosen <- fread(file.path(PS, "stability", "chosen_cluster_solution.tsv"))
STAB_THRESHOLD <- 0.6   # --stability-threshold default of 00.state_discovery.py

grid[, nn := factor(n_neighbors, labels = paste0(sort(unique(n_neighbors)),
                                                 " neighbours"))]
grid[, is_chosen := cluster_key == chosen$cluster_key]

## The y axis is expanded down to the pass gate on purpose: every setting clears
## 0.6 by a wide margin, so the selected solution was chosen among stable
## alternatives rather than rescued by the threshold.
pA <- ggplot(grid, aes(resolution, median_jaccard, colour = nn, group = nn)) +
    geom_hline(yintercept = STAB_THRESHOLD, linetype = "dashed",
               colour = "grey55", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.7) +
    geom_point(data = grid[is_chosen == TRUE], shape = 21, size = 4.4,
               stroke = 0.7, fill = NA, colour = "black", show.legend = FALSE) +
    ## The two series nearly touch at resolution 0.3, so the k labels are pushed
    ## apart by series rather than both sitting above their point.
    geom_text(aes(label = paste0("k=", n_clusters),
                  vjust = ifelse(n_neighbors == 30, -1.5, 2.1)),
              size = 2.0, show.legend = FALSE) +
    annotate("text", x = -Inf, y = STAB_THRESHOLD, label = " pass gate",
             hjust = 0, vjust = -0.5, size = 2.0, colour = "grey40") +
    scale_colour_manual(values = c(OKABE[6], OKABE[4]), name = NULL) +
    scale_x_continuous(breaks = sort(unique(grid$resolution))) +
    scale_y_continuous(limits = c(0.55, 1.03),
                       expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Leiden resolution", y = "Median bootstrap Jaccard",
         caption = sprintf("selected: %d neighbours, resolution %s, k = %d",
                           chosen$n_neighbors, chosen$resolution,
                           chosen$n_clusters)) +
    theme_ms() +
    ## Bottom-RIGHT: the gate annotation sits bottom-left, and the data occupy the
    ## top of the panel, so this is the only corner free of both.
    theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
          legend.background = element_blank(), legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0))

## ===== Panel B -- per-cluster stability of the selected solution ==========
## 00.state_discovery.py flattens the per-cluster Jaccards into one median per
## setting, so a median of 0.97 is compatible with one fragile cluster. Its
## --stability-only mode (step_0c.sh) re-runs the bootstrap at the stored solution
## and persists the per-cluster summary; it verifies that the recomputed labels
## reproduce `pericyte_state` for every cell before writing. Those values are a
## FRESH bootstrap rather than a replay of the sweep -- the sweep consumes one RNG
## stream across all eight settings, so a single-setting run necessarily draws a
## different resample sequence. Hence "fresh median" in the caption, not "the"
## median: panel A's 0.966 and this panel's median are two draws of the same
## quantity, agreeing to within the 0.05 tolerance the script enforces.
jac    <- fread(file.path(PS, "stability", "cluster_bootstrap_jaccard.tsv"))
counts <- fread(file.path(PS, "state_counts.tsv"))
setnames(counts, c("pericyte_state", "count"))

## Rows for the non-selected settings carry an empty pericyte_state; --stability-only
## writes only the selected one, but guard anyway so a full re-run cannot leak
## eight settings' worth of clusters into this panel.
jac <- jac[!is.na(pericyte_state)]
jac[, cluster := pclust(pericyte_state)]
jac <- merge(jac, counts[, .(cluster = pclust(pericyte_state), n = count)],
             by = "cluster")[order(cluster)]

## Cluster size on the axis, because the smallest clusters are the ones a
## stability sweep is there to protect against.
clab <- setNames(sprintf("%s\n(%s)", jac$cluster,
                         format(jac$n, big.mark = ",", trim = TRUE)),
                 as.character(jac$cluster))
fresh_median <- median(jac$jaccard_median)

pB <- ggplot(jac, aes(cluster, jaccard_median)) +
    geom_hline(yintercept = fresh_median, linetype = "dotted",
               colour = "grey40", linewidth = 0.35) +
    geom_linerange(aes(ymin = jaccard_lo, ymax = jaccard_hi), colour = "grey45",
                   linewidth = 0.4) +
    geom_point(aes(y = jaccard_min), shape = 4, size = 1.1, colour = "grey45",
               stroke = 0.5) +
    geom_point(aes(fill = cluster), shape = 21, size = 2.6, stroke = 0.35,
               colour = "grey20") +
    geom_text(aes(y = 1.055, label = sprintf("%.2f", jaccard_median)),
              size = 2.0, colour = "grey20") +
    scale_fill_manual(values = CLUST_COL, guide = "none") +
    scale_x_discrete(labels = clab) +
    scale_y_continuous(limits = c(NA, 1.08)) +
    labs(x = NULL, y = "Bootstrap Jaccard",
         caption = sprintf(paste("point = median over %d bootstraps at 80%% resampling,",
                                 "bar = 2.5-97.5%% range, cross = worst replicate.",
                                 "\ndotted line = fresh median across clusters (%.3f)"),
                           max(jac$n_boot), fresh_median)) +
    theme_ms() +
    theme(axis.text.x = element_text(size = 6.2, lineheight = 0.95),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0,
                                      lineheight = 1.1))

## ===== Panel C -- what the clusters express ===============================
## Two rows of six equally sized blocks: data-driven Wilcoxon markers per cluster
## on top, curated program panels underneath. Blocks are equal-sized (four genes
## each, set by 00b.annotation_support.py) so facet_wrap's uniform column widths
## leave every gene label legible -- the reason for facet_wrap here rather than
## facet_grid(space = "free_x").
##
## Fill is the z of mean expression ACROSS clusters within a gene, so a row can be
## read for specificity; dot size is the raw detection rate, which fill cannot show.
dot <- fread(file.path(PS, "annotations", "state_marker_dotplot.tsv"))
dot[, cluster := pclust(pericyte_state, rev(CLUST_ORDER))]
dot[, block_lab := fifelse(block_type == "marker",
                           paste0("P", block), PROG_LABS[block])]
BLOCK_ORDER <- c(paste0("P", 0:5), unname(PROG_LABS[PROG_LEVELS]))
dot[, block_lab := factor(block_lab, levels = BLOCK_ORDER)]
dot[, gene := factor(gene, levels = unique(gene[order(block_lab, block_rank)]))]
dot[, z := as.numeric(scale(mean_expr)), by = gene]

pC <- ggplot(dot, aes(gene, cluster)) +
    geom_point(aes(size = frac_expr, fill = z), shape = 21, colour = "grey30",
               stroke = 0.2) +
    facet_wrap(~ block_lab, nrow = 2, scales = "free_x") +
    scale_size_area(max_size = 3.2, name = "Fraction\nexpressing",
                    breaks = c(0.2, 0.5, 0.8)) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Mean expression\n(z across clusters)") +
    labs(x = NULL, y = NULL,
         caption = paste("top: Wilcoxon markers per cluster (ribosomal and mitochondrial genes excluded from display).",
                         "bottom: curated program panels, four most cluster-discriminating genes each.")) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.4,
                                     face = "italic"),
          axis.text.y = element_text(size = 6.5),
          panel.grid.major = element_line(linewidth = 0.2),
          strip.text = element_text(size = 6.0, face = "bold"),
          panel.spacing.x = unit(1.5, "pt"), panel.spacing.y = unit(4, "pt"),
          legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0,
                                      lineheight = 1.1))

## ===== Panel D -- relative program enrichment and the assignment ==========
pm <- fread(file.path(PS, "annotations", "state_program_map.tsv"))
rel <- melt(pm, id.vars = c("pericyte_state", "state_program"),
            measure.vars = paste0(PROG_LEVELS, "_relenrich"),
            variable.name = "program", value.name = "relenrich")
rel[, program := sub("_relenrich$", "", as.character(program))]
rel[, assigned := program == state_program]
rel[, cluster := pclust(pericyte_state, rev(CLUST_ORDER))]
rel[, program := factor(PROG_LABS[program], levels = unname(PROG_LABS[PROG_LEVELS]))]

pD <- ggplot(rel, aes(program, cluster)) +
    geom_tile(aes(fill = relenrich), colour = "white", linewidth = 0.5) +
    geom_tile(data = rel[assigned == TRUE], fill = NA, colour = "black",
              linewidth = 0.7) +
    geom_text(aes(label = sprintf("%.2f", relenrich)), size = 2.0,
              colour = "grey15") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Relative\nenrichment (z)") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         caption = "black outline = dominant program assigned to that cluster") +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 30, hjust = 0, size = 6.2),
          axis.text.y = element_text(size = 6.5),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "right", legend.key.size = unit(3.2, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.6),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0))

## ---- assemble ------------------------------------------------------------
## Panel C gets the most height: it is two stacked rows of blocks and its 45-degree
## gene labels need vertical room, whereas D is a 6 x 6 tile grid.
fig <- (pA | pB) / pC / pD +
    plot_layout(heights = c(1.0, 1.55, 0.95)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_state_annotation", fig, 7.5, 10.0)

cat("Wrote figureS_state_annotation to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
