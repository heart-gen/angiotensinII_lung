## Supplementary figure S6: an unsupervised factorisation of the same cells
## recovers the curated pericyte programs, at a rank the data chose.
##
## The curated 6-program NVU model is marker-panel driven, so it can only find what
## its panels encode. CoGAPS is given no panels: it factorises the HVG matrix into
## nP patterns de novo. The figure asks three separate questions and keeps them
## separate, because they can fail independently:
##
##   A  what rank does the data support?  Cross-seed reproducibility across
##      nP = 4-10, with the band over patterns. nP=8 (main) maximises the WEAKEST
##      pattern's reproducibility; nP=9 (sensitivity) is the largest rank that still
##      clears the 0.80 threshold.
##   B  do the de-novo patterns carry the curated programs?  Cell-level Spearman of
##      pattern weight against each curated score at nP=8.
##   C  is that agreement also visible in the GENES, not just the scores?  Top-50
##      PatternMarker overlap with the six curated panels (hypergeometric, BH).
##   D  is the nP=8 solution itself stable?  Per-pattern matched r across three
##      non-canonical seeds.
##   E  does any of this depend on picking 8 over 9?  Every nP=8 pattern matches an
##      nP=9 pattern, and matched pairs keep the same curated-program profile.
##
## B and C are deliberately both pattern x program grids: they are independent
## lines of evidence (continuous cell-level concordance vs discrete marker-gene
## membership) and a pattern can pass one and fail the other. Pattern_6 does exactly
## that -- fibroblast-like by score (rho 0.38) with zero marker overlap anywhere.
##
## WHAT THIS FIGURE DOES NOT CLAIM. It does not claim the CoGAPS patterns ARE the
## curated programs one-for-one. The mapping is many-to-one (two patterns lead on
## inflammatory, two on fibroblast-like, two on basement membrane), and Pattern_2
## and Pattern_8 correlate negatively with five of the six panels -- they behave as
## low-activity/baseline factors and their "basement_membrane" argmax label rests on
## a single modest positive rho. The claim is that the curated program AXES fall out
## of an unsupervised factorisation of the same cells, and that this does not depend
## on the rank chosen.
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

CG   <- P("pericyte_cogaps", "_m")
NP   <- 8L      # main rank
NPS  <- 9L      # sensitivity rank
THRESH <- 0.80  # stability threshold used by 02.select_rank.R

## Curated panels ordered stabilising -> synthetic -> matrix -> injury, so both
## grid panels read left-to-right along the same biological axis.
PROGRAMS <- c("vascular_stabilizing", "synthetic_contractile", "basement_membrane",
              "fibroblast_like", "activated_migratory", "inflammatory")
PROG_LABS <- c(vascular_stabilizing = "vascular\nstabilizing",
               synthetic_contractile = "synthetic\ncontractile",
               basement_membrane = "basement\nmembrane",
               fibroblast_like = "fibroblast-\nlike",
               activated_migratory = "activated\nmigratory",
               inflammatory = "inflammatory")
PROG_ONE  <- c(vascular_stabilizing = "vascular stabilizing",
               synthetic_contractile = "synthetic contractile",
               basement_membrane = "basement membrane",
               fibroblast_like = "fibroblast-like",
               activated_migratory = "activated migratory",
               inflammatory = "inflammatory")

## "Pattern_3" -> "P3". Used in every panel so a pattern can be tracked across them.
short_pat <- function(x) sub("^Pattern_", "P", x)
pat_factor <- function(x, n) factor(short_pat(x), levels = paste0("P", seq_len(n)))
## Passed to every scale_y_discrete explicitly. Relying on the factor's own levels
## is not enough: patchwork/ggplot resolves a discrete scale from the union of the
## layers, and panel C's zero-overlap layer does not contain every pattern, which
## silently reorders the axis.
PAT_LEVELS <- paste0("P", seq_len(NP))

sel   <- fread(file.path(CG, "cogaps_nP_selection.tsv"))
stab  <- fread(file.path(CG, "cogaps_seed_stability.tsv"))
sco   <- fread(file.path(CG, sprintf("validation_np%d", NP), "pattern_score_spearman.tsv.gz"))
ovl   <- fread(file.path(CG, sprintf("validation_np%d", NP), "pattern_panel_overlap.tsv.gz"))
conc  <- fread(file.path(CG, sprintf("cogaps_np%d_vs_np%d_concordance.tsv", NP, NPS)))

## ===== Panel A -- rank selection ==========================================
## The line is the mean over patterns of each pattern's mean matched r across the
## three non-canonical seeds; the band is the min-max over patterns. The band is
## the informative part: mean_r stays high everywhere (>= 0.92) because most
## patterns are trivially reproducible, and it is the FLOOR that separates the
## ranks. nP = 5, 6 and 10 fail on a single collapsing pattern.
##
## na.rm = TRUE reproduces 02.select_rank.R exactly, so the line and the floor of
## the band are the mean_r / min_r columns of cogaps_nP_selection.tsv rather than a
## second, differently-aggregated version of them. An NA arises when a pattern
## collapses in one seed and has no correlate to match against; those are counted
## and reported in the caption instead of being silently averaged away.
pat_r <- stab[, .(r = mean(r, na.rm = TRUE)), by = .(np, ref_pattern)]
band  <- pat_r[, .(mean_r = mean(r), lo = min(r), hi = max(r)), by = np][order(np)]
stopifnot(nrow(band) == nrow(sel), all(is.finite(band$mean_r)),
          isTRUE(all.equal(band$mean_r, sel[order(np), mean_r])),
          isTRUE(all.equal(band$lo, sel[order(np), min_r])))
marks <- band[np %in% c(NP, NPS)]
marks[, role := fifelse(np == NP, "main", "sensitivity")]

## Failed matches per nP, for the caption. nP = 8 is the only rank in the upper
## half of the sweep with none.
fails <- stab[, .(n = .N, n_na = sum(!is.finite(r))), by = np][order(np)][n_na > 0]
fail_txt <- paste(sprintf("%d of %d at nP=%d", fails$n_na, fails$n, fails$np),
                  collapse = ", ")

## Distributed CoGAPS matches patterns across data subsets to build the consensus,
## and that step is free to return a different count than requested: the nP=6 fit
## came back with 7 patterns. Read off the reference rather than hard-coded, so the
## caption cannot drift from the data.
ncol_np6 <- stab[np == 6L, uniqueN(ref_pattern)]

ROLE_COL <- c(main = "#0072B2", sensitivity = "#E69F00")

pA <- ggplot(band, aes(np, mean_r)) +
    geom_hline(yintercept = THRESH, linetype = "dotted", colour = "grey45",
               linewidth = 0.4) +
    annotate("text", x = 4, y = THRESH, hjust = 0, vjust = -0.55, size = 1.9,
             colour = "grey35", label = "stability threshold (0.80)") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
    geom_line(linewidth = 0.5, colour = "grey25") +
    geom_point(size = 1.4, colour = "grey25") +
    geom_point(data = marks, aes(colour = role), size = 2.6) +
    geom_text(data = marks, aes(y = lo, colour = role, label = role),
              vjust = 1.9, size = 2.0, show.legend = FALSE) +
    ## Reconstruction error is the other standard rank signal, and it does not
    ## adjudicate: it falls monotonically with no elbow, so the rank is chosen on
    ## reproducibility. Stated rather than plotted, since a flat monotone curve on
    ## a second axis would only invite an elbow to be read into it.
    annotate("text", x = Inf, y = -Inf, hjust = 1.03, vjust = -0.55, size = 1.9,
             colour = "grey35", lineheight = 1.15,
             label = sprintf(
                 "unexplained variance falls monotonically\n%.3f (nP=4) to %.3f (nP=10): no elbow",
                 sel[np == min(np), frac_unexplained],
                 sel[np == max(np), frac_unexplained])) +
    scale_colour_manual(values = ROLE_COL, name = NULL) +
    scale_x_continuous(breaks = sel$np) +
    scale_y_continuous(limits = c(0.3, 1.02), expand = expansion(mult = c(0.10, 0.02))) +
    labs(x = "Number of CoGAPS patterns (nP)",
         y = "Cross-seed pattern reproducibility, matched |r|\n(line = mean over patterns, band = min-max)",
         caption = sprintf(
             "Failed matches (a pattern collapsing in one seed) are dropped from the mean: %s; none at nP=8.\nDistributed CoGAPS returned %d consensus patterns for the requested nP=6.",
             fail_txt, ncol_np6)) +
    theme_ms() +
    theme(axis.title.y = element_text(size = 6.4),
          legend.position = "none",
          plot.caption = element_text(size = 4.9, colour = "grey35", hjust = 0))

## ===== Panel B -- pattern x curated program, cell level ===================
## Diverging fill capped at +/-0.65 rather than +/-1: the observed range is
## -0.47 to 0.62, and a full-range scale renders every cell near-white.
sco[, program := sub("_score$", "", score)]
b <- sco[program %in% PROGRAMS & is.finite(rho)]
b[, `:=`(pattern = pat_factor(pattern, NP),
         program = factor(program, levels = PROGRAMS))]

pB <- ggplot(b, aes(program, pattern, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", rho),
                  colour = abs(rho) > 0.42), size = 1.9, show.legend = FALSE) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-0.65, 0.65),
                         oob = scales::squish, name = "Spearman rho",
                         guide = guide_colourbar(barwidth = unit(22, "mm"),
                                                 barheight = unit(2.4, "mm"))) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey25")) +
    scale_x_discrete(labels = PROG_LABS) +
    scale_y_discrete(limits = rev(PAT_LEVELS)) +
    labs(x = NULL, y = sprintf("CoGAPS pattern (nP = %d)", NP)) +
    theme_ms() +
    theme(axis.text.x = element_text(size = 5.4, lineheight = 0.9),
          panel.grid = element_blank(), panel.border = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom", legend.title = element_text(size = 6),
          legend.text = element_text(size = 5.4),
          legend.margin = margin(-4, 0, 0, 0))

## ===== Panel C -- PatternMarker overlap with the curated panels ===========
## The gene-level counterpart of B, and a much harsher test: the panels are 7-16
## genes inside a 2,038-gene HVG background, so an overlap of 3 is already strongly
## non-random. Empty cells are drawn as open crosses so a zero overlap is visibly a
## measured zero and not a dropped row.
ovl[, `:=`(pattern = pat_factor(pattern, NP),
           program = factor(program, levels = PROGRAMS),
           sig = padj < 0.05)]
psize <- unique(ovl[, .(program, panel_size)])[order(program)]
## Panel sizes go in the caption, not on the axis: appending " (11)" to the second
## line of each two-line label widens them past the panel and runs
## "activated migratory" into "inflammatory". Keeping the axis identical to B's
## also makes the two grids directly comparable row by row.
size_txt <- paste(sprintf("%s %d", PROG_ONE[as.character(psize$program)],
                          psize$panel_size), collapse = ", ")
## The HVG background the hypergeometric test was run against, read off the
## loadings rather than restated, so the caption cannot drift from the analysis.
N_HVG <- nrow(fread(P("pericyte_cogaps", "_m",
                      sprintf("feature_loadings_np%d.tsv.gz", NP)), select = 1L))

pC <- ggplot(ovl[overlap > 0], aes(program, pattern)) +
    geom_point(data = ovl[overlap == 0], shape = 4, size = 0.9, colour = "grey78") +
    geom_point(aes(size = overlap, fill = sig), shape = 21, colour = "grey25",
               stroke = 0.3) +
    scale_size_area(max_size = 4.6, breaks = c(1, 2, 3, 5), name = "genes shared") +
    scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "white"),
                      labels = c(`TRUE` = "BH < 0.05", `FALSE` = "n.s."),
                      name = NULL) +
    scale_x_discrete(labels = PROG_LABS) +
    scale_y_discrete(limits = rev(PAT_LEVELS)) +
    labs(x = NULL, y = sprintf("CoGAPS pattern (nP = %d)", NP),
         caption = sprintf("top-%d PatternMarker genes per pattern, hypergeometric against the %s-gene HVG background.\nCurated panel sizes: %s.",
                           ovl$n_markers[1], format(N_HVG, big.mark = ","),
                           size_txt)) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 2.6)),
           size = guide_legend(order = 2, nrow = 1)) +
    theme_ms() +
    theme(axis.text.x = element_text(size = 5.4, lineheight = 0.9),
          legend.position = "bottom", legend.title = element_text(size = 5.6),
          legend.text = element_text(size = 5.4),
          legend.key.size = unit(3.4, "mm"),
          legend.box = "horizontal", legend.box.spacing = unit(1.5, "mm"),
          legend.margin = margin(-2, 0, 0, 0),
          plot.caption = element_text(size = 4.9, colour = "grey35", hjust = 0))

## ===== Panel D -- cross-seed reproducibility at the main rank =============
## The nP = 8 column of panel A, opened up. One row per pattern, one point per
## non-canonical seed, matched against the canonical (seed-13) fit.
d <- stab[np == NP]
d[, `:=`(pattern = pat_factor(ref_pattern, NP), seed = factor(seed))]
dsum <- d[, .(mean_r = mean(r), lo = min(r), hi = max(r)), by = pattern]

pD <- ggplot(d, aes(r, pattern)) +
    geom_linerange(data = dsum, inherit.aes = FALSE,
                   aes(y = pattern, xmin = lo, xmax = hi), colour = "grey70",
                   linewidth = 0.5) +
    geom_point(aes(colour = seed), size = 1.5, alpha = 0.9) +
    geom_point(data = dsum, inherit.aes = FALSE, aes(x = mean_r, y = pattern),
               shape = 124, size = 2.6, colour = "grey20") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.0,
             colour = "grey25",
             label = sprintf("all %d matches |r| >= %.2f\nweakest pattern mean %.3f",
                             nrow(d), min(d$r), min(dsum$mean_r))) +
    scale_colour_manual(values = c(`1` = "#0072B2", `42` = "#009E73",
                                   `2024` = "#CC79A7"), name = "seed") +
    scale_x_continuous(limits = c(0.94, 1.001),
                       expand = expansion(mult = c(0.03, 0.02))) +
    scale_y_discrete(limits = rev(PAT_LEVELS)) +
    labs(x = "Matched |r| to the canonical fit (gene loadings)",
         y = sprintf("CoGAPS pattern (nP = %d)", NP)) +
    theme_ms() +
    ## Outside, bottom: every point sits at |r| >= 0.96, hard against the right
    ## edge, so an inside legend can only go bottom-left -- where it reads as a
    ## fourth data row. Matching C's bottom legend also keeps the two panel regions
    ## in this band the same height.
    theme(legend.position = "bottom", legend.key = element_blank(),
          legend.title = element_text(size = 5.6),
          legend.text = element_text(size = 5.4),
          legend.key.size = unit(3.0, "mm"),
          legend.margin = margin(-2, 0, 0, 0))

## ===== Panel E -- does the rank choice matter? ============================
## Two different concordances per matched pair, and they answer different things:
##   r_loadings  -- is it the same factor?  (gene-loading correlation)
##   profile_r   -- does it MEAN the same thing?  (correlation of the two 6-program
##                  rho vectors from panel B)
## A pair could score high on the first and low on the second if a pattern kept its
## genes but changed which curated programs it tracked. None do.
m <- conc[matched == TRUE]
m[, lab := sprintf("P%s -> P%s   %s",
                   sub("^Pattern_", "", pattern_np_main),
                   sub("^Pattern_", "", pattern_np_sens),
                   PROG_ONE[top_program_main])]
m[, lab := factor(lab, levels = rev(lab))]

el <- melt(m[, .(lab, `same factor (gene loadings)` = r_loadings,
                 `same meaning (program profile)` = profile_r)],
           id.vars = "lab", variable.name = "metric", value.name = "value")
E_COL <- c(`same factor (gene loadings)` = "#0072B2",
           `same meaning (program profile)` = "#D55E00")

extra <- conc[matched == FALSE]
pE <- ggplot(el, aes(value, lab, colour = metric, shape = metric)) +
    geom_linerange(data = m, inherit.aes = FALSE,
                   aes(y = lab, xmin = pmin(r_loadings, profile_r),
                       xmax = pmax(r_loadings, profile_r)),
                   colour = "grey80", linewidth = 0.4) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = E_COL, name = NULL) +
    scale_shape_manual(values = c(`same factor (gene loadings)` = 16,
                                  `same meaning (program profile)` = 17),
                       name = NULL) +
    scale_x_continuous(limits = c(0.82, 1.005),
                       expand = expansion(mult = c(0.02, 0.02))) +
    ## The summary goes in the caption, not in the panel: the rows span the full
    ## height and an in-panel block lands on the P8 row whichever corner it takes.
    labs(x = sprintf("Concordance between the nP = %d and nP = %d solutions", NP, NPS),
         y = NULL,
         caption = sprintf(
             "%d of %d nP=%d patterns matched an nP=%d pattern; the argmax curated program agrees for %d of %d. The additional pattern at nP=%d (P%s) is unmatched\nand correlates weakly with every curated panel (max rho %.2f). Row labels give the leading curated program at nP=%d.",
             nrow(m), nrow(m), NP, NPS, sum(m$program_agree), nrow(m), NPS,
             sub("^Pattern_", "", extra$pattern_np_sens), extra$top_rho_sens, NP)) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 5.8, family = "mono"),
          legend.position = "top", legend.key = element_blank(),
          legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 5.8),
          legend.margin = margin(0, 0, -4, 0),
          plot.caption = element_text(size = 4.9, colour = "grey35", hjust = 0))

## ---- assemble ------------------------------------------------------------
## E gets its own band and is wrapped: its "P4 -> P7   fibroblast like" labels are
## the widest text in the figure, and patchwork would otherwise reserve that width
## in every row, shoving A and C to the right of their cells. Wrapping makes E
## opaque to the alignment machinery; A/B and C/D still align with each other.
fig <- (pA | pB) / (pC | pD) / wrap_elements(full = pE) +
    plot_layout(heights = c(1.00, 1.00, 0.85)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_cogaps_validation", fig, 7.5, 9.0)

cat("Wrote figureS_cogaps_validation to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
