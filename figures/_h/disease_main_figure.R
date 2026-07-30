## =============================================================================
## Disease main figure (Circ Research revision) -- CONTINUOUS-INJURY VERSION,
## rebuilt 2026-07-28.
##
## The narrative the four panels have to carry, in order:
##   disease is not associated with replacement of discrete pericyte states
##     -> it is associated with GRADED injury-stromal engagement            (A)
##     -> the signal is carried by SPECIFIC continuous programs             (B)
##     -> receptor dysregulation is more evident in fibroblasts than
##        in pericytes                                                      (C)
##     -> which receives QUALIFIED support in an independent, disease-
##        enriched dataset                                                  (D)
##
## The opening clause (no discrete-state replacement) is the null established in
## pericyte_states and shown compositionally in figure_mechanism_main / the state
## composition supplement; it is the premise this figure starts from, not a panel
## here.
##
## PANEL PROVENANCE:
##   A  was C of the pre-2026-07-28 version. Same donor-endpoint plot, but now
##      the THREE-group model (Healthy / Fibrotic-ILD / Other) from 03's
##      threegroup_* outputs rather than the two-group contrast, so "graded" is
##      something the panel shows rather than asserts.
##   B  was B, expanded from one contrast to both contrasts of the same
##      three-group model, so the component decomposition is on the same footing
##      as A.
##   C  AGTR1 disease effects across stromal cell types
##      (disease_association/_h/05).
##   D  Independent GSE136831 COPD/IPF evaluation of fibroblast AGTR1
##      (disease_association/agtr1_copd_ipf/_h/01).
##
## MOVED OUT 2026-07-29: the two robustness panels -- leave-one-study-out and the
## within-study random-effects forest -- were pulled together into a single
## two-panel supplement (figureS_disease_robustness, S16). They answer two
## different objections to the SAME panel-A estimate ("one cohort drives it" and
## "it is a between-study batch artifact"), so they belong beside each other; and
## the main figure's job is the biological chain A -> B -> C -> D, which they
## interrupt. Both are still built by this script.
##
## SCALE NOTE: A and B are read off the SAME three-group model, whose programs
## are z-standardised once over the three-group donor set, so they are mutually
## comparable in SD units. The supplement's forest panel standardises over
## Healthy+Fibrotic only and is on a DIFFERENT scale -- see the supplement block
## below, where the two panels' axis titles name their standardisation sets.
##
## No in-panel titles; interpretation belongs in the caption
## (figures/mechanism/README.md).
## =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(ggplot2); library(patchwork)
})

## ROOT/P/OUT, theme_ms(), save_fig(), fmt_p()
source("../_h/_fig_common.R")

SD <- P("disease_association", "_m", "mixed_model_forest")
MD <- P("disease_association", "_m", "mean_expr")
GD <- P("disease_association", "agtr1_copd_ipf", "_m", "stats_data")

## Three groups only. COPD is dropped upstream in 03 (its 12 HLCA donors are all
## from Kaminski_2020, so a COPD estimate there is inseparable from that study);
## COPD is evaluated instead in panel E's independent dataset.
TRI      <- c("Healthy", "Fibrotic_ILD", "Other")
TRI_COL  <- c(Healthy = "#0072B2", Fibrotic_ILD = "#D55E00", Other = "#999999")
TRI_LAB  <- c(Healthy = "Healthy", Fibrotic_ILD = "Fibrotic/\nILD", Other = "Other\ndisease")

## Contrast vocabulary shared by the contrast panels so one colour always means
## one thing.
CTR_LAB <- c("Fibrotic_ILD - Healthy" = "Fibrotic/ILD",
             "Other - Healthy"        = "Other disease")
CTR_COL <- c("Fibrotic/ILD" = "#D55E00", "Other disease" = "#999999")

## Compact, unique study labels for the leave-one-study-out panel: 23 rows of
## "Banovich_Kropski_2020" do not fit a half-width panel. Drop the year token,
## then put it back only for the names that would otherwise collide
## (Schiller_2020 / Schiller_2021).
short_study <- function(x) {
    base <- vapply(strsplit(x, "_"), function(p)
        paste(p[!grepl("^(19|20)[0-9]{2}$", p)], collapse = "_"), character(1))
    base <- sub("_unpubl$", "", base)
    dup <- base %in% base[duplicated(base)]
    year <- sub(".*_((?:19|20)[0-9]{2}).*", "\\1", x)
    ifelse(dup & grepl("(19|20)[0-9]{2}", x), paste0(base, " ", year), base)
}

## ===========================================================================
## Panel A -- graded injury-programme engagement across three disease groups
## ===========================================================================
donor <- fread(file.path(SD, "donor_endpoint_table_3group.tsv"))
emm   <- fread(file.path(SD, "threegroup_emmeans.tsv"))
eff   <- fread(file.path(SD, "threegroup_effects.tsv"))

donor[, dx := factor(disease_group, levels = TRI)]
emm[,   dx := factor(disease_group, levels = TRI)]
setkey(eff, contrast)

n_lab <- donor[, .N, by = dx][order(dx)]
x_lab <- setNames(sprintf("%s\n(n = %d)", TRI_LAB[as.character(n_lab$dx)], n_lab$N),
                  as.character(n_lab$dx))

## Two significance brackets against Healthy, stacked so they do not overlap.
y_top <- max(donor$injury_program_score, na.rm = TRUE)
y_rng <- diff(range(donor$injury_program_score, na.rm = TRUE))
brk <- data.table(
    contrast = c("Fibrotic_ILD - Healthy", "Other - Healthy"),
    x = 1, xend = c(2, 3), y = y_top + c(0.08, 0.24) * y_rng)
brk <- merge(brk, eff[, .(contrast, estimate, p.value)], by = "contrast")
## fmt_p() is scalar-only, so map it rather than calling it on the column.
brk[, lab := sprintf("%+.2f SD   %s", estimate, vapply(p.value, fmt_p, character(1)))]

pA <- ggplot(donor, aes(dx, injury_program_score)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.3) +
    geom_jitter(aes(colour = dx), width = 0.13, height = 0, alpha = 0.5, size = 1.3) +
    geom_errorbar(data = emm, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  width = 0.16, linewidth = 0.6, colour = "black", inherit.aes = TRUE) +
    geom_point(data = emm, aes(y = emmean), size = 2.4, colour = "black") +
    geom_segment(data = brk, aes(x = x, xend = xend, y = y, yend = y),
                 linewidth = 0.3, inherit.aes = FALSE) +
    geom_text(data = brk, aes(x = (x + xend) / 2, y = y + 0.05 * y_rng, label = lab),
              size = 2.3, fontface = "bold", inherit.aes = FALSE) +
    scale_colour_manual(values = TRI_COL, guide = "none") +
    scale_x_discrete(labels = x_lab) +
    expand_limits(y = y_top + 0.34 * y_rng) +
    labs(x = NULL, y = "Pericyte injury-program score\n(SD units, sex + study-adjusted)") +
    theme_ms()

## ===========================================================================
## Panel B -- which continuous programs carry the score (both contrasts)
## ===========================================================================
nice_prog <- c(z_fibroblast_like = "Fibrillar\nfibroblast-like",
               z_activated_migratory = "Activated/\nmigratory",
               z_inflammatory = "Inflammatory",
               z_vascular_stabilizing = "Vascular-\nstabilizing")
comp <- fread(file.path(SD, "component_effects_3group.tsv"))
comp <- comp[response %in% names(nice_prog)]
## order by the Fibrotic/ILD effect: the panel's job is to show WHICH programme
## dominates, so the dominant one belongs at the top.
ord <- comp[contrast == "Fibrotic_ILD - Healthy"][order(estimate), response]
comp[, `:=`(prog = factor(nice_prog[response], levels = nice_prog[ord]),
            ctr  = factor(CTR_LAB[contrast], levels = CTR_LAB),
            sig  = ifelse(p.value < 0.05, "*", ""))]

pB <- ggplot(comp, aes(estimate, prog, colour = ctr)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(width = 0.55), height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.4, position = position_dodge(width = 0.55)) +
    geom_text(aes(x = ci_hi, label = sig), position = position_dodge(width = 0.55),
              hjust = -0.4, vjust = 0.75, size = 4, show.legend = FALSE) +
    scale_colour_manual(values = CTR_COL, name = NULL) +
    ## headroom on the right so the significance asterisks -- drawn just past
    ## ci_hi -- are not clipped by the panel edge, and so the inset legend clears
    ## the widest interval.
    expand_limits(x = max(comp$ci_hi) * 1.28) +
    labs(x = "Difference vs Healthy (SD units)", y = NULL) +
    theme_ms() +
    theme(panel.grid.major.y = element_blank(),
          legend.position = "inside", legend.position.inside = c(0.98, 0.04),
          legend.justification = c(1, 0), legend.background = element_blank(),
          legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"))

## ===========================================================================
## Panel C -- AGTR1 disease effects across stromal cell types
## ===========================================================================
## The x-axis is the OMNIBUS disease effect (partial eta-squared from the 2-df
## joint test of disease group), not a single contrast. That is the statistic the
## claim needs: a pericyte contrast estimate near zero with a wide CI is
## imprecision, whereas partial eta-squared is the share of donor-level AGTR1
## variance attributable to disease and is comparable across cell types with
## different n. Direction and CI for the Fibrotic/ILD contrast are in
## disease_association/_m/mean_expr/agtr1_celltype_disease_effects.tsv.
## ROW ORDER: lineage BLOCK first, then partial eta-squared descending within the
## block. The claim is a class-level comparison (fibroblast populations vs mural
## populations), so lineage has to be structural -- faceting makes the two blocks
## a property of the layout that survives any future reordering, whereas a single
## sorted column would leave the reader unable to tell whether "all fibroblasts on
## top" is a grouping or a coincidence of the sort. The strip labels also do the
## job of the colour legend, so the legend is dropped and the colours are left as
## redundant encoding.
rank_dt <- fread(file.path(MD, "agtr1_celltype_disease_ranking.tsv"))
LIN_COL <- c(Fibroblast = "#009E73", Mural = "#0072B2")
## Mesothelium became testable in the 2026-07-30 rebuild of 05 (the age fix admitted
## 15 fibrotic donors where there had been 1) and it is NEITHER fibroblast nor mural,
## so it cannot join either block without corrupting the arm it lands in. Drawn as a
## third one-row block it clips its own rotated strip label -- "Mesothelial" is wider
## than a single row is tall -- and it contributes nothing to the fibroblast-vs-mural
## contrast this panel exists to make, being the lowest eta^2 of anything tested
## (0.006, P = 0.96). So it is left out of the PANEL and kept in the outputs:
## agtr1_celltype_disease_ranking.tsv and supplementary table S13 both carry it, and
## the figure legend says so. Excluded here by lineage, not by name, so any future
## non-fibroblast/non-mural cell type is handled the same way rather than silently
## appearing in the wrong block.
n_drop <- rank_dt[!lineage %in% names(LIN_COL), .N]
if (n_drop) message(sprintf("panel C: %d non-fibroblast/mural cell type(s) held out: %s",
                            n_drop, paste(rank_dt[!lineage %in% names(LIN_COL), cell_type],
                                          collapse = ", ")))
rank_dt <- rank_dt[lineage %in% names(LIN_COL)]
## global ascending eta^2 -> within each facet the strongest cell type sits at the
## top, and the facets themselves are ordered Fibroblast above Mural.
rank_dt[, ct := factor(cell_type, levels = cell_type[order(partial_eta_sq)])]
rank_dt[, lin := factor(lineage, levels = names(LIN_COL))]

pC <- ggplot(rank_dt, aes(partial_eta_sq, ct, colour = lin)) +
    geom_segment(aes(x = 0, xend = partial_eta_sq, yend = ct), linewidth = 0.5) +
    geom_point(size = 2.6) +
    geom_text(aes(label = sprintf("P = %.2f", p_omnibus)), hjust = -0.35,
              size = 2.1, colour = "grey25", show.legend = FALSE) +
    facet_grid(lin ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_colour_manual(values = LIN_COL, guide = "none") +
    ## right-hand headroom for the direct P labels; the facet strips took width
    ## off the panel, so 0.28 clipped the longest label.
    scale_x_continuous(expand = expansion(mult = c(0, 0.40))) +
    labs(x = expression("Disease-attributable " * italic("AGTR1") * " variance (partial " * eta^2 * ")"),
         y = NULL) +
    theme_ms() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, face = "bold", size = 6.5),
          panel.spacing.y = unit(1.6, "mm"))

## ===========================================================================
## Panel D -- independent GSE136831 (COPD/IPF) evaluation of AGTR1
## ===========================================================================
## Estimates are already signed disease-minus-Control by 01.agtr1_copd_stats.R.
## Compartments that fail the >=5-donors-per-arm floor are NOT silently dropped:
## the pericyte row is the reason this dataset gives only qualified support, so
## it is drawn as an explicit "not estimable" marker rather than as a gap.
gse <- fread(file.path(GD, "agtr1_copd_all_contrasts.tsv"))[gene == "AGTR1"]
gse[, ctr := factor(contrast, levels = c("COPD", "IPF"))]
## Row order, top to bottom: the pre-specified fibroblast-lineage compartments
## first (that is the claim being evaluated), then the other powered
## compartments, and the NON-estimable pericyte row last -- the panel ends on the
## limitation rather than opening with it.
PRIMARY_CMP <- c("Fibroblast", "Myofibroblast")
cmp_rank <- gse[, .(is_est = any(estimable)), by = compartment]
cmp_rank[, grp := fifelse(!is_est, 1L, fifelse(compartment %in% PRIMARY_CMP, 3L, 2L))]
setorder(cmp_rank, grp, -compartment)   # ggplot draws level 1 at the bottom
gse[, cmp := factor(compartment, levels = cmp_rank$compartment)]

est <- gse[estimable == TRUE]
nes <- gse[estimable == FALSE]
GSE_COL <- c(COPD = "#E69F00", IPF = "#D55E00")
x_lo <- min(c(est$ci_lo, 0), na.rm = TRUE); x_hi <- max(c(est$ci_hi, 0), na.rm = TRUE)

pD <- ggplot(est, aes(estimate, cmp, colour = ctr)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   position = position_dodge(width = 0.6), height = 0.18, linewidth = 0.5) +
    geom_point(aes(shape = family), size = 2.3, position = position_dodge(width = 0.6)) +
    ## non-estimable rows: open crosses on the zero line, direct-labelled
    geom_point(data = nes, aes(x = 0, y = cmp), shape = 4, size = 1.8,
               colour = "grey55", inherit.aes = FALSE) +
    geom_text(data = unique(nes[, .(cmp, n_control_donors)]),
              aes(x = 0, y = cmp,
                  label = sprintf("  not estimable (%d Control donor%s)",
                                  n_control_donors, ifelse(n_control_donors == 1, "", "s"))),
              hjust = 0, vjust = 0.5, size = 2.0, colour = "grey40",
              inherit.aes = FALSE) +
    scale_colour_manual(values = GSE_COL, name = NULL) +
    ## breaks pinned: `family` is a character column, so the default legend order
    ## is alphabetical and puts "exploratory" ahead of the pre-specified family
    ## the panel is actually testing.
    scale_shape_manual(values = c(primary = 16, exploratory = 1), name = NULL,
                       breaks = c("primary", "exploratory"),
                       labels = c(primary = "pre-specified", exploratory = "exploratory")) +
    ## `est` carries no Pericyte row, and ggplot2 drops that unused level from the
    ## main layer and then appends it from the `nes` layer at the END of the
    ## scale -- which silently moves the not-estimable row to the top. Pin the
    ## row order to the intended one instead of relying on the factor levels.
    scale_y_discrete(limits = cmp_rank$compartment) +
    coord_cartesian(xlim = c(x_lo, x_hi * 1.05)) +
    ## Half-width panel since the 2026-07-29 four-panel rebuild: the axis title is
    ## shortened ("pseudobulk" is implied by the donor-level CI) and the legend
    ## moved from the right margin to the bottom, where it costs height (which
    ## this 8-row panel has) instead of width (which it has not). Kept on ONE
    ## line -- plotmath's atop() is the only way to wrap an expression and it
    ## opens a gap between the lines wide enough to read as a layout error.
    labs(x = expression("Disease minus Control " * italic("AGTR1") *
                        " (log1p CP10K, 95% CI)"),
         y = NULL) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme_ms() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6.5),
          legend.position = "bottom", legend.text = element_text(size = 6),
          legend.box = "horizontal", legend.box.spacing = unit(1, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(3, "mm"), legend.spacing.x = unit(1, "mm"))

## ===========================================================================
## SUPPLEMENT S16 -- study-level robustness of the panel-A effect (two panels)
## ===========================================================================
## Both panels defend the SAME estimate against two different objections, which
## is why they now travel together rather than one sitting in the main figure and
## one in the supplement:
##   S16A  "a single cohort drives the effect"      -> leave-one-study-out refits
##   S16B  "it is a between-study batch artifact"   -> within-study RE forest
## Neither substitutes for the other: LOSO shows no one study creates the effect
## but says nothing about whether it is reproduced INSIDE a study, and the forest
## shows the latter but only for the three studies that sampled both arms.
##
## SCALE: the two panels are NOT on a common scale. S16A refits the three-group
## model (programs z-standardised over Healthy+Fibrotic+Other); S16B is the
## two-group primary model (standardised over Healthy+Fibrotic only). Sitting in
## one figure they invite a numeric comparison that is not valid, so each axis
## title names its own standardisation set and no cross-panel annotation is drawn.

## ---- S16A: leave-one-study-out ---------------------------------------------
loso <- fread(file.path(SD, "leave_one_study_out_3group.tsv"))
loso <- loso[grepl("Fibrotic", contrast)]
loso[, lab := factor(short_study(dropped_dataset),
                     levels = short_study(dropped_dataset)[order(estimate)])]
full_est <- eff["Fibrotic_ILD - Healthy", estimate]
n_sig <- loso[p.value < 0.05, .N]

pLoso <- ggplot(loso, aes(estimate, lab)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = full_est, colour = "#D55E00", linewidth = 0.4, alpha = 0.7) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.4,
                   colour = "#0072B2") +
    geom_point(size = 1.3, colour = "#0072B2") +
    annotate("text", x = full_est, y = nrow(loso) + 0.9,
             label = sprintf("full data %+.2f SD", full_est),
             size = 2.2, colour = "#D55E00", hjust = 0.5, fontface = "bold") +
    ## sits in the margin BELOW the lowest refit, not on top of it
    annotate("text", x = max(loso$ci_hi), y = -0.05,
             label = sprintf("%d / %d refits P < 0.05", n_sig, nrow(loso)),
             size = 2.2, hjust = 1, fontface = "bold") +
    coord_cartesian(clip = "off", ylim = c(-0.4, nrow(loso) + 1.2)) +
    labs(x = "Fibrotic/ILD minus Healthy after dropping each study\n(SD of the Healthy + Fibrotic + Other donor set)",
         y = NULL) +
    theme_ms() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 5),
          plot.margin = margin(9, 4, 3, 3))

## ---- S16B: within-study random-effects forest ------------------------------
study <- fread(file.path(SD, "forest_per_study.tsv"))
pool  <- fread(file.path(SD, "forest_pooled_RE.tsv"))
study_lab <- sprintf("%s  (%d H / %d F)", study$dataset, study$nH, study$nF)
pool_lab  <- sprintf("RE pooled  (I2 = %.0f%%)", pool$I2)
fp <- rbind(
    data.table(label = study_lab, y = study$yi, lo = study$ci_lo, hi = study$ci_hi,
               w = study$weight_pct, kind = "study"),
    data.table(label = pool_lab, y = pool$estimate, lo = pool$ci_lo, hi = pool$ci_hi,
               w = max(study$weight_pct), kind = "pooled"))
fp[, label := factor(label, levels = rev(c(study_lab[order(study$yi)], pool_lab)))]

pFor <- ggplot(fp, aes(y = label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = kind), height = 0.2, linewidth = 0.5) +
    geom_point(aes(x = y, colour = kind, size = w, shape = kind)) +
    annotate("text", x = pool$estimate, y = 0.45,
             label = sprintf("%+.2f SD   %s", pool$estimate, fmt_p(pool$p_pooled)),
             size = 2.5, fontface = "bold", colour = "#D55E00") +
    scale_colour_manual(values = c(study = "#0072B2", pooled = "#D55E00"), guide = "none") +
    scale_shape_manual(values = c(study = 16, pooled = 18), guide = "none") +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    coord_cartesian(clip = "off", ylim = c(0.8, length(levels(fp$label)))) +
    labs(x = "Fibrotic/ILD minus Healthy within study\n(SD of the Healthy + Fibrotic donor set)",
         y = NULL) +
    theme_ms() +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size = 6.5))

## Stacked, not side by side: A has 23 rows and B has 4, so a single row would
## force one shared panel height and leave B's four rows floating in whitespace.
## The height weights keep the row pitch comparable between the two panels.
figS <- pLoso / pFor +
    plot_layout(heights = c(2.5, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10),
          ## one axis-title size for the whole figure. Both panels carry a long
          ## two-line title naming its standardisation set, which does not fit at
          ## theme_ms()'s 8 pt; setting it here rather than per panel keeps the
          ## two from drifting apart.
          axis.title = element_text(size = 7))
save_fig("figureS_disease_robustness", figS, 5.2, 6.0)

## ---- assemble --------------------------------------------------------------
## Reading order A,B / C,D: the graded effect and the programs that carry it on
## top, then the receptor-level result and its independent evaluation below.
## D gets the extra width -- 8 compartment rows, long "not estimable" in-panel
## text and a bottom legend, against C's 5 faceted rows. The widths are set on
## the inner row, since the outer `/` only controls the two row heights.
fig <- (pA | pB) / ((pC | pD) + plot_layout(widths = c(0.88, 1.12))) +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10),
          ## Axis-title size is set ONCE for all four panels rather than per
          ## panel. C and D need to be below theme_ms()'s 8 pt for their long
          ## titles to fit a half-width panel, and leaving A and B at 8 pt made
          ## the same element render at two visibly different sizes in one
          ## figure.
          axis.title = element_text(size = 7))
save_fig("figure_disease_main", fig, 7.5, 6.0)

cat("wrote", file.path(OUT, "figure_disease_main.{pdf,svg,png}"), "\n")
cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
