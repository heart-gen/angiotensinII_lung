## Integrated pericyte-layer figure (Circ Research revision): ties the localization
## "where" (AGTR1 across the pericyte compartment) to the state-scoring "what" and
## the DPT "why" on ONE shared UMAP embedding (the embedding localization built and
## pericyte_states reused unchanged). No in-panel titles; captions carry meaning.
##
## Main  (figure_pericyte_layer): A subcluster UMAP, B AGTR1 expr UMAP, C dominant
##   state-program UMAP, D AGTR1 four-lens reversal across the six subclusters
##   (raw/detection/denoised/count-model arbiter; the linchpin that AGTR1 is a
##   compartment label, not a state marker), E DPT pseudotime UMAP, F donor-level
##   continuum trends.
## Supp  (figureS_pericyte_layer): per-program score UMAPs, ACTA2 expr + AGTR1
##   detection overlays (the contractile benchmark / dropout visual), and AGTR1 vs
##   ACTA2 donor-mean by subcluster.
##
## Reads: figures/_m/pericyte_umap_coords.tsv.gz (00.export_pericyte_umap.py),
##   pericyte_states/_m/pericytes_states_metadata.tsv.gz, continuum_metadata.tsv.gz,
##   pericyte_states/_m/stats_data/pseudotime_trend_correlations.tsv,
##   basement_membrane/_m/stats_data/agtr1_lens_by_cluster_emmeans.tsv
##     (19.agtr1_lens_by_cluster.R) and agtr1_count_by_cluster.tsv
##     (10.agtr1_count_models.R) -- panel D's four series, all at one unit.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

## ROOT/P/OUT, OKABE, DISEASE_*, theme_ms(), save_fig()
source("../_h/_fig_common.R")

## ---- visual language specific to this script ----------------------------
## PROG_ORDER = the three DOMINANT stable-cluster states (state_program). With the
## basement-membrane panel, clusters 1/3/5 relabel fibroblast_like -> basement_membrane
## (36% of pericytes), so BM replaces fibroblast_like here -- otherwise panels C/D/E
## (which filter state_program %in% PROG_ORDER) would silently drop those 4,200 cells.
## STATE_LEVELS/STATE_LABS1 keep fibroblast_like because panel F plots the continuous
## program SCORES (the fibrillar-ECM panel is a real, separate score); BM is added there too.
PROG_ORDER <- c("vascular_stabilizing", "basement_membrane", "activated_migratory")
STATE_LEVELS <- c("vascular_stabilizing", "synthetic_contractile",
                  "activated_migratory", "inflammatory", "fibroblast_like",
                  "basement_membrane")
STATE_LABS1 <- c(vascular_stabilizing = "Vascular-stabilizing",
                 synthetic_contractile = "Synthetic/contractile",
                 activated_migratory = "Activated/migratory",
                 inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
                 basement_membrane = "Basement-membrane")
PROG_COL <- c(vascular_stabilizing = "#0072B2", basement_membrane = "#D55E00",
              activated_migratory = "#CC79A7")
## subcluster hues grouped into program families (matches figureS_alluvial)
CLUST_COL <- c(P0 = "#08519C", P2 = "#6BAED6", P1 = "#D94801",
               P3 = "#FD8D3C", P5 = "#FDD0A2", P4 = "#CC79A7")
CLUST_ORDER <- c("P0", "P2", "P1", "P3", "P5", "P4")

## blank-axis UMAP base (the embedding axes carry no quantitative meaning)
umap_base <- function(df) ggplot(df, aes(UMAP1, UMAP2)) +
    theme_ms() + labs(x = "UMAP 1", y = "UMAP 2") +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
umap_cont <- function(df, fillvar, lab, option = "viridis", dir = 1) {
    umap_base(df) +
        geom_point(aes(colour = .data[[fillvar]]), size = 0.28, alpha = 0.7, stroke = 0) +
        scale_colour_viridis_c(option = option, direction = dir, name = lab) +
        theme(legend.position = "right", legend.key.height = unit(4, "mm"),
              legend.key.width = unit(2.5, "mm"), legend.title = element_text(size = 6.5),
              legend.text = element_text(size = 6))
}

## ---- load + join per-cell data on barcode ------------------------------
coords <- fread(P("figures", "_m", "pericyte_umap_coords.tsv.gz"))
meta <- fread(P("pericyte_states", "_m", "pericytes_states_metadata.tsv.gz"))
setnames(meta, 1, "barcode")
df <- merge(coords, meta, by = "barcode")
if (!"dpt_pseudotime" %in% names(df)) {
    cm <- fread(P("pericyte_states", "_m", "continuum_metadata.tsv.gz"))
    setnames(cm, 1, "barcode")
    df <- merge(df, cm[, .(barcode, dpt_pseudotime)], by = "barcode")
}
df <- df %>%
    mutate(cluster = factor(paste0("P", pericyte_state), levels = CLUST_ORDER),
           program = factor(state_program, levels = PROG_ORDER))

## ===== MAIN: figure_pericyte_layer ======================================
## A: subcluster UMAP (WHERE -- the population, six stable clusters)
pA <- umap_base(df) +
    geom_point(aes(colour = cluster), size = 0.28, alpha = 0.7, stroke = 0) +
    scale_colour_manual(values = CLUST_COL, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1)) +
    theme(legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6))

## B: AGTR1 expression on the SAME embedding (WHERE -- diffuse across compartment)
pB <- umap_cont(df, "AGTR1_expr", "AGTR1\n(log)", option = "viridis")

## C: dominant state-program (WHAT -- the three-program reframing)
pC <- umap_base(df) +
    geom_point(aes(colour = program), size = 0.28, alpha = 0.7, stroke = 0) +
    scale_colour_manual(values = PROG_COL, labels = STATE_LABS1[PROG_ORDER], name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1)) +
    theme(legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6))

## D: AGTR1 across the six stable subclusters under FOUR measurement lenses
## (LINCHPIN). Centered within lens so the across-cluster PATTERN is comparable
## despite different native scales.
##
## REBUILT 2026-09-10. Two changes, both about comparing like with like.
##
## (1) FOUR lenses, not three -- the count-model arbiter is now drawn, not merely
##     cited in the legend. It is the only readout that imputes nothing: AGTR1
##     integer counts as an NB-GLMM response with a library-size offset
##     (10.agtr1_count_models.R). The standing rule since 2026-09-02 is that the
##     count model, not the denoiser, arbitrates group contrasts, and a panel
##     whose whole claim is "the raw ordering is a dropout artifact" should show
##     the arbiter rather than assert it.
##
## (2) The x-axis is `pericyte_state` (P0-P5), NOT `state_program`. Programs are
##     assigned by a marker-panel argmax; the Leiden clusters come from 2,000 HVGs
##     that exclude AGTR1 (highly_variable = FALSE), so the grouping is provably
##     independent of the readout being compared across it. The count model's own
##     README names the by-cluster table as "the arbiter for any AGTR1-across-
##     clusters claim" for that reason. It is also where the evidence is: at
##     program level only 1 of 5 count-model specs separates basement-membrane
##     from vascular-stabilizing; at cluster level all twelve pseudobulk BM-vs-VS
##     contrasts across the three pseudobulk specs agree in sign, 10 of 12
##     significant. The cluster axis additionally ties D back to panel A.
##
## ALL FOUR SERIES ARE FIT AT THE SAME UNIT (214 donor x cluster pseudobulks,
## 95 donors, (1|study) + (1|donor_id) + depth covariate). That is the point of
## 19.agtr1_lens_by_cluster.R. The previous panel drew the denoised lens from a
## CELL-level lmer on 11,680 cells with only (1|donor) -- its SEs were 0.061-0.076
## on the log-rate scale against the count model's 0.131-0.204, a 2.0-2.7x gap
## that was entirely unit-of-analysis and would have read as the denoiser being
## the more precise measurement. Refit at the shared unit the same lens gives
## 0.110-0.221 and the two are within 5-43% of each other.
##
## TWO FACETS, because two of the four series are not on the log-rate scale and
## must never be read against the ones that are. Facet labels name the units.
LENS_LABS <- c(`raw AGTR1_expr` = "AGTR1 (raw)",
               `AGTR1_detect` = "AGTR1 (detection)",
               `denoised (retrained)` = "AGTR1 (denoised)",
               `AGTR1_count` = "AGTR1 (count model)")
LENS_COL  <- c("AGTR1 (raw)" = "#56B4E9", "AGTR1 (detection)" = "#999999",
               "AGTR1 (denoised)" = "#D55E00", "AGTR1 (count model)" = "#000000")
## Strip text names the UNITS so a reader cannot compare two series that are not
## comparable. Kept short: at this panel width (3 in) anything longer is clipped
## mid-word by the strip, which is worse than a terse label.
## "10k", not "10\u2074". The glyph itself is fine -- figureS_acta2_control renders
## it -- but that panel is 3.5 in wide and this one is 3.0, and at 3.0 the strip
## clipped mid-label to "log AGTR1 / 10... transcripts". Spelling the unit out is
## cheaper than widening the panel, which would squeeze E and F.
SCALE_OF  <- c("AGTR1 (raw)" = "log expr / detected fraction",
               "AGTR1 (detection)" = "log expr / detected fraction",
               "AGTR1 (denoised)" = "log AGTR1 per 10k transcripts",
               "AGTR1 (count model)" = "log AGTR1 per 10k transcripts")
SCALE_LEVELS <- c("log expr / detected fraction",
                  "log AGTR1 per 10k transcripts")

lens_dt <- norm_ci(fread(P("basement_membrane", "_m", "stats_data",
                           "agtr1_lens_by_cluster_emmeans.tsv")))[
    , .(cl = as.character(cl), lens, emmean, SE, n_units, underpowered)]

## The count arbiter: `spec = with_offset` is the concentration estimand (AGTR1
## per transcript sampled), the one whose units match the denoised lens. An
## unconverged fit is not a result and must not be drawn as one.
cnt_dt <- fread(P("basement_membrane", "_m", "stats_data",
                  "agtr1_count_by_cluster.tsv"))[
    level == "pseudobulk" & spec == "with_offset"]
if (nrow(cnt_dt) && !all(cnt_dt$converged)) {
    warning("figure_pericyte_layer panel D: dropping ", sum(!cnt_dt$converged),
            " unconverged count-model row(s)")
    cnt_dt <- cnt_dt[converged == TRUE]
}
## estimate is log(AGTR1 per transcript); + log(1e4) puts it on the per-10^4
## scale. Centering removes the constant anyway, but it keeps the facet honest.
cnt_dt <- cnt_dt[, .(cl = as.character(pericyte_state), lens = "AGTR1_count",
                     emmean = estimate + log(1e4), SE, n_units,
                     underpowered = as.logical(underpowered))]

emm <- rbind(lens_dt, cnt_dt, fill = TRUE) %>%
    mutate(readout = factor(LENS_LABS[lens], levels = LENS_LABS),
           cluster = factor(paste0("P", cl), levels = CLUST_ORDER),
           scale_grp = factor(SCALE_OF[as.character(readout)], levels = SCALE_LEVELS))
## Same failure mode the bm_relabel guard was written for: a lens or a cluster
## that vanished during the join must error, not be plotted around silently.
require_programs(emm$readout, unname(LENS_LABS), "figure_pericyte_layer panel D")
require_programs(emm$cluster, CLUST_ORDER, "figure_pericyte_layer panel D (clusters)")
emm <- emm %>%
    group_by(readout) %>% mutate(centered = emmean - mean(emmean)) %>% ungroup()

## P4 (13 donors, 134 cells) and P5 (4 donors, 44 cells) carry the count model's
## `underpowered` flag. Marked from the DATA rather than hardcoded, so a change
## upstream cannot leave the mark on the wrong cluster. A plain asterisk, not a
## dagger: the dagger widened the tick label past the axis clip, so P5/P4 rendered
## as "P5..."/"P4...".
low <- emm %>% filter(underpowered) %>% pull(cluster) %>% unique() %>% as.character()
clust_lab <- setNames(ifelse(CLUST_ORDER %in% low, paste0(CLUST_ORDER, "*"),
                            CLUST_ORDER), CLUST_ORDER)

pD <- ggplot(emm, aes(cluster, centered, colour = readout, group = readout)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(aes(ymin = centered - SE, ymax = centered + SE), width = 0.12,
                  linewidth = 0.4) +
    geom_point(size = 1.9) +
    facet_wrap(~ scale_grp, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = LENS_COL, name = NULL) +
    ## Four readouts on one row overflow the panel width and clip the leading key
    ## glyph; two rows fit.
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_x_discrete(labels = clust_lab) +
    labs(x = NULL, y = "Centered pseudobulk marginal mean") +
    theme_ms() + theme(legend.position = "bottom", legend.background = element_blank(),
                       legend.text = element_text(size = 6),
                       legend.key.size = unit(3, "mm"), legend.margin = margin(t = -4),
                       strip.text = element_text(size = 6))

## E: DPT pseudotime on the SAME embedding (WHY -- stabilizing<->basement-membrane axis)
pE <- umap_cont(df, "dpt_pseudotime", "Pseudotime", option = "magma", dir = -1)

## F: donor-level continuum trends (programs & AGTR1 vs pseudotime)
## BOTH AGTR1 lenses are plotted. Raw AGTR1 and the program scores share a
## sequencing-depth gradient, which is what 03.agtr1_lenses.R showed manufactures
## AGTR1's apparent program bias; showing the raw trend alone here would invite
## exactly the reading panel D exists to refute. The two lens rows carry the panel
## D lens colours as a ring so the reader connects the panels.
LENS_FEAT <- c(AGTR1_expr = "AGTR1 (raw)", AGTR1_scvi = "AGTR1 (denoised)")
SIG_COL <- c("p < 0.05" = "#D55E00", "n.s." = "grey60")
nice_feat <- c(vascular_stabilizing = "Vascular-stabilizing",
               synthetic_contractile = "Synthetic/contractile",
               activated_migratory = "Activated/migratory",
               inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
               basement_membrane = "Basement-membrane", LENS_FEAT)
trend <- fread(P("pericyte_states", "_m", "pseudotime_trend_correlations.tsv")) %>%
    filter(level == "donor") %>%
    mutate(feature = sub("_score$", "", feature),
           feature = recode(feature, !!!nice_feat),
           sig = ifelse(p_value < 0.05, "p < 0.05", "n.s.")) %>%
    arrange(spearman_rho)
## Same failure mode as panel D: a missing lens would silently drop a row rather
## than error, and the panel's whole point is the raw-vs-denoised comparison.
require_programs(trend$feature, unname(LENS_FEAT),
                 "figure_pericyte_layer panel F (AGTR1 lenses)")
trend$feature <- factor(trend$feature, levels = trend$feature)
## Colour encodes significance only, exactly as on every other row -- the two lens
## rows are identified by their axis labels, which match panel D's legend text
## verbatim. An earlier version outlined them in the panel D lens colours, but the
## denoised ring re-used the significance orange, so a non-significant denoised
## point acquired an orange outline and read as significant.
pF <- ggplot(trend, aes(spearman_rho, feature, colour = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = spearman_rho, yend = feature), linewidth = 0.5) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = SIG_COL, name = NULL) +
    labs(x = "Spearman correlation (donor)", y = NULL) +
    theme_ms() +
    ## Legend at the bottom: with eight rows there is no interior gap that stays
    ## clear as the estimates move between runs.
    theme(legend.position = "bottom", legend.background = element_blank(),
          legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"),
          legend.margin = margin(t = -4))

main <- (pA | pB | pC) / (pD | pE | pF) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figure_pericyte_layer", main, 9.0, 6.0)

## ===== SUPPLEMENT: figureS_pericyte_layer ===============================
## sA: per-program score UMAPs (small multiples, shared embedding)
score_cols <- paste0(STATE_LEVELS, "_score")
sl <- df %>%
    select(UMAP1, UMAP2, all_of(score_cols)) %>%
    pivot_longer(all_of(score_cols), names_to = "program", values_to = "score") %>%
    mutate(program = factor(STATE_LABS1[sub("_score$", "", program)],
                            levels = STATE_LABS1[STATE_LEVELS]))
sA <- ggplot(sl, aes(UMAP1, UMAP2, colour = score)) +
    geom_point(size = 0.18, alpha = 0.7, stroke = 0) +
    facet_wrap(~ program, nrow = 1) +
    scale_colour_viridis_c(option = "viridis", name = "Score") +
    labs(x = "UMAP 1", y = "UMAP 2") + theme_ms() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), legend.position = "right",
          legend.key.height = unit(4, "mm"), legend.key.width = unit(2.5, "mm"),
          legend.title = element_text(size = 6.5), legend.text = element_text(size = 6),
          strip.text = element_text(size = 6.5))

## sB: ACTA2 expression (contractile benchmark); sC: AGTR1 detection (dropout visual)
sB <- umap_cont(df, "ACTA2_expr", "ACTA2\n(log)", option = "viridis")
sC <- umap_cont(df %>% mutate(AGTR1_detect = as.numeric(AGTR1_detect)),
                "AGTR1_detect", "AGTR1\ndetected", option = "cividis")

## sD: AGTR1 vs ACTA2 donor-mean by subcluster (donor-aware distribution)
ds <- df %>%
    group_by(donor_id, cluster) %>%
    summarise(AGTR1 = mean(AGTR1_expr, na.rm = TRUE),
              ACTA2 = mean(ACTA2_expr, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>%
    pivot_longer(c(AGTR1, ACTA2), names_to = "gene", values_to = "expr")
sD <- ggplot(ds, aes(cluster, expr, fill = gene)) +
    geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85, linewidth = 0.3,
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = c(AGTR1 = "#56B4E9", ACTA2 = "#009E73"), name = NULL) +
    labs(x = NULL, y = "Donor-mean expression (log)") + theme_ms() +
    theme(legend.position = c(0.85, 0.85), legend.background = element_blank(),
          legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"))

supp <- sA / (sB | sC | sD) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_pericyte_layer", supp, 9.0, 5.6)

cat("pericyte-layer figures written to", OUT, "\n")
sessioninfo::session_info()
