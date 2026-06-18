## Integrated pericyte-layer figure (Circ Research revision): ties the localization
## "where" (AGTR1 across the pericyte compartment) to the state-scoring "what" and
## the DPT "why" on ONE shared UMAP embedding (the embedding localization built and
## pericyte_states reused unchanged). No in-panel titles; captions carry meaning.
##
## Main  (figure_pericyte_layer): A subcluster UMAP, B AGTR1 expr UMAP, C dominant
##   state-program UMAP, D AGTR1 three-lens reversal (raw/detection/denoised; the
##   linchpin that AGTR1 is a compartment label, not a state marker), E DPT
##   pseudotime UMAP, F donor-level continuum trends.
## Supp  (figureS_pericyte_layer): per-program score UMAPs, ACTA2 expr + AGTR1
##   detection overlays (the contractile benchmark / dropout visual), and AGTR1 vs
##   ACTA2 donor-mean by subcluster.
##
## Reads: figures/_m/pericyte_umap_coords.tsv.gz (00.export_pericyte_umap.py),
##   pericyte_states/_m/pericytes_states_metadata.tsv.gz, continuum_metadata.tsv.gz,
##   stats_data/agtr1_lenses_by_program_emmeans.tsv, pseudotime_trend_correlations.tsv.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P <- function(...) file.path(ROOT, ...)
OUT <- P("figures", "mechanism"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- shared visual language (matches manuscript_mechanism_figure.R) -----
PROG_ORDER <- c("vascular_stabilizing", "fibroblast_like", "activated_migratory")
STATE_LEVELS <- c("vascular_stabilizing", "synthetic_contractile",
                  "activated_migratory", "inflammatory", "fibroblast_like")
STATE_LABS1 <- c(vascular_stabilizing = "Vascular-stabilizing",
                 synthetic_contractile = "Synthetic/contractile",
                 activated_migratory = "Activated/migratory",
                 inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like")
PROG_COL <- c(vascular_stabilizing = "#0072B2", fibroblast_like = "#D55E00",
              activated_migratory = "#CC79A7")
## subcluster hues grouped into program families (matches figureS_alluvial)
CLUST_COL <- c(P0 = "#08519C", P2 = "#6BAED6", P1 = "#D94801",
               P3 = "#FD8D3C", P5 = "#FDD0A2", P4 = "#CC79A7")
CLUST_ORDER <- c("P0", "P2", "P1", "P3", "P5", "P4")

theme_ms <- function(base = 8) {
    theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              strip.background = element_blank(),
              legend.position = "none")
}
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}

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

## D: AGTR1 three-lens reversal across programs (LINCHPIN). Centered within lens so
## the across-program PATTERN is comparable despite different native scales.
LENS_LABS <- c(AGTR1_expr = "AGTR1 (raw)", AGTR1_detect = "AGTR1 (detection)",
               AGTR1_scvi = "AGTR1 (denoised)")
LENS_COL  <- c("AGTR1 (raw)" = "#56B4E9", "AGTR1 (detection)" = "#999999",
               "AGTR1 (denoised)" = "#D55E00")
emm <- fread(P("pericyte_states", "_m", "stats_data", "agtr1_lenses_by_program_emmeans.tsv")) %>%
    filter(lens %in% names(LENS_LABS), state_program %in% PROG_ORDER) %>%
    mutate(lens = factor(LENS_LABS[lens], levels = LENS_LABS),
           program = factor(state_program, levels = PROG_ORDER)) %>%
    group_by(lens) %>% mutate(centered = emmean - mean(emmean)) %>% ungroup()
pD <- ggplot(emm, aes(program, centered, colour = lens, group = lens)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(aes(ymin = centered - SE, ymax = centered + SE), width = 0.12, linewidth = 0.4) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = LENS_COL, name = NULL) +
    scale_x_discrete(labels = STATE_LABS1[PROG_ORDER]) +
    labs(x = NULL, y = "Centered donor-aware emmean") +
    theme_ms() + theme(legend.position = c(0.5, 0.13), legend.background = element_blank(),
                       legend.text = element_text(size = 6), legend.key.size = unit(3, "mm"),
                       axis.text.x = element_text(angle = 20, hjust = 1))

## E: DPT pseudotime on the SAME embedding (WHY -- stabilizing<->injury continuum)
pE <- umap_cont(df, "dpt_pseudotime", "Pseudotime", option = "magma", dir = -1)

## F: donor-level continuum trends (programs & AGTR1 vs pseudotime)
nice_feat <- c(vascular_stabilizing = "Vascular-stabilizing",
               synthetic_contractile = "Synthetic/contractile",
               activated_migratory = "Activated/migratory",
               inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
               AGTR1 = "AGTR1")
trend <- fread(P("pericyte_states", "_m", "pseudotime_trend_correlations.tsv")) %>%
    filter(level == "donor") %>%
    mutate(feature = sub("_score$", "", feature),
           feature = recode(feature, AGTR1_expr = "AGTR1"),
           feature = recode(feature, !!!nice_feat),
           sig = ifelse(p_value < 0.05, "p < 0.05", "n.s.")) %>%
    arrange(spearman_rho)
trend$feature <- factor(trend$feature, levels = trend$feature)
pF <- ggplot(trend, aes(spearman_rho, feature, colour = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = spearman_rho, yend = feature), linewidth = 0.5) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("p < 0.05" = "#D55E00", "n.s." = "grey60"), name = NULL) +
    labs(x = "Spearman correlation (donor)", y = NULL) +
    theme_ms() + theme(legend.position = c(0.72, 0.18), legend.background = element_blank(),
                       legend.text = element_text(size = 6))

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
