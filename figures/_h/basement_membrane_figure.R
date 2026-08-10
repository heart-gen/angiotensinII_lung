## Basement-membrane program and AGT ligand-axis figures (Circ Research revision).
##
## Main: which basement-membrane components pericytes selectively deposit, how the
## BM axis separates from fibrillar ECM, and the state-model consequence.
## Supplements: the local RAS landscape, and the COPD contrast with its power bound.
## No in-panel titles; interpretation belongs in the caption.
##
## 2026-08-10 revision (collaborator request). The main figure now separates
## basement-membrane matrix from fibril-forming interstitial collagen explicitly:
## the dot plot is blocked into BM / fibrillar I-III / fibrillar V-XI / ambient
## tracer, the endpoint panel contrasts BM against the fibril-forming collagens
## rather than the mixed fibrillar panel, and two new panels carry the question
## the split raises -- is the residual pericyte fibrillar signal transcribed
## (collagen I stoichiometry, panel D) and is it above soup (ambient controls,
## panel E). Grew from 4 panels to 7; the role palette is new and shared by
## panels C-E.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

## ROOT/P/OUT, OKABE, theme_ms(), save_fig(), tag()
source("../_h/_fig_common.R")

## This script keeps ggplot's default legend placement -- its heatmap/dotplot
## panels need their colour and size guides.
theme_ms <- function(base = 8) .theme_ms(base = base, legend = NULL)

BM_M <- P("basement_membrane", "_m")
AGT_M <- P("agt_axis", "_m")
rd <- function(f) if (file.exists(f)) fread(f) else NULL

## ---------------------------------------------------------------- vocabulary ----
## Gene display blocks come from basement_membrane/_m/bm_panel_genes.tsv (written
## by bm_panels.py), so the figure never re-declares which gene is which -- the
## 2026-08-10 revision exists precisely because COL4A1 had been living in a
## fibrillar panel, and a second hard-coded copy here would invite the same bug.
panel_tbl <- rd(file.path(BM_M, "bm_panel_genes.tsv"))
if (is.null(panel_tbl))
    stop("bm_panel_genes.tsv missing; run basement_membrane step_0 first")
BLOCK_LEVELS <- c("basement_membrane", "fibrillar_core", "fibrillar_minor",
                  "ambient_tracer")
blk <- unique(panel_tbl[, .(gene, block, block_index)])[order(block_index)]
## `fibrillar_other` -- FN1/POSTN/LUM/DCN/FBN1/BGN, the non-collagen members of
## the frozen fibrillar_ecm panel -- is deliberately NOT displayed. Those genes
## still drive the frozen score, but the collaborator's question is specifically
## about collagen, and showing them would put a fifth unlabelled block on the
## dot plot without adding evidence about fibril formation.
blk <- blk[block %in% BLOCK_LEVELS]
GENE_ORDER <- blk$gene
BLOCK_OF <- setNames(blk$block, blk$gene)

BLOCK_LABS <- c(basement_membrane = "Basement membrane",
                fibrillar_core    = "Fibrillar I/III",
                fibrillar_minor   = "Fibrillar V/XI",
                ambient_tracer    = "Ambient tracer")
BLOCK_COL <- c(basement_membrane = "#0072B2", fibrillar_core = "#D55E00",
               fibrillar_minor = "#E69F00", ambient_tracer = "#999999")
blk_factor <- function(x) factor(x, levels = BLOCK_LEVELS, labels = BLOCK_LABS)

## Cell-type roles. Mirrors the definition in basement_membrane/_h/08.fibrillar_ambient.R;
## panels that read that script's outputs take `role` from the data, and only the
## panels fed by other scripts derive it here.
AMBIENT_REFERENCE <- c("EC general capillary", "EC aerocyte capillary",
                       "EC venous pulmonary", "EC venous systemic",
                       "Alveolar macrophages", "Interstitial macrophages")
role_of <- function(g) fifelse(
    g == "Pericytes", "Pericytes",
    fifelse(grepl("fibroblast|Myofibro", g, ignore.case = TRUE), "Fibroblasts",
    fifelse(g %in% AMBIENT_REFERENCE, "Ambient reference", "Other")))
ROLE_LEVELS <- c("Pericytes", "Fibroblasts", "Ambient reference", "Other")
ROLE_COL <- c(Pericytes = "#009E73", Fibroblasts = "#CC79A7",
              `Ambient reference` = "#999999", Other = "grey80")
role_factor <- function(x) factor(x, levels = ROLE_LEVELS)

## Cell types shown in the dot plot. The full 22-group profile drives every
## statistic; the panel shows the interpretable subset -- pericytes, the mural
## neighbour, every fibroblast class, and the ambient-reference lineages -- so the
## block inversion is legible instead of being buried in 22 columns of immune cells.
FOCUS_GROUPS <- c("Pericytes", "Vascular smooth muscle",
                  "Adventitial fibroblasts", "Alveolar fibroblasts",
                  "Peribronchial fibroblasts", "Myofibroblasts",
                  "Subpleural fibroblasts",
                  "EC general capillary", "EC aerocyte capillary",
                  "Alveolar macrophages", "AT2_AGTR2det", "AT1")

## ============ A: matrix-block x cell-type dot plot (the inversion) ==========
prof <- rd(file.path(BM_M, "stats_data", "bm_celltype_profile.tsv"))
pbk <- rd(file.path(BM_M, "bm_pseudobulk_celltype.tsv.gz"))
pA <- NULL
if (!is.null(prof)) {
    d <- prof[gene %in% GENE_ORDER & value_type == "expr"]
    ## z within gene across ALL profiled cell types, then subset for display, so
    ## the colour scale reports each gene's position in the whole lung rather
    ## than within the curated subset.
    d[, z := (emmean - mean(emmean)) / (sd(emmean) + 1e-9), by = gene]
    if (!is.null(pbk)) {
        dcols <- intersect(paste0(GENE_ORDER, "__detect"), names(pbk))
        det <- pbk[n_cells >= 5, lapply(.SD, mean, na.rm = TRUE),
                   by = ccc_group, .SDcols = dcols]
        det <- melt(det, id.vars = "ccc_group", variable.name = "gene",
                    value.name = "detect")
        det[, gene := sub("__detect$", "", gene)]
        d <- merge(d, det, by = c("ccc_group", "gene"), all.x = TRUE)
    } else {
        d[, detect := 0.5]
    }
    d <- d[ccc_group %in% FOCUS_GROUPS]
    d[, gene := factor(gene, levels = GENE_ORDER)]
    d[, block := blk_factor(BLOCK_OF[as.character(gene)])]
    ## Fixed semantic row order, shared with panels C-E, so the reader can carry
    ## one mental ordering across the figure.
    d[, ccc_group := factor(ccc_group, levels = rev(FOCUS_GROUPS))]
    pA <- ggplot(d, aes(gene, ccc_group)) +
        geom_point(aes(size = detect, fill = z), shape = 21, stroke = 0.15,
                   colour = "grey30") +
        facet_grid(. ~ block, scales = "free_x", space = "free_x") +
        scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                             midpoint = 0, name = "Expression\n(z within gene)") +
        scale_size_continuous(range = c(0.3, 3.4), name = "Detected\nfraction",
                              limits = c(0, 1)) +
        labs(x = NULL, y = NULL) +
        theme_ms(7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right", legend.key.size = unit(0.3, "cm"),
              panel.spacing.x = unit(0.18, "cm"),
              strip.text = element_text(face = "bold", size = 6.5))
}

## ============ B: selectivity, pericyte vs next-highest cell type ============
tau <- rd(file.path(BM_M, "stats_data", "bm_tau_specificity.tsv"))
pB <- NULL
if (!is.null(tau)) {
    ## Tracers are excluded here: they are off-lineage by construction, so their
    ## ratios are uninformative and would compress the axis. They carry the
    ## ambient argument in panel E instead.
    d <- tau[value_type == "expr" & gene %in% GENE_ORDER &
             BLOCK_OF[gene] != "ambient_tracer"]
    setorder(d, log2_pericyte_over_next)
    d[, gene := factor(gene, levels = d$gene)]
    d[, block := blk_factor(BLOCK_OF[as.character(gene)])]
    pB <- ggplot(d, aes(log2_pericyte_over_next, gene, fill = block)) +
        geom_col(width = 0.7, colour = "grey30", linewidth = 0.15) +
        geom_vline(xintercept = 0, linewidth = 0.3) +
        geom_text(aes(label = sprintf("#%d", pericyte_rank),
                      hjust = ifelse(log2_pericyte_over_next > 0, -0.25, 1.25)),
                  size = 2.1) +
        ## drop = TRUE: the tracer genes are excluded from this panel, so their
        ## key would advertise evidence the panel does not contain.
        scale_fill_manual(values = setNames(BLOCK_COL[BLOCK_LEVELS],
                                            BLOCK_LABS[BLOCK_LEVELS]),
                          name = NULL, drop = TRUE) +
        guides(fill = guide_legend(title = NULL)) +
        scale_x_continuous(expand = expansion(mult = 0.18)) +
        labs(x = expression(log[2]~"(pericyte / next-highest cell type)"), y = NULL) +
        theme_ms(7) +
        theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))
}

## ====== C: BM minus fibril-forming collagen, by cell type (the endpoint) =====
## The frozen `fibrillar_ecm` endpoint mixes collagens with proteoglycans and so
## answers "is this cell fibroblast-like"; the block endpoint asks the narrower
## question the revision is about -- does this cell build interstitial fibrils.
blkend <- rd(file.path(BM_M, "stats_data", "bm_block_endpoint_emmeans.tsv"))
pe <- rd(file.path(BM_M, "stats_data", "bm_primary_endpoint_emmeans.tsv"))
pC <- NULL
d <- NULL
if (!is.null(blkend) && "endpoint" %in% names(blkend) &&
    any(blkend$endpoint == "bm_minus_fibrillar_core")) {
    d <- as.data.table(blkend)[endpoint == "bm_minus_fibrillar_core"]
    xlab_C <- "BM - fibrillar I/III score (panel z difference)"
} else if (!is.null(pe)) {
    ## Falls back to the frozen endpoint if step_2 has not been re-run yet.
    d <- as.data.table(pe)
    xlab_C <- "BM - fibrillar score (panel z difference)"
}
if (!is.null(d)) {
    setorder(d, emmean)
    d[, ccc_group := factor(ccc_group, levels = d$ccc_group)]
    d[, role := role_factor(role_of(as.character(ccc_group)))]
    pC <- ggplot(d, aes(emmean, ccc_group, colour = role)) +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
        geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0,
                      orientation = "y", linewidth = 0.4) +
        geom_point(size = 1.5) +
        scale_colour_manual(values = ROLE_COL, name = NULL, drop = FALSE) +
        labs(x = xlab_C, y = NULL) +
        theme_ms(7) +
        theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))
}

## ======= D: collagen I heterotrimer stoichiometry (COL1A1 vs COL1A2) ========
## Collagen I is alpha1(I)2 alpha2(I)1. Ambient transfer preserves the source
## cell's ratio, so a pericyte-specific departure from the fibroblast ratio is
## evidence that the signal is transcribed -- and transcribed incompletely.
st_emm <- rd(file.path(BM_M, "stats_data", "fibrillar_stoichiometry_emmeans.tsv"))
st_don <- rd(file.path(BM_M, "stats_data", "fibrillar_stoichiometry_donor.tsv"))
pD <- NULL
if (!is.null(st_emm) && nrow(st_emm)) {
    e <- as.data.table(st_emm)
    e[, role := role_factor(if ("role" %in% names(e)) role else
                            role_of(as.character(ccc_group)))]
    setorder(e, emmean)
    lev <- as.character(e$ccc_group)
    e[, ccc_group := factor(ccc_group, levels = lev)]
    pD <- ggplot(e, aes(emmean, ccc_group, colour = role))
    if (!is.null(st_don) && nrow(st_don)) {
        dd <- as.data.table(st_don)[ccc_group %in% lev]
        dd[, ccc_group := factor(ccc_group, levels = lev)]
        dd[, role := role_factor(role)]
        pD <- pD + geom_point(data = dd, aes(stoich, ccc_group), size = 0.35,
                              alpha = 0.28, stroke = 0,
                              position = position_jitter(height = 0.18, width = 0))
    }
    pD <- pD +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey40") +
        geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0,
                      orientation = "y", linewidth = 0.45) +
        geom_point(size = 1.7) +
        annotate("text", x = 0, y = length(lev) + 0.75, label = "balanced",
                 size = 1.9, colour = "grey40", vjust = 0) +
        scale_colour_manual(values = ROLE_COL, name = NULL, drop = FALSE) +
        scale_y_discrete(expand = expansion(add = c(0.6, 1.4))) +
        labs(x = "COL1A1 - COL1A2 within unit (log1p CP10K)", y = NULL) +
        theme_ms(7) +
        theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))
}

## ================== E: ambient controls, two independent tests ==============
## Left facet: how far each gene sits above the co-dissociated non-collagen
## lineages (the soup floor). Right facet: whether the pericyte value tracks that
## donor's measured soup burden, which is what ambient contamination requires.
## The tracer genes calibrate both: they are pure soup, so they mark where a
## contaminated gene should land.
flr <- rd(file.path(BM_M, "stats_data", "fibrillar_ambient_floor.tsv"))
areg <- rd(file.path(BM_M, "stats_data", "fibrillar_ambient_regression.tsv"))
pE <- NULL
if (!is.null(flr) && nrow(flr)) {
    f <- as.data.table(flr)[value_type == "expr"]
    f <- if ("is_pericyte_contrast" %in% names(f)) f[is_pericyte_contrast == TRUE] else
        f[grepl("Pericytes", contrast)]
    f <- f[block %in% c("fibrillar_core", "fibrillar_minor", "ambient_tracer")]
    f <- f[, .(gene, block, estimate, lo = estimate - 1.96 * SE,
               hi = estimate + 1.96 * SE,
               facet = "Above ambient floor\n(pericyte - EC/macrophage)")]
    parts <- list(f)
    if (!is.null(areg) && nrow(areg)) {
        a <- as.data.table(areg)[term == "soup_burden"]
        if (nrow(a))
            parts <- c(parts, list(a[, .(gene, block, estimate,
                                         lo = estimate - 1.96 * SE,
                                         hi = estimate + 1.96 * SE,
                                         facet = "Tracks donor soup burden\n(slope)")]))
    }
    g <- rbindlist(parts, fill = TRUE)
    ## One gene order across both facets, set by the floor gap: the reader should
    ## be able to run a finger horizontally from "high above floor" to "does not
    ## track soup" without re-reading the axis.
    ord <- f[order(estimate), gene]
    g[, gene := factor(gene, levels = ord)]
    g[, block := blk_factor(block)]
    g <- g[!is.na(gene)]
    pE <- ggplot(g, aes(estimate, gene, colour = block)) +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
        geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y",
                      linewidth = 0.4) +
        geom_point(size = 1.5) +
        facet_wrap(~ facet, scales = "free_x") +
        ## drop = TRUE here, unlike panels A/B: no basement-membrane gene appears
        ## in this panel, and keeping the level would render a BM key for a panel
        ## that carries no BM evidence.
        scale_colour_manual(values = setNames(BLOCK_COL[BLOCK_LEVELS],
                                              BLOCK_LABS[BLOCK_LEVELS]),
                            name = NULL, drop = TRUE) +
        guides(colour = guide_legend(title = NULL)) +
        labs(x = "Estimate (log1p CP10K)", y = NULL) +
        theme_ms(7) +
        theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
              strip.text = element_text(size = 6.2))
}

## ================ F: the two matrix programs along the continuum ============
crho <- rd(file.path(BM_M, "stats_data", "bm_continuum_donor_rho.tsv"))
csum <- rd(file.path(BM_M, "stats_data", "bm_continuum_summary.tsv"))
pF <- NULL
if (!is.null(crho) && nrow(crho)) {
    METRIC_LABS <- c(rho_bm = "Basement\nmembrane", rho_core = "Fibrillar\nI/III",
                     rho_minor = "Fibrillar\nV/XI",
                     rho_switch_core = "BM - fibrillar\n(switch)",
                     rho_tracer = "Ambient\ntracer")
    mm <- intersect(names(METRIC_LABS), names(crho))
    d <- melt(as.data.table(crho), id.vars = "donor_id", measure.vars = mm,
              variable.name = "metric", value.name = "rho")[is.finite(rho)]
    d[, metric := factor(metric, levels = mm, labels = METRIC_LABS[mm])]
    med <- d[, .(rho = median(rho)), by = metric]
    pF <- ggplot(d, aes(metric, rho)) +
        geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
        geom_violin(fill = "grey92", colour = NA, width = 0.85) +
        geom_jitter(width = 0.13, height = 0, size = 0.3, alpha = 0.35, stroke = 0) +
        geom_point(data = med, size = 1.9, colour = "#D55E00") +
        labs(x = NULL, y = expression("Per-donor Spearman"~rho~"vs pseudotime")) +
        theme_ms(7)
    ## BH stars from the donor-level Wilcoxon, placed above each violin.
    if (!is.null(csum) && "p_BH" %in% names(csum)) {
        s <- as.data.table(csum)[metric %in% mm]
        s[, lab := fifelse(p_BH < 0.001, "***",
                   fifelse(p_BH < 0.01, "**",
                   fifelse(p_BH < 0.05, "*", "n.s.")))]
        s[, metric := factor(metric, levels = mm, labels = METRIC_LABS[mm])]
        s[, y := max(d$rho, na.rm = TRUE) * 1.06]
        pF <- pF + geom_text(data = s, aes(metric, y, label = lab), size = 2.2,
                             inherit.aes = FALSE)
    }
}

## ============= G: the state-model consequence (the gate result) =============
gate <- rd(file.path(BM_M, "state_gate_relenrich.tsv"))
pG <- NULL
if (!is.null(gate)) {
    d <- melt(gate, id.vars = "pericyte_state", variable.name = "program",
              value.name = "relenrich")
    d[, program := sub("_relenrich$", "", program)]
    d[, program := factor(program, levels = c("vascular_stabilizing", "inflammatory",
                                              "synthetic_contractile",
                                              "activated_migratory",
                                              "fibroblast_like", "basement_membrane"))]
    d[, pericyte_state := factor(pericyte_state)]
    ## Mark the winning program per cluster under the 6-panel model.
    win <- d[, .SD[which.max(relenrich)], by = pericyte_state]
    pG <- ggplot(d, aes(program, pericyte_state, fill = relenrich)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_point(data = win, shape = 8, size = 1.4, colour = "black") +
        scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                             midpoint = 0, name = "Relative\nenrichment") +
        scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
        labs(x = NULL, y = "Pericyte cluster") +
        theme_ms(7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right", legend.key.size = unit(0.3, "cm"),
              panel.grid = element_blank())
}

## ------------------------------------------------------------- assemble ----
## Narrative order: what the compartment makes (A), how selectively (B), the
## endpoint that survives capture confounding (C), the stoichiometric evidence
## that the fibrillar residue is transcribed (D), the ambient controls behind
## that claim (E), the behaviour along the continuum (F), and the consequence for
## the state model (G). A spans the width because the block inversion is the
## figure's single most important read.
have <- !vapply(list(pA, pB, pC, pD, pE, pF, pG), is.null, logical(1))
if (any(have)) {
    labs_all <- LETTERS[1:7]
    plots <- Filter(Negate(is.null), list(pA, pB, pC, pD, pE, pF, pG))
    plots <- Map(tag, plots, labs_all[have])
    if (all(have)) {
        design <- "AA
                   BC
                   DE
                   FG"
        fig <- wrap_plots(plots, design = design, heights = c(1.15, 1.25, 1, 0.95))
        save_fig("figure_basement_membrane", fig, 11, 13.5)
    } else {
        ## Partial run (an upstream step has not landed yet): keep it renderable
        ## rather than failing, but the layout above is the publication one.
        warning("only ", sum(have), "/7 panels available; emitting a fallback grid",
                call. = FALSE)
        fig <- wrap_plots(plots, ncol = 2)
        save_fig("figure_basement_membrane", fig, 11, 3.4 * ceiling(sum(have) / 2))
    }
    message("wrote figure_basement_membrane (", sum(have), " panels)")
}

## ================== Supplement: local RAS landscape =======================
ras <- rd(file.path(AGT_M, "stats_data", "ras_celltype_profile.tsv"))
comp <- rd(file.path(AGT_M, "stats_data", "ras_circuit_completeness.tsv"))
if (!is.null(ras)) {
    RAS_ORDER <- c("AGT", "REN", "ACE", "ACE2", "CMA1", "CTSG", "CTSD", "ENPEP",
                   "MME", "AGTR1", "AGTR2", "LRP2", "MAS1", "TGFB1", "TGFB2")
    d <- as.data.table(ras)[gene %in% RAS_ORDER]
    d[, gene := factor(gene, levels = rev(RAS_ORDER))]
    ord <- d[gene == "AGT"][order(-detect), ccc_group]
    d[, ccc_group := factor(ccc_group, levels = ord)]
    s1 <- ggplot(d, aes(ccc_group, gene, fill = detect)) +
        geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_viridis_c(option = "magma", direction = -1,
                             name = "Detected\nfraction", trans = "sqrt") +
        labs(x = NULL, y = NULL) +
        theme_ms(7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right", legend.key.size = unit(0.3, "cm"),
              panel.grid = element_blank())
    s2 <- NULL
    if (!is.null(comp) && all(c("has_substrate", "has_protease", "has_receptor")
                              %in% names(comp))) {
        ## A count-of-steps bar is uninformative here (almost every cell type
        ## carries exactly one step). Showing the three circuit REQUIREMENTS
        ## side by side makes the actual point legible: no row is complete.
        cc <- as.data.table(comp)
        req <- melt(cc[, .(ccc_group, Substrate = has_substrate,
                           Protease = has_protease, `AT1R` = has_receptor)],
                    id.vars = "ccc_group", variable.name = "requirement",
                    value.name = "present")
        ord <- req[, .(n = sum(present)), by = ccc_group][order(-n), ccc_group]
        req[, ccc_group := factor(ccc_group, levels = rev(ord))]
        s2 <- ggplot(req, aes(requirement, ccc_group, fill = present)) +
            geom_tile(colour = "white", linewidth = 0.5) +
            scale_fill_manual(values = c(`TRUE` = OKABE[4], `FALSE` = "grey88"),
                              name = "Present", labels = c("no", "yes")) +
            labs(x = NULL, y = NULL) +
            theme_ms(7) +
            theme(legend.position = "right", legend.key.size = unit(0.3, "cm"),
                  panel.grid = element_blank())
    }
    sup <- if (is.null(s2)) tag(s1, "A") else
        wrap_plots(list(tag(s1, "A"), tag(s2, "B")), nrow = 1, widths = c(1.4, 1))
    save_fig("figureS_ras_landscape", sup, 11, 4.6)
    message("wrote figureS_ras_landscape")
}

## ================== Supplement: COPD contrast + power bound ================
prim <- rd(file.path(BM_M, "stats_data", "bm_copd_primary.tsv"))
pw <- rd(file.path(BM_M, "stats_data", "bm_pericyte_power.tsv"))
if (!is.null(prim) && nrow(prim)) {
    d <- as.data.table(prim)[grepl("COPD", contrast)]
    if (nrow(d)) {
        d[, lab := paste0(gene, " - ", compartment)]
        setorder(d, estimate)
        d[, lab := factor(lab, levels = d$lab)]
        d[, sig := !is.na(p_BH) & p_BH < 0.05]
        q1 <- ggplot(d, aes(estimate, lab, colour = sig)) +
            geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3,
                       colour = "grey50") +
            geom_errorbar(aes(xmin = estimate - 1.96 * SE,
                               xmax = estimate + 1.96 * SE), width = 0, orientation = "y",
                           linewidth = 0.4) +
            geom_point(size = 1.6) +
            scale_colour_manual(values = c(`TRUE` = OKABE[4], `FALSE` = "grey55"),
                                name = "BH < 0.05") +
            labs(x = "COPD - Control (log1p CP10K)", y = NULL) +
            theme_ms(7) +
            theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))
        q2 <- NULL
        if (!is.null(pw) && "mde_log1p_cp10k" %in% names(pw)) {
            m <- as.data.table(pw)[is.finite(mde_log1p_cp10k)]
            if (nrow(m)) {
                setorder(m, mde_log1p_cp10k)
                m[, gene := factor(gene, levels = m$gene)]
                q2 <- ggplot(m, aes(mde_log1p_cp10k, gene)) +
                    geom_col(width = 0.7, fill = "grey70", colour = "grey30",
                             linewidth = 0.15) +
                    labs(x = "Minimum detectable effect\n(pericytes, 80% power)",
                         y = NULL) +
                    theme_ms(7)
            }
        }
        sup2 <- if (is.null(q2)) tag(q1, "A") else
            wrap_plots(list(tag(q1, "A"), tag(q2, "B")), nrow = 1, widths = c(1.5, 1))
        save_fig("figureS_bm_copd", sup2, 10, 5)
        message("wrote figureS_bm_copd")
    }
}

sessioninfo::session_info()
