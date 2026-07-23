## Basement-membrane program and AGT ligand-axis figures (Circ Research revision).
##
## Main: which basement-membrane components pericytes selectively deposit, how the
## BM axis separates from fibrillar ECM, and the state-model consequence.
## Supplements: the local RAS landscape, and the COPD contrast with its power bound.
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P <- function(...) file.path(ROOT, ...)
OUT <- P("figures", "mechanism"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- shared visual language (matches manuscript_mechanism_figure.R) ------
OKABE <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9",
           "#F0E442", "#999999")

theme_ms <- function(base = 8) {
    theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              strip.background = element_blank())
}
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}
tag <- function(p, lab) p + labs(tag = lab) +
    theme(plot.tag = element_text(face = "bold", size = 10))

BM_M <- P("basement_membrane", "_m")
AGT_M <- P("agt_axis", "_m")
rd <- function(f) if (file.exists(f)) fread(f) else NULL

## Gene order: structural grouping, so the reader sees collagen IV / laminins /
## linkers as blocks rather than an alphabetical jumble.
GENE_ORDER <- c("COL4A1", "COL4A2", "COL18A1",
                "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1",
                "NID1", "NID2", "HSPG2", "AGRN")

## ======================= A: BM gene x cell-type dot plot ==================
prof <- rd(file.path(BM_M, "stats_data", "bm_celltype_profile.tsv"))
pbk <- rd(file.path(BM_M, "bm_pseudobulk_celltype.tsv.gz"))
pA <- NULL
if (!is.null(prof)) {
    d <- prof[gene %in% GENE_ORDER & value_type == "expr"]
    ## z within gene so cell types are comparable across genes of very different
    ## absolute abundance; dot size carries the detection fraction, which lives in
    ## the pseudobulk table rather than the emmeans profile.
    d[, z := (emmean - mean(emmean)) / (sd(emmean) + 1e-9), by = gene]
    if (!is.null(pbk)) {
        dcols <- paste0(GENE_ORDER, "__detect")
        dcols <- intersect(dcols, names(pbk))
        det <- pbk[n_cells >= 5, lapply(.SD, mean, na.rm = TRUE),
                   by = ccc_group, .SDcols = dcols]
        det <- melt(det, id.vars = "ccc_group", variable.name = "gene",
                    value.name = "detect")
        det[, gene := sub("__detect$", "", gene)]
        d <- merge(d, det, by = c("ccc_group", "gene"), all.x = TRUE)
    } else {
        d[, detect := 0.5]
    }
    d[, gene := factor(gene, levels = rev(GENE_ORDER))]
    ord <- d[, .(m = mean(z)), by = ccc_group][order(-m), ccc_group]
    d[, ccc_group := factor(ccc_group, levels = ord)]
    pA <- ggplot(d, aes(ccc_group, gene)) +
        geom_point(aes(size = detect, fill = z), shape = 21, stroke = 0.15,
                   colour = "grey30") +
        scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
                             midpoint = 0, name = "Expression\n(z within gene)") +
        scale_size_continuous(range = c(0.3, 3.4), name = "Detected\nfraction") +
        labs(x = NULL, y = NULL) +
        theme_ms(7) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "right", legend.key.size = unit(0.3, "cm"))
}

## ============ B: selectivity forest, pericyte vs next-highest =============
tau <- rd(file.path(BM_M, "stats_data", "bm_tau_specificity.tsv"))
pB <- NULL
if (!is.null(tau)) {
    d <- tau[value_type == "expr" & gene %in% GENE_ORDER]
    setorder(d, log2_pericyte_over_next)
    d[, gene := factor(gene, levels = d$gene)]
    d[, enriched := log2_pericyte_over_next > 0]
    pB <- ggplot(d, aes(log2_pericyte_over_next, gene, fill = enriched)) +
        geom_col(width = 0.7, colour = "grey30", linewidth = 0.15) +
        geom_vline(xintercept = 0, linewidth = 0.3) +
        geom_text(aes(label = sprintf("#%d", pericyte_rank),
                      hjust = ifelse(enriched, -0.25, 1.25)), size = 2.1) +
        scale_fill_manual(values = c(`TRUE` = OKABE[1], `FALSE` = "grey75"),
                          guide = "none") +
        scale_x_continuous(expand = expansion(mult = 0.18)) +
        labs(x = expression(log[2]~"(pericyte / next-highest cell type)"), y = NULL) +
        theme_ms(7)
}

## ======= C: primary endpoint, BM-vs-fibrillar shift by cell type ==========
pe <- rd(file.path(BM_M, "stats_data", "bm_primary_endpoint_emmeans.tsv"))
pC <- NULL
if (!is.null(pe)) {
    d <- as.data.table(pe)
    setorder(d, emmean)
    d[, ccc_group := factor(ccc_group, levels = d$ccc_group)]
    d[, grp := fifelse(ccc_group == "Pericytes", "Pericytes",
               fifelse(grepl("fibroblast|Myofibro", ccc_group, ignore.case = TRUE),
                       "Fibroblasts",
               fifelse(grepl("^EC |Lymphatic", ccc_group), "Endothelium", "Other")))]
    pC <- ggplot(d, aes(emmean, ccc_group, colour = grp)) +
        geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3, colour = "grey50") +
        geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0, orientation = "y",
                       linewidth = 0.4) +
        geom_point(size = 1.5) +
        scale_colour_manual(values = c(Pericytes = OKABE[4], Fibroblasts = OKABE[2],
                                       Endothelium = OKABE[1], Other = "grey60"),
                            name = NULL) +
        labs(x = "BM - fibrillar score (panel z difference)", y = NULL) +
        theme_ms(7) +
        theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))
}

## ============= D: the state-model consequence (the gate result) ===========
gate <- rd(file.path(BM_M, "state_gate_relenrich.tsv"))
ct <- rd(file.path(BM_M, "state_gate_crosstab.tsv"))
pD <- NULL
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
    pD <- ggplot(d, aes(program, pericyte_state, fill = relenrich)) +
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

panels <- Filter(Negate(is.null), list(tag(pA, "A"), tag(pB, "B"),
                                       tag(pC, "C"), tag(pD, "D")))
if (length(panels)) {
    fig <- wrap_plots(panels, ncol = 2, widths = c(1.15, 1))
    save_fig("figure_basement_membrane", fig, 11, 8.5)
    message("wrote figure_basement_membrane")
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
