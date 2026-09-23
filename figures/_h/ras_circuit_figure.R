## Figure 5 (figure_ras_circuit) and Figure S18 (figureS_ras_circuit_robustness):
## the pericyte as the AT1R-responsive node of a distributed lung RAS.
##
## Reads ras_circuit/_m only. Every annotation is DERIVED from those tables; a
## missing input stops the script (read_req) rather than silently dropping a
## panel -- the failure that cost assemble_mechanism_figures.R a panel on
## 2026-09-10. No in-panel titles; interpretation lives in the README legends.
## ASCII only inside text (cairo_pdf mangles rho/beta glyphs).
##
## Layout: A spans the top row (the schematic needs the width), then B|C, D|E and
## F. Panel A's coordinates are RECOMPUTED here, because the stored x/y are a logical
## ordering and collide when drawn at page width; the drawing order of the nodes is a
## layout decision and lives here, while the graph itself and the estimates that
## annotate its edges come from 07 and 11 respectively.

suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})
source("../_h/_fig_common.R")

RC <- function(...) P("ras_circuit", "_m", ...)
SD <- function(f) RC("stats_data", f)
read_req <- function(path) {
    if (!file.exists(path))
        stop("missing input: ", path, "\n  Run ras_circuit steps 1-4 first.", call. = FALSE)
    fread(path)
}
theme_ms <- function(base = 8) .theme_ms(base = base, legend = NULL)
COMP_LAB <- c(AT1 = "AT1", AT2 = "AT2", EC_aerocyte = "Aerocyte EC", EC_gcap = "General cap. EC")
CLASS_COL <- c(substrate = OKABE[2], processing = OKABE[1], receptor = OKABE[4],
               counter_regulatory = "grey75", response = OKABE[5], response_branch = OKABE[5],
               matrix = OKABE[6], neighbor = OKABE[3], latent = "white")

## ================================ A: DAG ==========================================
## Edge grammar, stated in the legend: an ARROW head is a promoting / mass-flow link,
## a BAR head is an antagonistic one (prior biology, from edge_sign). A SOLID edge
## carries a supported donor-level estimate, a DOTTED edge one that was tested and came
## back null, and a THIN GREY edge is a biochemical conversion with no estimate at all
## (the peptides are latent). Numbers come from 11.dag_annotations.R, never from here.
nodes <- read_req(SD("ras_dag_nodes.tsv"))
edges <- read_req(SD("ras_dag_edges.tsv"))[drawn == TRUE]
scal  <- read_req(SD("ras_dag_scalars.tsv"))
dann  <- read_req(SD("ras_dag_annotations.tsv"))
dnote <- read_req(SD("ras_dag_notes.tsv"))

## Compact labels: gene, cell type abbreviated, detection.
CELL_SHORT <- c("Vascular smooth muscle" = "VSMC", "Alveolar fibroblasts" = "Alv fib",
                "Adventitial fibroblasts" = "Adv fib", "EC aerocyte capillary" = "aCap EC",
                "EC general capillary" = "gCap EC", "Alveolar macrophages" = "Alv mac",
                "Mast cells" = "Mast", "Pericytes" = "Pericyte")
NODE_LAB <- c(Peri_AT1R_response = "AT1R\nresponse", contractile = "Contractile",
              inflammatory = "Inflammatory", activated_migratory = "Activated /\nmigratory",
              matrix = "Matrix\nremodelling",
              BM = "Basement\nmembrane", FIB = "Fibrillar\nECM",
              out_EC = "Capillary\nEC", out_FIB = "Fibroblasts", out_EPI = "AT1 / AT2",
              AngI = "Ang I", AngII = "Ang II", Ang17 = "Ang 1-7")
nodes[, lab := fifelse(
    !is.na(gene) & !is.na(cell_type),
    sprintf("%s\n%s\n%.2f", gene, CELL_SHORT[cell_type], detect),
    fifelse(node %in% names(NODE_LAB), NODE_LAB[node], label))]
## Coordinates are set node by node rather than spread evenly inside a tier: the
## outgoing receivers have to sit level with the node that feeds them, or the
## neighbour edges cross the matrix branch.
## Columns widen towards the right, where the edges are short and every one of them
## carries a label.
XY <- data.table(matrix(byrow = TRUE, ncol = 3, c(
    "VSMC_AGT", 1, 1.7,     "AlvFib_AGT", 1, 0.2,     "AdvFib_AGT", 1, -1.3,
    "AngI", 2.1, 0.2,
    "ECaero_ACE", 3.2, 2.2, "ECgcap_ACE", 3.2, 1.0,   "AlvMac_ACE", 3.2, -0.2,
    "Mast_CMA1", 3.2, -1.4, "Mast_CTSG", 3.2, -2.6,
    "AngII", 4.4, 0.6,      "Peri_ACE2", 4.4, -2.9,
    "Peri_AGTR1", 5.8, 1.9, "Ang17", 5.8, -2.9,
    "Peri_AT1R_response", 7.0, 0.4,                   "Peri_MAS1", 7.0, -2.9,
    "contractile", 8.6, 2.7, "inflammatory", 8.6, 1.7, "activated_migratory", 8.6, 0.7,
    "matrix", 8.6, -0.7,
    "BM", 10.4, 0.0,        "FIB", 10.4, -1.6,
    "out_EPI", 12.4, 2.6,   "out_EC", 12.4, 0.0,      "out_FIB", 12.4, -1.6)))
setnames(XY, c("node", "xx", "yy")); XY[, c("xx", "yy") := lapply(.SD, as.numeric),
                                        .SDcols = c("xx", "yy")]
nodes <- merge(nodes, XY, by = "node")
stopifnot(!anyNA(nodes$xx))
eg <- merge(merge(edges, XY, by.x = "from", by.y = "node"),
            XY, by.x = "to", by.y = "node", suffixes = c("", "end"))
eg <- merge(eg, dann[, .(from, to, stat_label = label, stat_kind, measured_sign, show_in_panel)],
            by = c("from", "to"), all.x = TRUE)
## Solid where an estimate is supported, dotted where it was tested and null, thin grey
## where the link is a conversion no model can reach.
eg[, evidence := fifelse(stat_kind == "NicheNet z", "inferred",
                  fifelse(is.na(measured_sign), "prior",
                   fifelse(measured_sign == "null", "tested_null", "supported")))]
eg[is.na(evidence), evidence := "prior"]
eg[, head_shape := fifelse(edge_sign == "negative", "bar", "arrow")]
## Shorten every edge at both ends so the head lands outside the label box; the
## curved ones are trimmed less, because the arc runs longer than the chord.
shrink <- function(a, b, f) a + (b - a) * f
eg[, `:=`(dx = xxend - xx, dy = yyend - yy)]
eg[, len := pmax(sqrt(dx^2 + dy^2), 1e-6)]
eg[, `:=`(f0 = 0.14, f1 = 0.84)]
## The outgoing edges run between short and long, so trim them by a fixed distance
## rather than a fraction: a fraction either buries the head in the receiver box or
## leaves it hanging in white space.
eg[evidence == "inferred", `:=`(f0 = 0.30 / len, f1 = 1 - 0.55 / len)]
eg[, `:=`(x0 = shrink(xx, xxend, f0), y0 = shrink(yy, yyend, f0),
          x1 = shrink(xx, xxend, f1), y1 = shrink(yy, yyend, f1))]
## Labels sit beside their edge, not on it: offset perpendicular to the segment, and
## for a curved edge far enough out to clear the apex of the arc.
eg[, off := fifelse(evidence == "inferred", 0.25 + 0.15 * len, 0.30)]
eg[, lpos := fifelse(evidence == "inferred", 0.62, 0.50)]
eg[, `:=`(lx = shrink(xx, xxend, lpos) - off * dy / len,
          ly = shrink(yy, yyend, lpos) + off * dx / len)]
elab <- eg[show_in_panel == TRUE & nzchar(stat_label)]
elab[, stat_label := gsub("; ", "\n", stat_label)]
TIER_CAP <- data.table(
    xx = XY[match(c("VSMC_AGT", "ECaero_ACE", "Peri_AGTR1", "Peri_AT1R_response",
                    "contractile", "BM", "out_EC"), node), xx], yy = 3.5,
    cap = c("AGT source", "processing", "receptor", "response", "programmes",
            "matrix", "neighbours"))
max_ren <- scal[quantity == "max_REN_detect", value]
notes <- dnote[show_in_panel == TRUE]
arrow_h <- arrow(length = unit(1.6, "mm"), type = "closed")
bar_h   <- arrow(length = unit(1.8, "mm"), angle = 90, type = "open")
pA <- ggplot() +
    geom_segment(data = eg[evidence == "prior" & head_shape == "arrow"],
                 aes(x0, y0, xend = x1, yend = y1), colour = "grey72", linewidth = 0.22,
                 arrow = arrow_h) +
    geom_segment(data = eg[evidence == "supported" & head_shape == "arrow"],
                 aes(x0, y0, xend = x1, yend = y1), colour = "grey35", linewidth = 0.42,
                 arrow = arrow_h) +
    geom_segment(data = eg[evidence == "tested_null" & head_shape == "arrow"],
                 aes(x0, y0, xend = x1, yend = y1), colour = "grey55", linewidth = 0.3,
                 linetype = "dotted", arrow = arrow_h) +
    ## Ligand-receptor inference (NicheNet), not a donor-level estimate: S18E carries
    ## the donor-level test, which is null for every receiver but adventitial fibroblasts.
    ## Curved, so the long edge to the epithelium arcs over the programme column
    ## instead of cutting through it.
    geom_curve(data = eg[evidence == "inferred" & len > 3],
               aes(x0, y0, xend = x1, yend = y1), colour = OKABE[3], linewidth = 0.35,
               linetype = "dashed", curvature = -0.30, arrow = arrow_h) +
    ## The short outgoing edges take a gentler arc; the same curvature over a short
    ## chord turns the arrow head sideways.
    geom_curve(data = eg[evidence == "inferred" & len <= 3],
               aes(x0, y0, xend = x1, yend = y1), colour = OKABE[3], linewidth = 0.35,
               linetype = "dashed", curvature = -0.10, arrow = arrow_h) +
    ## Antagonistic links: bar head, and here the estimate behind it is null.
    geom_segment(data = eg[head_shape == "bar"],
                 aes(x0, y0, xend = x1, yend = y1), colour = "grey35", linewidth = 0.38,
                 linetype = "dotted", arrow = bar_h) +
    geom_label(data = elab, aes(lx, ly, label = stat_label), size = 1.5, colour = "grey15",
               fill = "white", label.size = 0, lineheight = 0.9,
               label.padding = unit(0.3, "mm"), alpha = 0.85) +
    geom_label(data = nodes[node_class != "latent"],
               aes(xx, yy, label = lab, fill = node_class), size = 1.75, label.size = 0.12,
               lineheight = 0.9, label.padding = unit(0.6, "mm"), alpha = 0.9) +
    geom_text(data = nodes[node_class == "latent"], aes(xx, yy, label = lab),
              size = 1.9, fontface = "italic", colour = "grey20") +
    geom_text(data = TIER_CAP, aes(xx, yy, label = cap), size = 1.75, fontface = "italic",
              colour = "grey45") +
    annotate("text", x = 0.55, y = -3.55 - 0.42 * seq_len(nrow(notes)), hjust = 0,
             size = 1.6, colour = "grey30", label = notes$text) +
    annotate("text", x = 0.55, y = -3.55, hjust = 0, size = 1.6, colour = "grey30",
             label = sprintf("REN max detection %s; no cell type carries >1 step", max_ren)) +
    scale_fill_manual(values = CLASS_COL, guide = "none") +
    scale_x_continuous(expand = expansion(add = 0.62)) +
    scale_y_continuous(expand = expansion(add = 0.30)) +
    theme_void(base_size = 8)

## ======================= B: niche-affinity decomposition =========================
sl <- read_req(SD("niche_affinity_agtr1_models.tsv"))[arm == "primary"]
gl <- read_req(SD("niche_affinity_global.tsv"))[arm == "primary"]
sl[, comp := factor(COMP_LAB[compartment], levels = rev(COMP_LAB))]
pB <- ggplot(sl, aes(slope, comp)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = null_mean - 2 * null_sd, xmax = null_mean + 2 * null_sd),
                   height = 0, linewidth = 2.2, colour = "grey88") +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.18, linewidth = 0.4) +
    geom_point(size = 1.7, colour = OKABE[1]) +
    labs(x = "Affinity slope per SD AGTR1", y = NULL,
         caption = sprintf("interaction LRT %s, n = %d donors", fmt_p(gl$lrt_p), gl$n_donors)) +
    theme_ms() + theme(plot.caption = element_text(size = 5.5, hjust = 0))

## ============================ C: RAS network =====================================
ed <- read_req(SD("ras_network_edges.tsv"))[arm == "primary" &
                                              family %in% c("allowed", "forbidden")]
PRETTY <- c(AGT_source_index = "AGT source", processing_index = "ACE/chymase",
            Peri_AGTR1 = "Pericyte AGTR1", Peri_ACE2 = "Pericyte ACE2",
            Peri_AT1R_response = "AT1R response", BM_minus_FIB = "BM - fibrillar",
            EC_readout = "EC barrier", activated_migratory = "Activated/migratory",
            contractile = "Contractile", inflammatory = "Inflammatory",
            VSMC_AGT = "VSMC AGT", Mast_CMA1 = "Mast CMA1")
nm <- function(x) ifelse(is.na(PRETTY[x]), x, PRETTY[x])
ed[, lab := sprintf("%s -> %s", nm(from), nm(to))]
setorder(ed, family, partial_rho)
ed[, lab := factor(lab, levels = unique(lab))]
ed[, fam := factor(fifelse(family == "allowed", "DAG edges", "Expected null"),
                   levels = c("DAG edges", "Expected null"))]
pC <- ggplot(ed, aes(partial_rho, lab)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = boot_lo, xmax = boot_hi), height = 0.18, linewidth = 0.35) +
    geom_point(aes(colour = p_BH < 0.05), size = 1.5) +
    scale_colour_manual(values = c(`TRUE` = OKABE[4], `FALSE` = "grey55"),
                        labels = c(`TRUE` = "BH < 0.05", `FALSE` = "n.s."), name = NULL) +
    facet_grid(fam ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Partial Spearman (95% bootstrap)", y = NULL,
         caption = sprintf("n = %d-%d donors per edge",
                           min(ed$n_donors, na.rm = TRUE), max(ed$n_donors, na.rm = TRUE))) +
    theme_ms() + theme(legend.position = "bottom", strip.text = element_text(size = 5.5),
                       axis.text.y = element_text(size = 5.5),
                       legend.text = element_text(size = 6),
                       legend.key.size = unit(3, "mm"),
                       plot.caption = element_text(size = 5.5, hjust = 0))

## =========================== D: AT1R response ====================================
arb <- read_req(SD("at1r_vs_agtr1_count.tsv"))[model %in% c("NB GLMM", "Poisson+OLRE") &
                                                 spec == "primary"]
cs  <- read_req(SD("at1r_continuum_summary.tsv"))[kind == "raw"]
keep_fam <- c("signature", "legacy_curated", "progeny_prespecified", "tf")
nice <- function(s) sub("^progeny_", "PROGENy ", sub("^tf_", "TF ", sub(
    "at1r_response_score", "AngII signature", sub("at1r_legacy_curated_score", "Curated AT1R", s))))
dA <- arb[family %in% keep_fam, .(score, family, metric = "AGTR1 count model",
                                  est = estimate, lo = estimate - 1.96 * SE,
                                  hi = estimate + 1.96 * SE, null_mean, null_sd)]
dB <- cs[family %in% keep_fam, .(score, family, metric = "Pseudotime rho",
                                 est = median_rho, lo = q25, hi = q75,
                                 null_mean = if ("null_mean" %in% names(cs)) null_mean else NA_real_,
                                 null_sd = if ("null_sd" %in% names(cs)) null_sd else NA_real_)]
dd <- rbind(dA, dB, fill = TRUE)
ord <- unique(dd[order(match(family, keep_fam), score), score])
dd[, lab := factor(nice(score), levels = rev(nice(ord)))]
sig_emp <- arb[score == "at1r_response_score", p_emp][1]
pD <- ggplot(dd, aes(est, lab)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbarh(data = dd[is.finite(null_mean)],
                   aes(xmin = null_mean - 2 * null_sd, xmax = null_mean + 2 * null_sd),
                   height = 0, linewidth = 2, colour = "grey88") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.32) +
    geom_point(aes(colour = family == "signature"), size = 1.4) +
    scale_colour_manual(values = c(`TRUE` = OKABE[4], `FALSE` = "grey35"), guide = "none") +
    facet_wrap(~ metric, scales = "free_x") +
    labs(x = NULL, y = NULL,
         caption = sprintf("signature vs %d matched panels, p_emp = %s",
                           arb[score == "at1r_response_score", n_null][1],
                           formatC(sig_emp, digits = 2, format = "g"))) +
    theme_ms() + theme(axis.text.y = element_text(size = 5.5),
                       strip.text = element_text(size = 6),
                       plot.caption = element_text(size = 5.5, hjust = 0))

## ========================== E: matrix consequence ================================
mx <- read_req(SD("matrix_vs_at1r.tsv"))
SPEC_LAB <- c(primary = "BM - fibrillar collagen", `_tracer_adj` = "+ ambient tracer",
              `_depth_spline` = "+ depth spline", `_fibrillar_ecm_contrast` = "BM - fibrillar ECM",
              `_bm_alone` = "BM alone (n.c.)", `_fib_alone` = "Fibrillar alone (n.c.)")
mx[, lab := factor(SPEC_LAB[spec], levels = rev(SPEC_LAB))]
pE <- ggplot(mx, aes(estimate, lab)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbarh(data = mx[is.finite(null_mean)],
                   aes(xmin = null_mean - 2 * null_sd, xmax = null_mean + 2 * null_sd),
                   height = 0, linewidth = 2, colour = "grey88") +
    geom_errorbarh(aes(xmin = estimate - 1.96 * SE, xmax = estimate + 1.96 * SE),
                   height = 0.18, linewidth = 0.35) +
    geom_point(aes(shape = claimable), size = 1.7, colour = OKABE[6]) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    labs(x = "Matrix outcome per SD of signature", y = NULL,
         caption = sprintf("primary vs matched null, p_emp = %s; n.c. = not claimable",
                           formatC(mx[spec == "primary", p_emp], digits = 2, format = "g"))) +
    theme_ms() + theme(axis.text.y = element_text(size = 6),
                       plot.caption = element_text(size = 5.5, hjust = 0))

## =========================== F: communication ====================================
ds <- read_req(RC("liana_direction_summary.tsv"))[disease_group == "Healthy"]
## Keep the recombined AT2 row only, so each partner appears once.
ds <- ds[!partner %in% c("AT2_AGTR2det", "AT2_AGTR2undet")]
ds[, partner := sub(" \\(either stratum\\)", "", partner)]
ds[, n := fifelse(direction == "outgoing", n_supported, -n_supported)]
ord_p <- ds[direction == "outgoing"][order(n_supported), partner]
ds[, partner := factor(partner, levels = ord_p)]
f1 <- ggplot(ds, aes(n, partner, fill = direction)) +
    geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(outgoing = OKABE[3], incoming = OKABE[1]), name = NULL) +
    scale_x_continuous(labels = abs) +
    labs(x = "Supported edges (Healthy)", y = NULL) +
    theme_ms() + theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5),
                       legend.text = element_text(size = 6),
                       legend.key.size = unit(3, "mm"))
nn <- read_req(RC("nichenet_outgoing", "outgoing_nichenet_summary.tsv"))[status == "run"]
SHORT <- c("EC aerocyte capillary" = "Aerocyte EC", "EC general capillary" = "Gen. cap. EC",
           "AT1" = "AT1", "AT2" = "AT2", "Alveolar fibroblasts" = "Alv. fibroblast",
           "Adventitial fibroblasts" = "Adv. fibroblast")
tp <- rbindlist(lapply(nn$target, function(t) {
    f <- RC("nichenet_outgoing", paste0("specificity_", gsub("[^A-Za-z0-9]+", "_", t), ".tsv"))
    x <- head(read_req(f)[!not_a_ligand & !substrate_not_ligand][order(rank_by_z)], 3)
    x[, target := t]
}))
## One small panel per receiver, each on its OWN x scale: the permutation z runs from
## ~1 in the capillary EC programme to ~75 in the epithelial one, so a shared axis
## flattens four of the six receivers into a single tick. Six stacked rows in a
## half-width column also collided; a 2 x 3 grid of three ligands each does not.
tp[, tgt := factor(SHORT[target], levels = SHORT)]
tp[, row := factor(paste(test_ligand, target), levels = paste(test_ligand, target)[order(tgt, z)])]
f2 <- ggplot(tp, aes(z, row)) +
    geom_segment(aes(x = 0, xend = z, yend = row), colour = "grey75", linewidth = 0.3) +
    geom_point(size = 1.3, colour = OKABE[3]) +
    scale_y_discrete(labels = function(s) sub(" .*$", "", s),
                     expand = expansion(add = 0.62)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18)), n.breaks = 4) +
    facet_wrap(~ tgt, ncol = 2, scales = "free") +
    labs(x = "Pericyte ligand -> receiver programme (permutation z)", y = NULL) +
    theme_ms(base = 7) +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 5.5),
          strip.text = element_text(size = 6), panel.spacing.x = unit(3, "mm"),
          panel.spacing.y = unit(1.6, "mm"))
pF <- f1 + f2 + plot_layout(widths = c(1, 1.45))

## ============================ assemble Figure 5 ==================================
## pF is two plots; wrap it so patchwork tags it once as panel F, not F and G.
fig <- wrap_elements(full = pA) / (pB | pC) / (pD | pE) / wrap_elements(full = pF) +
    plot_layout(heights = c(1.4, 0.95, 0.9, 1.45)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figure_ras_circuit", fig, 7.2, 9.8)

## ========================== Figure S18: robustness ===============================
sa <- read_req(SD("niche_affinity_agtr1_models.tsv"))
sa[, comp := factor(COMP_LAB[compartment], levels = rev(COMP_LAB))]
sA <- ggplot(sa, aes(slope, comp, colour = arm)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_pointrange(aes(xmin = lower.CL, xmax = upper.CL),
                    position = position_dodge(width = 0.7), size = 0.12, linewidth = 0.3) +
    scale_colour_manual(values = OKABE, name = NULL) +
    labs(x = "Affinity slope on AGTR1, by arm", y = NULL) +
    theme_ms() + theme(legend.position = "bottom", legend.text = element_text(size = 5),
                       legend.key.size = unit(2.5, "mm"))
nb <- read_req(SD("at1r_vs_agtr1_count_null.tsv"))
obs_b <- arb[score == "at1r_response_score", estimate][1]
sB <- ggplot(nb, aes(estimate)) + geom_histogram(bins = 40, fill = "grey75") +
    geom_vline(xintercept = obs_b, colour = OKABE[4], linewidth = 0.6) +
    labs(x = "Count-model beta, null panels", y = "Panels") + theme_ms()
nmx <- read_req(SD("matrix_vs_at1r_null.tsv"))[outcome == "bm_minus_fib_z"]
sC <- ggplot(nmx, aes(estimate)) + geom_histogram(bins = 40, fill = "grey75") +
    geom_vline(xintercept = mx[spec == "primary", estimate], colour = OKABE[6], linewidth = 0.6) +
    labs(x = "BM - fibrillar beta, null panels", y = "Panels") + theme_ms()
lo <- read_req(SD("ras_network_lodo.tsv"))[family == "allowed"]
lo[, lab := sprintf("%s -> %s", nm(from), nm(to))]
edp <- ed[family == "allowed"][, .(lab = sprintf("%s -> %s", nm(from), nm(to)), partial_rho)]
sD <- ggplot(lo, aes(partial_rho, lab)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_point(size = 0.5, alpha = 0.5, colour = "grey45") +
    geom_point(data = edp, aes(partial_rho, lab), colour = OKABE[4], size = 1.5) +
    labs(x = "Partial rho, each dataset left out", y = NULL) +
    theme_ms() + theme(axis.text.y = element_text(size = 5))
dva <- read_req(SD("outgoing_donor_validation.tsv"))
dva[, tgt := factor(SHORT[target], levels = rev(SHORT))]
sE <- ggplot(dva, aes(estimate, tgt, colour = spec)) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_errorbarh(data = dva[spec == "primary" & is.finite(null_mean)],
                   aes(xmin = null_mean - 2 * null_sd, xmax = null_mean + 2 * null_sd),
                   height = 0, linewidth = 2, colour = "grey88") +
    geom_pointrange(aes(xmin = estimate - 1.96 * SE, xmax = estimate + 1.96 * SE),
                    position = position_dodge(width = 0.6), size = 0.12, linewidth = 0.3) +
    scale_colour_manual(values = OKABE, name = NULL) +
    labs(x = "Receiver programme per SD of ligand composite", y = NULL) +
    theme_ms() + theme(legend.position = "bottom", legend.text = element_text(size = 5),
                       axis.text.y = element_text(size = 5.5),
                       legend.key.size = unit(2.5, "mm"))
au <- read_req(SD("at1r_signature_audit.tsv"))
au[, status := fifelse(kept %in% c("True", "TRUE"), "kept",
                fifelse(pruned_reason == "absent", "absent from lung object",
                        sub("\\..*$", "", pruned_reason)))]
sF <- ggplot(au[, .N, by = .(direction, status)], aes(N, status, fill = direction)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(up = OKABE[4], down = OKABE[1]), name = NULL) +
    labs(x = "Signature genes", y = NULL) +
    theme_ms() + theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5),
                       legend.text = element_text(size = 6),
                       legend.key.size = unit(3, "mm"))
figS <- (sA | sB | sC) / (sD | sE | sF) +
    plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_ras_circuit_robustness", figS, 9.0, 7.0)
message("figure_ras_circuit + figureS_ras_circuit_robustness written")
