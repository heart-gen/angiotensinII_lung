## Supplementary Tables S7 (donor- and cell-level correlations), S8 (cross-cell-type
## basement-membrane expression and pericyte selectivity) and S9 (cluster-level BM
## estimates, orthogonal models, panel correlations, subpanel robustness).
##
## Three gaps are filled here rather than assembled:
##   S7  state_gate_axis_correlation.tsv stores a bare pearson_r with no n, CI or P.
##       n and a Fisher-z CI are computed from the cell-level score matrix, and the
##       pseudotime correlations get the BH column they never had.
##   S8C the Figure 3A colour (z within gene), dot size (detection fraction) and
##       per-cell-type denominators are computed INLINE in
##       figures/_h/basement_membrane_figure.R:32-66 and written nowhere. The same
##       arithmetic is reproduced here so the table is genuinely the figure's source.
##   S9  the posthoc files carry p.value only, so the "12 contrasts significant at
##       BH < 0.05" claim has no BH column behind it. BH is added.
##
## Outputs: tsv/tableS07.tsv, tableS08A.tsv, tableS08B.tsv, tableS08C.tsv,
##          tableS09A.tsv .. tableS09D.tsv

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

BM  <- function(...) P("basement_membrane", "_m", ...)
BMS <- function(f) BM("stats_data", f)
PS  <- function(f) P("pericyte_states", "_m", f)
PSS <- function(f) P("pericyte_states", "_m", "stats_data", f)

## Gene order comes from basement_membrane/_m/bm_panel_genes.tsv (written by
## bm_panels.py), NOT from a literal list here. It was a literal until 2026-09-01,
## when the BM panel grew from 13 to 20 genes and this table silently kept
## reporting the old 13 -- the same duplicated-panel failure the BM module exists
## to prevent. Falls back to the frozen 13 only if the file is missing.
pan0 <- read_src(BM("bm_panel_genes.tsv"))
GENE_ORDER <- if (!is.null(pan0) && "block_index" %in% names(pan0))
    unique(pan0[block == "basement_membrane"][order(block_index), gene]) else
    c("COL4A1", "COL4A2", "COL18A1",
      "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1",
      "NID1", "NID2", "HSPG2", "AGRN")

## =========================================================================
## S7 -- donor- and cell-level correlations
## =========================================================================
bits <- list()

pt <- read_src(PS("pseudotime_trend_correlations.tsv"))
if (!is.null(pt)) {
    ## BH within level: cell-level and donor-level are separate families.
    pt[, p_BH := bh(p_value), by = level]
    if ("partial_p" %in% names(pt)) pt[, partial_p_BH := bh(partial_p), by = level]
    bits$pt <- pt[, .(analysis = "pseudotime trend", level, feature,
                      estimate = spearman_rho, method = "Spearman",
                      p_value, p_BH, n,
                      partial_estimate = if ("partial_rho" %in% names(pt)) partial_rho else NA_real_,
                      partial_p = if ("partial_p" %in% names(pt)) partial_p else NA_real_)]
}

ac <- read_src(PSS("acta2_control_correlations.tsv"))
if (!is.null(ac)) bits$acta <- ac[, .(
    analysis = "ACTA2 control", level = "donor x program",
    feature = paste(x, "vs", y), estimate = pearson_r, method = "Pearson",
    p_value = pearson_p, p_BH = pearson_p_BH, n)]

ap <- read_src(PSS("AGTR1_program_correlations.tsv"))
if (!is.null(ap)) bits$agtr1 <- ap[, .(
    analysis = "AGTR1 vs program score", level = "donor", feature,
    estimate = spearman_rho, method = "Spearman", p_value, p_BH, n)]

bf <- read_src(BMS("bm_vs_fibrillar_corr.tsv"))
if (!is.null(bf)) bits$bmfib <- bf[, .(
    analysis = "BM vs fibrillar panel", level = "donor x cluster pseudobulk",
    feature = comparison, estimate = pearson_r, method = "Pearson",
    ci_lower = pearson_ci_lo, ci_upper = pearson_ci_hi,
    p_value = pearson_p, n = n_units,
    spearman_rho, spearman_p)]

bc <- read_src(BMS("bm_continuum_summary.tsv"))
if (!is.null(bc)) bits$bmcont <- bc[, .(
    analysis = "BM score vs continuum (per-donor rho, Wilcoxon)", level = "donor",
    feature = metric, estimate = median_rho, method = "Spearman (median of donors)",
    p_value = p_wilcox, p_BH = p_BH, n = n_donors)]

sub <- read_src(P("localization", "pericyte_analysis", "_m", "stats_data",
                  "correlation_by_subcluster.tsv"))
if (!is.null(sub)) bits$sub <- sub[, .(
    analysis = "AGTR1 vs airspace within subcluster", level = "donor",
    feature = paste("subcluster", leiden_pericytes), estimate = spearman_rho,
    method = "Spearman", p_value, n = n_donors)]

## Cell-level BM axis correlations: pearson_r only in the source. Recover n from
## the cell-level score table and attach a Fisher-z interval, so the r = 0.052 /
## 0.068 / 0.209 values quoted in BM_SUMMARY.md finally have an interval and a P.
ax <- read_src(BM("state_gate_axis_correlation.tsv"))
if (!is.null(ax)) {
    ## n is the cell-level score table these correlations were computed on
    ## (basement_membrane/_h/01.state_gate.py operates on bm_metadata).
    n_cells <- NA_integer_
    sg <- BM("bm_metadata.tsv.gz")
    if (file.exists(sg)) n_cells <- nrow(fread(sg, select = 1L))
    ci <- fisher_ci(ax$pearson_r, n_cells)
    ## Two-sided P for a Pearson r at this n.
    ## stats:: qualified because `pt` is bound to the pseudotime table above.
    tstat <- ax$pearson_r * sqrt((n_cells - 2) / pmax(1 - ax$pearson_r^2, 1e-12))
    pv <- 2 * stats::pt(abs(tstat), n_cells - 2, lower.tail = FALSE)
    bits$axis <- data.table(
        analysis = "BM vs contrast panel (cell level)", level = "cell",
        feature = paste(ax$panel_a, "vs", ax$panel_b),
        estimate = ax$pearson_r, method = "Pearson",
        ci_lower = ci$lo, ci_upper = ci$hi, p_value = pv, n = n_cells)
}

if (length(bits)) {
    s7 <- rbindlist(bits, fill = TRUE)
    if ("p_BH" %in% names(s7)) s7[is.na(p_BH) & !is.na(p_value), p_BH := bh(p_value)]
    s7[, p_formatted := fmt_p(p_value)]
    setcolorder(s7, intersect(c("analysis", "level", "feature", "method", "estimate",
                                "ci_lower", "ci_upper", "p_value", "p_BH", "n"),
                              names(s7)))
    write_part(s7, "07",
        "Donor- and cell-level correlations across the pericyte analyses",
        supports = "Figures 1-3; Figures S1, S4, S7",
        sources = c("pericyte_states/_m/pseudotime_trend_correlations.tsv",
                    "pericyte_states/_m/stats_data/{acta2_control,AGTR1_program}_correlations.tsv",
                    "basement_membrane/_m/stats_data/{bm_vs_fibrillar_corr,bm_continuum_summary}.tsv",
                    "basement_membrane/_m/state_gate_axis_correlation.tsv",
                    "localization/pericyte_analysis/_m/stats_data/correlation_by_subcluster.tsv"),
        notes = paste("BH ADDED HERE for the pseudotime correlations (adjusted",
                      "within level) and for any block whose source stored only a",
                      "nominal P. The cell-level BM-axis rows carry a Fisher-z 95%",
                      "CI and a P computed for this table: the source file stored",
                      "a bare Pearson r with no n, interval or test. P values",
                      "shown as '<2.2e-16' underflowed in the source."))
}

## =========================================================================
## S8 -- cross-cell-type BM expression and pericyte selectivity
## =========================================================================
## ---- A: panel definition ------------------------------------------------
pan <- pan0
det_per <- read_src(BM("bm_gene_detection_pericytes.tsv"))
## Structural subclass, covering all 20 genes of the 2026-09-01 panel. Anything
## not listed falls back to the sub-panel membership recorded in the panel table.
SUBCLASS <- c(COL4A1 = "collagen IV", COL4A2 = "collagen IV",
              COL15A1 = "collagen XV (multiplexin)",
              COL18A1 = "collagen XVIII (multiplexin)",
              LAMA1 = "laminin alpha", LAMA2 = "laminin alpha",
              LAMA3 = "laminin alpha", LAMA4 = "laminin alpha",
              LAMA5 = "laminin alpha",
              LAMB1 = "laminin beta", LAMB2 = "laminin beta",
              LAMB3 = "laminin beta", LAMB4 = "laminin beta",
              LAMC1 = "laminin gamma", LAMC2 = "laminin gamma",
              LAMC3 = "laminin gamma",
              NID1 = "linker/proteoglycan", NID2 = "linker/proteoglycan",
              HSPG2 = "linker/proteoglycan", AGRN = "linker/proteoglycan")
if (!is.null(pan)) {
    a <- pan[panel == "basement_membrane"]
    a[, structural_subclass := SUBCLASS[gene]]
    a[, n_panels := pan[, .N, by = gene][match(a$gene, gene), N]]
    a[, other_panels := sapply(a$gene, function(g)
        paste(setdiff(pan[gene == g, panel], "basement_membrane"), collapse = "; "))]
    if (!is.null(det_per)) a <- merge(a, det_per, by = "gene", all.x = TRUE)
    a[, gene := factor(gene, levels = GENE_ORDER)]; setorder(a, gene)
    write_part(a, "08A",
        "Basement-membrane panel definition and pericyte detection",
        supports = "Figure 3A",
        sources = c("basement_membrane/_m/bm_panel_genes.tsv",
                    "basement_membrane/_m/bm_gene_detection_pericytes.tsv"),
        notes = paste("Panel provenance: basement_membrane/_h/bm_panels.py.",
                      "Detection statistics are within pericytes only."))
}

## ---- B: per-gene pericyte selectivity -----------------------------------
tau <- read_src(BMS("bm_tau_specificity.tsv"))
sel <- read_src(BMS("bm_selectivity_emmeans.tsv"))
if (!is.null(tau) || !is.null(sel)) {
    b <- list()
    if (!is.null(tau)) b$tau <- tau[, block := "tau specificity"]
    if (!is.null(sel)) b$sel <- sel[gene %in% GENE_ORDER][, block := "pericyte vs cell-type contrast"]
    s8b <- rbindlist(b, fill = TRUE)
    if ("p.value" %in% names(s8b)) s8b[, p_formatted := fmt_p(p.value)]
    write_part(s8b, "08B",
        "Per-gene pericyte selectivity of basement-membrane components",
        supports = "Figure 3C",
        sources = c("basement_membrane/_m/stats_data/bm_tau_specificity.tsv",
                    "basement_membrane/_m/stats_data/bm_selectivity_emmeans.tsv"),
        notes = paste("Contrast rows are Pericytes minus each other population.",
                      "`p_BH_across_genes` is the source script's BH adjustment.",
                      "Orthogonal replicates without the depth covariate and with a",
                      "detection GLMM are in bm_selectivity_emmeans_nodepth.tsv and",
                      "bm_detection_glmm.tsv."))
}

## ---- C: full BM-panel x 22-cell-type source data for Figure 3A/3B -------
prof <- read_src(BMS("bm_celltype_profile.tsv"))
pbk  <- read_src(BM("bm_pseudobulk_celltype.tsv.gz"))
if (!is.null(prof)) {
    d <- prof[gene %in% GENE_ORDER & value_type == "expr"]
    ## Identical arithmetic to basement_membrane_figure.R:37 -- z within gene.
    d[, z_within_gene := (emmean - mean(emmean)) / (sd(emmean) + 1e-9), by = gene]
    if (!is.null(pbk)) {
        dcols <- intersect(paste0(GENE_ORDER, "__detect"), names(pbk))
        det <- pbk[n_cells >= 5, lapply(.SD, mean, na.rm = TRUE),
                   by = ccc_group, .SDcols = dcols]
        det <- melt(det, id.vars = "ccc_group", variable.name = "gene",
                    value.name = "detect_frac")
        det[, gene := sub("__detect$", "", gene)]
        d <- merge(d, det, by = c("ccc_group", "gene"), all.x = TRUE)
        ## Per-cell-type denominators under EXACTLY the filter the model used
        ## (03.bm_selectivity_stats.R:70-83): n_cells >= 5, HEALTHY DONORS ONLY,
        ## then cell types with >= 5 donors. Omitting the healthy filter inflates
        ## this to 4,376 profiles / 353 donors -- which are the RAS module's
        ## all-donor numbers, not the BM cohort's 2,329 / 220.
        pb_h <- pbk[n_cells >= 5]
        if ("disease_group" %in% names(pb_h)) pb_h <- pb_h[disease_group == "Healthy"]
        den <- pb_h[, .(n_donors = uniqueN(donor_id), n_profiles = .N,
                        n_cells_total = sum(n_cells)), by = ccc_group]
        den <- den[n_donors >= 5]
        cat(sprintf("BM cohort: %d profiles, %d donors, %d cell types\n",
                    sum(den$n_profiles), uniqueN(pb_h[ccc_group %in% den$ccc_group,
                                                      donor_id]), nrow(den)))
        d <- merge(d, den, by = "ccc_group", all.x = TRUE)
    }
    d[, within_gene_rank := frank(-emmean, ties.method = "min"), by = gene]
    d[, gene := factor(gene, levels = GENE_ORDER)]
    setorder(d, gene, within_gene_rank)
    write_part(d, "08C",
        "Figure 3A/3B source data: matrix-panel expression across 22 lung cell types",
        supports = "Figure 3A; Figure 3B",
        sources = c("basement_membrane/_m/stats_data/bm_celltype_profile.tsv",
                    "basement_membrane/_m/bm_pseudobulk_celltype.tsv.gz"),
        notes = paste("Model: expr_z ~ ccc_group + mean_log10_total_counts +",
                      "(1|donor_id) + (1|study), healthy donors only.",
                      "`z_within_gene` is the panel colour and `detect_frac` the",
                      "dot size; both, and the per-cell-type denominators, were",
                      "computed inline in the figure script and are written to a",
                      "file for the first time here. Cohort totals: 2,329 eligible",
                      "donor x cell-type profiles, 220 healthy donors,",
                      "22 populations."))
}

## =========================================================================
## S9 -- cluster-level BM estimates, orthogonal models, robustness
## =========================================================================
## The BM posthoc files were written with `pairs(emm, adjust = "BH")`
## (04.bm_state_stats.R:106, 141), so their `p.value` is ALREADY BH-adjusted
## within each score. Re-adjusting would apply BH twice and silently drop the
## significant-contrast count from 12 to 11. Relabel instead of adjusting.
label_adjusted <- function(dt) {
    if (is.null(dt) || !"p.value" %in% names(dt)) return(dt)
    setnames(dt, "p.value", "p_BH_within_score")
    dt[, adjustment := "BH within score (applied by emmeans::pairs)"]
    dt[, p_formatted := fmt_p(p_BH_within_score)][]
}

cl_emm <- read_src(BMS("bm_by_cluster_emmeans.tsv"))
cl_ph  <- read_src(BMS("bm_by_cluster_posthoc.tsv"))
## Per-cluster denominators. The headline "214 donor x cluster units, 95 donors"
## exists only as prose in bm_pericyte_axis_README.txt and was never broken down
## per cluster. Rebuild the same pseudobulk the model used: donor x cluster with
## >= 5 cells (04.bm_state_stats.R:74-76, --min-cells default 5).
bm_meta <- read_src(BM("bm_metadata.tsv.gz"),
                    select = c("donor_id", "pericyte_state"))
if (!is.null(cl_emm)) {
    if (!is.null(bm_meta)) {
        pbu <- bm_meta[, .(n_cells = .N), by = .(donor_id, pericyte_state)][n_cells >= 5]
        den <- pbu[, .(n_units = .N, n_donors = uniqueN(donor_id),
                       n_cells_total = sum(n_cells)), by = pericyte_state]
        den[, pericyte_state := as.character(pericyte_state)]
        cl_emm[, pericyte_state := as.character(pericyte_state)]
        cl_emm <- merge(cl_emm, den, by = "pericyte_state", all.x = TRUE)
        cat(sprintf("donor x cluster units (>=5 cells): %d; donors: %d\n",
                    nrow(pbu), uniqueN(pbu$donor_id)))
    }
    write_part(cl_emm, "09A",
        "Cluster-level basement-membrane and contrast-panel scores (marginal means)",
        supports = "Figure S17A",
        sources = "basement_membrane/_m/stats_data/bm_by_cluster_emmeans.tsv",
        notes = paste0(uniqueN(cl_emm$score), " panel scores across the ",
                       uniqueN(cl_emm$pericyte_state),
                       " stable pericyte clusters. Grouping is `pericyte_state` ",
                       "(Leiden on X_pca_harmony), not `state_program`: the latter ",
                       "is assigned by an argmax that includes the BM panel, so ",
                       "testing a BM score against it would be circular."))
}
if (!is.null(cl_ph)) {
    cl_ph <- label_adjusted(cl_ph)
    n_sig <- cl_ph[score == "basement_membrane_score_z", sum(p_BH_within_score < 0.05)]
    write_part(cl_ph, "09B",
        "Cluster-level basement-membrane scores: pairwise contrasts",
        supports = "Figure S17A",
        sources = "basement_membrane/_m/stats_data/bm_by_cluster_posthoc.tsv",
        notes = paste0("`p_BH_within_score` is the source file's `p.value`, ",
                       "renamed: emmeans::pairs was called with adjust = 'BH', so ",
                       "it is already adjusted across the 15 contrasts within each ",
                       "score and must NOT be adjusted again. The basement-membrane ",
                       "score has ", n_sig, " contrasts significant at BH < 0.05, ",
                       "which is the figure quoted in BM_SUMMARY.md."))
}

orth_e <- read_src(BMS("bm_orthogonalized_emmeans.tsv"))
orth_p <- read_src(BMS("bm_orthogonalized_posthoc.tsv"))
if (!is.null(orth_e) || !is.null(orth_p)) {
    o <- list()
    if (!is.null(orth_e)) o$emm <- orth_e[, block := "marginal means"]
    if (!is.null(orth_p)) o$ph  <- label_adjusted(orth_p)[, block := "pairwise contrasts"]
    write_part(rbindlist(o, fill = TRUE), "09C",
        "Basement-membrane score adjusted for the fibrillar-ECM panel (orthogonal model)",
        supports = "Figure S17A",
        sources = c("basement_membrane/_m/stats_data/bm_orthogonalized_emmeans.tsv",
                    "basement_membrane/_m/stats_data/bm_orthogonalized_posthoc.tsv"),
        notes = paste("Tests whether the cluster-level BM signal survives",
                      "adjustment for the correlated fibrillar matrix program.",
                      "Contrast P values are already BH-adjusted by the source",
                      "script (emmeans::pairs adjust = 'BH')."))
}

## ---- E: matrix-vs-predictor associations (Figure S17B/C source data) -----
## Added 2026-09-01 with the AGTR1 / TGF-beta association analysis. Both matrix
## categories and their difference are outcomes; only the difference is invariant
## to a per-unit multiplicative capture constant, so it is the row that licenses a
## directional statement.
assoc <- rbindlist(Filter(Negate(is.null), list(
    read_src(BMS("bm_vs_agtr1_models.tsv")),
    read_src(BMS("bm_vs_tgfb_models.tsv")))), fill = TRUE)
within <- rbindlist(Filter(Negate(is.null), list(
    read_src(BMS("bm_vs_agtr1_within_donor.tsv")),
    read_src(BMS("bm_vs_tgfb_within_donor.tsv")))), fill = TRUE)
if (nrow(assoc) || nrow(within)) {
    e_bits <- list()
    if (nrow(assoc)) e_bits$pb <- assoc[, block := "donor x cluster pseudobulk (mixed model)"]
    if (nrow(within)) e_bits$wd <- within[, block := "within-donor cell-level Spearman"]
    write_part(rbindlist(e_bits, fill = TRUE), "09E",
        "Associations of the two matrix programs with AGTR1 and TGF-beta signalling",
        supports = "Figure S17B; Figure S17C",
        sources = c("basement_membrane/_m/stats_data/bm_vs_agtr1_{models,within_donor}.tsv",
                    "basement_membrane/_m/stats_data/bm_vs_tgfb_{models,within_donor}.tsv"),
        notes = paste("Outcomes are the basement-membrane score, the fibrillar-",
                      "collagen score, their difference, and the BM score",
                      "orthogonalized against the frozen fibrillar_ecm panel.",
                      "The scVI-denoised lens is the AGTR1 readout; the raw and",
                      "detection lenses are carried as the dropout-sensitive",
                      "comparison and must not be read as the answer.",
                      "The TGF-beta response panel shares no gene with either",
                      "matrix panel or with any pericyte_states program",
                      "(asserted in basement_membrane/_h/bm_panels.py), so the",
                      "association is not arithmetic."))
}

rob <- read_src(BM("state_gate_robustness.tsv"))
xt  <- read_src(BM("state_gate_crosstab.tsv"))
re  <- read_src(BM("state_gate_relenrich.tsv"))
gs  <- read_src(BM("state_gate_summary.tsv"))
d_bits <- list()
if (!is.null(rob)) d_bits$rob <- rob[, block := "subpanel variant"]
if (!is.null(xt))  d_bits$xt  <- xt[,  block := "5-panel vs 6-panel assignment"]
if (!is.null(re))  d_bits$re  <- re[,  block := "relative enrichment"]
if (!is.null(gs))  d_bits$gs  <- gs[,  block := "gate summary"]
if (length(d_bits))
    write_part(rbindlist(d_bits, fill = TRUE), "09D",
        "Subpanel robustness of the basement-membrane program assignment",
        supports = "Figure 3H",
        sources = c("basement_membrane/_m/state_gate_robustness.tsv",
                    "basement_membrane/_m/state_gate_{crosstab,relenrich,summary}.tsv"),
        notes = paste("Variants drop each BM subpanel in turn (collagen IV only,",
                      "laminin only, linker only) and remove the shared COL4A1",
                      "from the fibroblast-like contrast panel. `flipped` marks a",
                      "cluster whose dominant program changes under the variant."))

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
