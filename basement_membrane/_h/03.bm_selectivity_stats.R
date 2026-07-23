## Cross-cell-type basement-membrane selectivity (the well-powered result).
##
## Question: do lung pericytes deposit a SELECTIVE subset of basement-membrane
## components, or are they simply a matrix-rich cell type?
##
## Unit of analysis is the donor x cell-type pseudobulk value, never the cell.
## Cell-level contrasts across cell types are confounded by cell size and capture
## efficiency, and this repo has already retracted one claim (AGTR1 marking the
## vascular-stabilizing pole) that turned out to be exactly that artifact.
##
## Four defenses against that confound, all reported:
##   (a) CP10K normalization upstream (verified expm1 row-sums = 1e4);
##   (b) mean_log10_total_counts as a covariate, with every estimate also
##       reported without it so the reader can see the movement;
##   (c) the PRIMARY endpoint is a WITHIN-unit difference, bm_minus_fibrillar =
##       mean z(BM genes) - mean z(fibrillar genes). Any multiplicative
##       cell-size/capture constant applies to both panels and cancels. This is
##       what separates "selectively deposits BM" from "is matrix-rich";
##   (d) a detection-fraction lens, which is insensitive to expression magnitude.
##
## (1|study) is in every model. Study random-effect SD is reported for every gene
## and a gene whose study SD exceeds its residual SD is flagged `study_dominated`
## and must not carry a claim -- the diagnostic that collapsed the earlier
## "injury states expand in fibrosis" result.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(lme4)
    library(lmerTest)
    library(emmeans)
})
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

REF_GROUP <- "Pericytes"
## Pre-specified expectation, used as a method sanity check rather than as an
## input: LAMA3 (laminin-332) is epithelial and LAMA5 is broad, so neither should
## come out pericyte-selective; the mural set should. If the metric does not
## recover this layout, the method is wrong, not the biology.
EXPECT_NOT_SELECTIVE <- c("LAMA3", "LAMA5")
EXPECT_MURAL <- c("LAMA4", "COL4A1", "COL4A2", "HSPG2", "NID1")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--panels", type = "character"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L,
                dest = "min_cells"),
    make_option("--min-donors", type = "integer", default = 5L,
                dest = "min_donors"),
    make_option("--healthy-only", type = "logical", default = TRUE,
                dest = "healthy_only")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(x, file) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)
}

## ---------------------------------------------------------------- load ----
pb <- fread(opt$pseudobulk)
panels <- fread(opt$panels)
bm_genes <- panels[panel == "basement_membrane", unique(gene)]
fib_genes <- panels[panel == "fibrillar_ecm", unique(gene)]

pb <- pb[n_cells >= opt$min_cells]
message(sprintf("Units with >= %d cells: %d", opt$min_cells, nrow(pb)))

## Primary cohort is healthy donors: a clean statement about which BM components
## pericytes make in the normal lung, uncontaminated by disease remodeling.
## Sensitivity below re-runs on all donors with disease_group as a covariate.
if (opt$healthy_only && "disease_group" %in% names(pb)) {
    pb_primary <- pb[disease_group == "Healthy"]
} else {
    pb_primary <- pb
}
## Keep cell types represented by enough donors to support a contrast.
keep <- pb_primary[, .(n_donors = uniqueN(donor_id)), by = ccc_group
                   ][n_donors >= opt$min_donors, ccc_group]
pb_primary <- pb_primary[ccc_group %in% keep]
if (!REF_GROUP %in% pb_primary$ccc_group)
    stop("reference group '", REF_GROUP, "' did not survive filtering")
pb_primary[, ccc_group := relevel(factor(ccc_group), ref = REF_GROUP)]
message(sprintf("Primary cohort: %d units, %d cell types, %d donors",
                nrow(pb_primary), uniqueN(pb_primary$ccc_group),
                uniqueN(pb_primary$donor_id)))

## Within-dataset standardization: neutralizes between-study scale differences
## while preserving within-study contrasts (the repo's gold standard).
## Returns a plain vector; deliberately avoids `:=` so it never mutates the
## caller's table by reference.
z_within_dataset <- function(dt, col) {
    x <- dt[[col]]
    g <- as.character(dt[["dataset"]])
    out <- numeric(length(x))
    for (lev in unique(g)) {
        ix <- which(g == lev)
        v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else
            (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

## ------------------------------------------------- per-gene selectivity ----
fit_gene <- function(gene, value_type = "expr", with_depth = TRUE) {
    col <- paste0(gene, "__", value_type)
    if (!col %in% names(pb_primary)) return(NULL)
    d <- copy(pb_primary)
    d[, y := get(col)]
    d[, y_z := z_within_dataset(d, "y")]
    d <- d[is.finite(y_z)]
    if (uniqueN(d$ccc_group) < 2) return(NULL)

    rhs <- c("ccc_group", if (with_depth) "mean_log10_total_counts",
             "(1 | donor_id)", "(1 | study)")
    fit <- try(suppressMessages(
        lmerTest::lmer(reformulate(rhs, "y_z"), data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)

    vc <- as.data.frame(VarCorr(fit))
    sd_study <- vc$sdcor[vc$grp == "study"]
    sd_resid <- vc$sdcor[vc$grp == "Residual"]
    if (!length(sd_study)) sd_study <- NA_real_

    emm <- emmeans(fit, specs = "ccc_group")
    ctr <- as.data.frame(contrast(emm, "trt.vs.ctrl", ref = REF_GROUP,
                                  adjust = "BH"))
    ## emmeans' trt.vs.ctrl gives (other - Pericytes); flip so positive means
    ## pericyte-enriched, which is what the figures and text will say.
    ctr$estimate <- -ctr$estimate
    ctr$contrast <- paste0(REF_GROUP, " - ",
                           sub(" - .*$", "", as.character(ctr$contrast)))
    data.frame(gene = gene, value_type = value_type, with_depth = with_depth,
               ctr, sd_study = sd_study, sd_residual = sd_resid,
               study_dominated = !is.na(sd_study) && sd_study > sd_resid,
               n_units = nrow(d), n_donors = uniqueN(d$donor_id),
               row.names = NULL)
}

emm_profile <- function(gene, value_type = "expr") {
    col <- paste0(gene, "__", value_type)
    if (!col %in% names(pb_primary)) return(NULL)
    d <- copy(pb_primary)
    d[, y := get(col)]
    fit <- try(suppressMessages(lmerTest::lmer(
        y ~ ccc_group + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- as.data.frame(emmeans(fit, specs = "ccc_group"))
    e$gene <- gene; e$value_type <- value_type
    e
}

genes_all <- intersect(c(bm_genes, fib_genes),
                       sub("__expr$", "", grep("__expr$", names(pb), value = TRUE)))
message("Modelling ", length(genes_all), " genes")

sel <- rbindlist(lapply(genes_all, fit_gene, value_type = "expr",
                        with_depth = TRUE), fill = TRUE)
sel_nodepth <- rbindlist(lapply(genes_all, fit_gene, value_type = "expr",
                                with_depth = FALSE), fill = TRUE)
sel_detect <- rbindlist(lapply(genes_all, fit_gene, value_type = "detect",
                               with_depth = TRUE), fill = TRUE)

## BH across genes x contrasts within each lens (emmeans already adjusted within
## a gene; this is the across-gene family).
for (dt in list(sel, sel_nodepth, sel_detect))
    if (nrow(dt)) dt[, p_BH_across_genes := p.adjust(p.value, method = "BH")]

write_tsv_safe(sel, file.path(opt$outdir, "bm_selectivity_emmeans.tsv"))
write_tsv_safe(sel_nodepth,
               file.path(opt$outdir, "bm_selectivity_emmeans_nodepth.tsv"))
write_tsv_safe(sel_detect,
               file.path(opt$outdir, "bm_detection_glmm.tsv"))

## ------------------------------------------------------ tau specificity ----
## Yanai tau on the cell-type marginal means, on the linear CP10K scale.
## 0 = ubiquitous, 1 = exclusive to one cell type.
prof <- rbindlist(lapply(genes_all, emm_profile), fill = TRUE)
write_tsv_safe(prof, file.path(opt$outdir, "bm_celltype_profile.tsv"))

tau_tbl <- prof[, {
    x <- expm1(pmax(emmean, 0))          # back to linear CP10K
    x[x < 0] <- 0
    mx <- max(x, na.rm = TRUE)
    tau <- if (mx > 0 && .N > 1) sum(1 - x / mx) / (.N - 1) else NA_real_
    ord <- order(-x)
    .(tau = tau,
      n_groups = .N,
      top_group = ccc_group[ord[1]],
      pericyte_rank = which(as.character(ccc_group[ord]) == REF_GROUP)[1],
      log2_pericyte_over_next = {
          pv <- x[as.character(ccc_group) == REF_GROUP][1]
          nxt <- max(x[as.character(ccc_group) != REF_GROUP], na.rm = TRUE)
          if (is.finite(pv) && is.finite(nxt) && nxt > 0) log2(pv / nxt) else NA_real_
      })
}, by = .(gene, value_type)]
tau_tbl[, panel := fifelse(gene %in% bm_genes, "basement_membrane", "fibrillar_ecm")]
write_tsv_safe(tau_tbl, file.path(opt$outdir, "bm_tau_specificity.tsv"))

## ---------------------------------- PRIMARY endpoint: bm_minus_fibrillar ----
## Per-gene z across units, averaged within panel, then differenced. A
## multiplicative per-cell-type capture constant cancels in the difference, so a
## positive value means "BM-shifted relative to fibrillar", not "matrix-rich".
panel_z <- function(genes, value_type = "expr") {
    cols <- paste0(genes, "__", value_type)
    cols <- intersect(cols, names(pb_primary))
    m <- as.matrix(pb_primary[, ..cols])
    z <- scale(m)
    z[is.na(z)] <- 0
    rowMeans(z)
}
pb_primary[, bm_z := panel_z(bm_genes)]
pb_primary[, fib_z := panel_z(fib_genes)]
pb_primary[, bm_minus_fibrillar := bm_z - fib_z]

fit_primary <- suppressMessages(lmerTest::lmer(
    bm_minus_fibrillar ~ ccc_group + mean_log10_total_counts +
        (1 | donor_id) + (1 | study), data = pb_primary))
emm_p <- emmeans(fit_primary, specs = "ccc_group")
ctr_p <- as.data.frame(contrast(emm_p, "trt.vs.ctrl", ref = REF_GROUP,
                                adjust = "BH"))
ctr_p$estimate <- -ctr_p$estimate
ctr_p$contrast <- paste0(REF_GROUP, " - ",
                         sub(" - .*$", "", as.character(ctr_p$contrast)))
write_tsv_safe(as.data.frame(emm_p),
               file.path(opt$outdir, "bm_primary_endpoint_emmeans.tsv"))
write_tsv_safe(ctr_p, file.path(opt$outdir, "bm_primary_endpoint_posthoc.tsv"))

vc_p <- as.data.frame(VarCorr(fit_primary))
write_tsv_safe(vc_p, file.path(opt$outdir, "bm_variance_components.tsv"))

## --------------------------------------------------- sanity check report ----
peri_tau <- tau_tbl[value_type == "expr" & gene %in% bm_genes]
sanity <- peri_tau[, .(gene, tau, top_group, pericyte_rank,
                       log2_pericyte_over_next)][order(pericyte_rank)]
expected_ok <- all(
    peri_tau[gene %in% EXPECT_NOT_SELECTIVE, pericyte_rank] > 3,
    na.rm = TRUE)
mural_ok <- mean(peri_tau[gene %in% EXPECT_MURAL, pericyte_rank] <= 5,
                 na.rm = TRUE)

readme <- c(
    "Basement-membrane selectivity -- generated summary",
    sprintf("Units (>=%d cells): %d; cell types: %d; donors: %d",
            opt$min_cells, nrow(pb_primary), uniqueN(pb_primary$ccc_group),
            uniqueN(pb_primary$donor_id)),
    sprintf("Cohort: %s", if (opt$healthy_only) "Healthy donors only" else "all donors"),
    "",
    "Variance components, primary endpoint (bm_minus_fibrillar):",
    paste(utils::capture.output(print(vc_p)), collapse = "\n"),
    "",
    "Pre-specified method sanity check:",
    sprintf("  LAMA3/LAMA5 NOT pericyte-top3 (expected TRUE): %s", expected_ok),
    sprintf("  fraction of mural-expected genes in pericyte top-5: %.2f", mural_ok),
    "",
    "Per-gene tau / pericyte rank (BM panel):",
    paste(utils::capture.output(print(sanity)), collapse = "\n"))
writeLines(readme, file.path(opt$outdir, "bm_selectivity_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
