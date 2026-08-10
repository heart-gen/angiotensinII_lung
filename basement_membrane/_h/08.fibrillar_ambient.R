## Do pericytes TRANSCRIBE fibrillar collagen, or are they carrying fibroblast soup?
##
## The cross-cell-type selectivity in 03 establishes that pericytes sit below
## every fibroblast class for COL1A1/COL3A1/COL5A*. It does not establish that
## the residual pericyte signal is real. Ambient RNA in lung scRNA-seq is
## dominated by whatever the abundant secretory cells release, and fibroblast
## collagen is a textbook contributor -- so "pericytes express COL1A2 in 60% of
## cells" is exactly the kind of statement that a reviewer will attribute to soup.
##
## Three independent controls, none of which requires an ambient-corrected layer
## (the `soupX` layer in pericyte_states.h5ad is bit-identical to `counts` --
## verified 2026-08-10, max abs diff 0 over all 11,680 cells x 55,329 genes -- so
## it carries no correction and cannot be used for this):
##
##   (1) AMBIENT FLOOR. Capillary endothelium and macrophages are in the same
##       dissociations, are at least as abundant as pericytes, and do not
##       transcribe fibrillar collagen. Whatever they show IS the soup level.
##       A pericyte gene indistinguishable from that floor is not evidence of
##       transcription; a gene an order of magnitude above it is.
##
##   (2) AMBIENT TRACER REGRESSION. Off-lineage transcripts (SFTPC, SFTPB,
##       SCGB1A1, SCGB3A2, PTPRC) measure each donor's soup burden directly.
##       Under the ambient hypothesis a donor whose pericytes carry more soup
##       must also carry more pericyte COL1A2, and a donor with more/higher-
##       expressing fibroblasts must carry more still. Both slopes are estimated.
##       A null slope with a high intercept is the transcription signature.
##
##   (3) STOICHIOMETRY. Collagen I is an obligate heterotrimer, alpha1(I)2
##       alpha2(I)1, so a cell that transcribes it carries COL1A1 and COL1A2
##       together. Ambient contamination transfers the SOURCE cell's ratio, so
##       soup-driven pericyte signal must reproduce the fibroblast COL1A1:COL1A2
##       relationship. A pericyte-specific departure from that ratio cannot be
##       produced by soup, and is the positive evidence that the pericyte
##       fibrillar signal is transcribed but incomplete.
##
## Unit of analysis is the donor x cell-type pseudobulk throughout, the same unit
## and the same (1|study) + depth model as 03. The stoichiometry endpoint is a
## WITHIN-unit difference of two genes, so per-cell-type capture constants cancel
## exactly as they do in bm_minus_fibrillar.
##
## Outputs (stats_data/):
##   fibrillar_ambient_floor.tsv       per-gene pericyte vs ambient-reference contrast
##   fibrillar_ambient_regression.tsv  soup-burden and fibroblast-burden slopes
##   fibrillar_stoichiometry_emmeans.tsv / _posthoc.tsv   COL1A1 - COL1A2 by cell type
##   fibrillar_stoichiometry_donor.tsv donor-level values behind the figure panel
##   fibrillar_ambient_README.txt      generated summary

suppressPackageStartupMessages({
    library(optparse)
    library(data.table)
    library(lme4)
    library(lmerTest)
    library(emmeans)
})
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

REF_GROUP <- "Pericytes"

## Cell types that do not transcribe fibrillar collagen and are co-dissociated
## with pericytes: the empirical soup floor. Endothelium and macrophages only --
## epithelium is deliberately excluded because AT2/club cells are the dominant
## SOURCE of lung soup, so their own tracer values are not a background estimate.
AMBIENT_REFERENCE <- c("EC general capillary", "EC aerocyte capillary",
                       "EC venous pulmonary", "EC venous systemic",
                       "Alveolar macrophages", "Interstitial macrophages")
FIBROBLAST_PAT <- "fibroblast|Myofibro"

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--panels", type = "character"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--min-donors", type = "integer", default = 5L, dest = "min_donors"),
    ## Matches 03: the identity claim is made in the normal lung.
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
core_genes <- panels[panel == "fibrillar_core", unique(gene)]
minor_genes <- panels[panel == "fibrillar_minor", unique(gene)]
tracer_genes <- panels[panel == "ambient_tracer", unique(gene)]
block_map <- unique(panels[, .(gene, block)])

have <- sub("__expr$", "", grep("__expr$", names(pb), value = TRUE))
fib_genes <- intersect(c(core_genes, minor_genes), have)
tracer_genes <- intersect(tracer_genes, have)
bm_genes <- intersect(bm_genes, have)
if (!length(tracer_genes))
    stop("no ambient tracer genes in the pseudobulk; re-run step_2 after ",
         "regenerating bm_panel_genes.tsv from the current bm_panels.py")

pb <- pb[n_cells >= opt$min_cells]
if (opt$healthy_only && "disease_group" %in% names(pb))
    pb <- pb[disease_group == "Healthy"]
keep <- pb[, .(n_donors = uniqueN(donor_id)), by = ccc_group
           ][n_donors >= opt$min_donors, ccc_group]
pb <- pb[ccc_group %in% keep]
if (!REF_GROUP %in% pb$ccc_group)
    stop("reference group '", REF_GROUP, "' did not survive filtering")

## `role` travels with every output so the figure never has to re-derive which
## cell types are the ambient reference -- one definition, here.
pb[, role := fifelse(ccc_group == REF_GROUP, "Pericytes",
             fifelse(grepl(FIBROBLAST_PAT, ccc_group, ignore.case = TRUE),
                     "Fibroblasts",
             fifelse(ccc_group %in% AMBIENT_REFERENCE, "Ambient reference",
                     "Other")))]
message(sprintf("Units: %d | cell types: %d | donors: %d | cohort: %s",
                nrow(pb), uniqueN(pb$ccc_group), uniqueN(pb$donor_id),
                if (opt$healthy_only) "Healthy only" else "all donors"))
print(unique(pb[, .(ccc_group, role)])[order(role, ccc_group)])
if (!any(pb$role == "Ambient reference"))
    stop("no ambient-reference cell type survived filtering; the floor control ",
         "cannot be computed and the fibrillar claim must not be made without it")

## ============================================================================
## (1) AMBIENT FLOOR: pericytes vs the co-dissociated non-collagen cell types
## ============================================================================
## Fitted on the raw log1p(CP10K) scale, NOT z-scored: the question is the size
## of the gap above background, which a within-gene z-score would destroy.
floor_gene <- function(gene, value_type = "expr") {
    col <- paste0(gene, "__", value_type)
    if (!col %in% names(pb)) return(NULL)
    d <- copy(pb)[role %in% c("Pericytes", "Ambient reference", "Fibroblasts")]
    d[, y := get(col)]
    d[, role := relevel(factor(role), ref = "Ambient reference")]
    d <- d[is.finite(y)]
    if (uniqueN(d$role) < 2) return(NULL)
    fit <- try(suppressMessages(lmerTest::lmer(
        y ~ role + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- emmeans(fit, specs = "role")
    ct <- as.data.frame(contrast(e, "trt.vs.ctrl", ref = "Ambient reference",
                                 adjust = "BH"))
    vc <- as.data.frame(VarCorr(fit))
    sd_study <- vc$sdcor[vc$grp == "study"]
    if (!length(sd_study)) sd_study <- NA_real_
    data.frame(gene = gene, value_type = value_type,
               block = block_map$block[match(gene, block_map$gene)],
               ct, sd_study = sd_study,
               sd_residual = vc$sdcor[vc$grp == "Residual"],
               n_units = nrow(d), n_donors = uniqueN(d$donor_id),
               row.names = NULL)
}

floor_genes <- unique(c(fib_genes, tracer_genes, bm_genes))
flr <- rbindlist(c(lapply(floor_genes, floor_gene, value_type = "expr"),
                   lapply(floor_genes, floor_gene, value_type = "detect")),
                 fill = TRUE)
if (nrow(flr)) {
    ## emmeans parenthesises levels containing spaces, so the contrast label is
    ## "(Pericytes) - (Ambient reference)" here and "Pericytes - Ambient
    ## reference" elsewhere depending on version. Carry an explicit flag rather
    ## than anchoring a regex to the start of the string.
    flr[, is_pericyte_contrast := grepl("Pericytes", contrast, fixed = TRUE)]
    flr[, study_dominated := is.finite(sd_study) & sd_study > sd_residual]
    ## BH within (lens x block): the fibrillar block carries the claim, the
    ## tracer block is the calibration, and they are separate questions.
    flr[, p_BH := p.adjust(p.value, method = "BH"), by = .(value_type, block)]
    write_tsv_safe(flr, file.path(opt$outdir, "fibrillar_ambient_floor.tsv"))
    message("Ambient-floor contrasts (expr lens, pericyte - ambient reference):")
    print(flr[value_type == "expr" & is_pericyte_contrast,
              .(gene, block, estimate = round(estimate, 3),
                p_BH = signif(p_BH, 2))][order(-estimate)])
}

## ============================================================================
## (2) AMBIENT TRACER + FIBROBLAST-BURDEN REGRESSION, donor level
## ============================================================================
## For each donor: the pericyte value, that donor's pericyte soup burden (mean
## z of the tracer panel measured IN pericytes), and that donor's fibroblast
## burden (fibroblast expression of the same gene, and the fibroblast share of
## the donor's niche). Ambient predicts positive slopes on all three.
niche_tot <- pb[, .(n_niche = sum(n_cells)), by = donor_id]
fib_share <- pb[role == "Fibroblasts", .(n_fib = sum(n_cells)), by = donor_id]
share <- merge(niche_tot, fib_share, by = "donor_id", all.x = TRUE)
share[is.na(n_fib), n_fib := 0]
share[, fib_frac := n_fib / n_niche]

peri <- copy(pb[ccc_group == REF_GROUP])
tracer_cols <- paste0(tracer_genes, "__expr")
tz <- scale(as.matrix(peri[, ..tracer_cols]))
tz[is.na(tz)] <- 0
peri[, soup_burden := rowMeans(tz)]

## Donor fibroblast expression, weighted by how many fibroblasts of each class
## that donor contributed -- an unweighted mean would let a 5-cell subpleural
## unit count as much as a 400-cell alveolar one in setting the soup source.
fib_cols <- paste0(fib_genes, "__expr")
fib_expr <- pb[role == "Fibroblasts", c(
    lapply(.SD, function(v) {
        ok <- is.finite(v) & is.finite(n_cells)
        if (!any(ok)) NA_real_ else stats::weighted.mean(v[ok], n_cells[ok])
    }), .(n_fib_cells = sum(n_cells))),
    by = donor_id, .SDcols = fib_cols]
setnames(fib_expr, fib_cols, paste0(fib_genes, "__fibexpr"))

dd <- merge(peri, share[, .(donor_id, fib_frac, n_niche)], by = "donor_id")
dd <- merge(dd, fib_expr, by = "donor_id", all.x = TRUE)
message(sprintf("Ambient regression: %d donors with pericytes (%d with fibroblasts)",
                nrow(dd), sum(is.finite(dd[[paste0(fib_genes[1], "__fibexpr")]]))))

amb_reg <- rbindlist(lapply(fib_genes, function(g) {
    y <- paste0(g, "__expr"); xf <- paste0(g, "__fibexpr")
    d <- copy(dd[is.finite(get(y)) & is.finite(get(xf)) & is.finite(fib_frac)])
    if (nrow(d) < 15 || uniqueN(d$study) < 2) return(NULL)
    d[, `:=`(y_ = get(y), fibexpr_ = get(xf))]
    fit <- try(suppressMessages(lmerTest::lmer(
        y_ ~ soup_burden + fibexpr_ + fib_frac + mean_log10_total_counts +
            (1 | study), data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    co <- as.data.frame(coef(summary(fit)))
    co$term <- rownames(co); rownames(co) <- NULL
    setDT(co)
    setnames(co, c("Estimate", "Std. Error", "Pr(>|t|)"),
             c("estimate", "SE", "p_value"), skip_absent = TRUE)
    co <- co[term %in% c("(Intercept)", "soup_burden", "fibexpr_", "fib_frac")]
    co[, `:=`(gene = g, block = block_map$block[match(g, block_map$gene)],
              n_donors = nrow(d),
              ## Partial R2 of the ambient block, likelihood-ratio style: how
              ## much of the between-donor variance the soup terms explain at all.
              pericyte_mean = mean(d$y_), fibroblast_mean = mean(d$fibexpr_))]
    co[]
}), fill = TRUE)
if (nrow(amb_reg)) {
    amb_reg[term != "(Intercept)",
            p_BH := p.adjust(p_value, method = "BH"), by = term]
    write_tsv_safe(amb_reg, file.path(opt$outdir, "fibrillar_ambient_regression.tsv"))
    message("Ambient-burden slopes (positive => consistent with soup):")
    print(amb_reg[term %in% c("soup_burden", "fib_frac"),
                  .(gene, term, estimate = round(estimate, 3),
                    p_BH = signif(p_BH, 2))][order(term, -estimate)])
}

## ============================================================================
## (3) STOICHIOMETRY: COL1A1 - COL1A2 within unit
## ============================================================================
## Collagen I is alpha1(I)2 alpha2(I)1. In a transcribing cell COL1A1 runs at or
## above COL1A2; ambient transfer preserves the source ratio. So the contrast of
## interest is pericytes vs FIBROBLASTS on this difference, and the ambient
## reference is shown alongside to bound what soup alone would produce.
STOICH <- c("COL1A1", "COL1A2")
if (!all(paste0(STOICH, "__expr") %in% names(pb))) {
    warning("COL1A1/COL1A2 not both present; skipping stoichiometry", call. = FALSE)
} else {
    d <- copy(pb)
    d[, stoich := COL1A1__expr - COL1A2__expr]
    d[, stoich_detect := COL1A1__detect - COL1A2__detect]
    ## Units where neither chain is seen carry no ratio information and would
    ## pile a spurious zero onto every cell type; drop them and say how many.
    n_all <- nrow(d)
    d <- d[COL1A1__expr > 0 | COL1A2__expr > 0]
    message(sprintf("Stoichiometry: %d/%d units carry >=1 collagen I chain",
                    nrow(d), n_all))
    d[, ccc_group := relevel(factor(ccc_group), ref = REF_GROUP)]

    fit_s <- suppressMessages(lmerTest::lmer(
        stoich ~ ccc_group + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d))
    e_s <- emmeans(fit_s, specs = "ccc_group")
    ct_s <- as.data.frame(contrast(e_s, "trt.vs.ctrl", ref = REF_GROUP,
                                   adjust = "BH"))
    ct_s$estimate <- -ct_s$estimate
    ct_s$contrast <- paste0(REF_GROUP, " - ",
                            sub(" - .*$", "", as.character(ct_s$contrast)))
    emm_s <- as.data.frame(e_s)
    emm_s$role <- d$role[match(emm_s$ccc_group, d$ccc_group)]
    emm_s$n_units <- as.integer(table(d$ccc_group)[as.character(emm_s$ccc_group)])
    write_tsv_safe(emm_s,
                   file.path(opt$outdir, "fibrillar_stoichiometry_emmeans.tsv"))
    write_tsv_safe(ct_s,
                   file.path(opt$outdir, "fibrillar_stoichiometry_posthoc.tsv"))

    vc_s <- as.data.frame(VarCorr(fit_s))
    write_tsv_safe(data.frame(endpoint = "COL1A1_minus_COL1A2", vc_s),
                   file.path(opt$outdir, "fibrillar_stoichiometry_varcomp.tsv"))

    ## Donor-level source data for the figure panel: one point per unit, so the
    ## panel shows the spread rather than only the model estimate.
    write_tsv_safe(
        d[, .(donor_id, ccc_group, role, study, n_cells,
              COL1A1 = COL1A1__expr, COL1A2 = COL1A2__expr,
              COL1A1_detect = COL1A1__detect, COL1A2_detect = COL1A2__detect,
              stoich, stoich_detect)],
        file.path(opt$outdir, "fibrillar_stoichiometry_donor.tsv"))

    message("COL1A1 - COL1A2 by role (emmeans):")
    print(as.data.table(emm_s)[, .(ccc_group, role, emmean = round(emmean, 3),
                                   n_units)][order(emmean)])
}

## --------------------------------------------------------------- summary ----
peri_floor <- if (exists("flr") && nrow(flr))
    flr[value_type == "expr" & is_pericyte_contrast &
        block %in% c("fibrillar_core", "fibrillar_minor", "ambient_tracer"),
        .(gene, block, estimate = round(estimate, 3), p_BH = signif(p_BH, 2))
        ][order(block, -estimate)] else NULL

readme <- c(
    "Fibrillar-collagen ambient controls -- generated summary",
    sprintf("Cohort: %s | units: %d | cell types: %d | donors: %d",
            if (opt$healthy_only) "Healthy donors only" else "all donors",
            nrow(pb), uniqueN(pb$ccc_group), uniqueN(pb$donor_id)),
    sprintf("Ambient reference groups: %s",
            paste(sort(unique(pb[role == "Ambient reference", ccc_group])),
                  collapse = ", ")),
    "",
    "NOTE: the soupX layer in pericyte_states.h5ad is bit-identical to counts",
    "(verified 2026-08-10). It carries no ambient correction and is NOT used.",
    "",
    "(1) Pericyte minus ambient-floor, log1p(CP10K), expr lens:",
    if (!is.null(peri_floor))
        paste(utils::capture.output(print(peri_floor)), collapse = "\n") else
        "  not computed",
    "",
    "(2) Ambient-burden slopes (soup_burden = pericyte off-lineage tracer z;",
    "    fibexpr_ = same gene in that donor's fibroblasts; fib_frac = fibroblast",
    "    share of the niche). Positive and significant => consistent with soup.",
    if (exists("amb_reg") && nrow(amb_reg))
        paste(utils::capture.output(print(
            amb_reg[term != "(Intercept)",
                    .(gene, term, estimate = round(estimate, 3),
                      p_BH = signif(p_BH, 2))][order(gene, term)])),
            collapse = "\n") else "  not computed",
    "",
    "(3) COL1A1 - COL1A2 within unit, by cell type:",
    if (exists("emm_s"))
        paste(utils::capture.output(print(
            as.data.table(emm_s)[order(emmean),
                                 .(ccc_group, role, emmean = round(emmean, 3),
                                   n_units)])), collapse = "\n") else
        "  not computed")
writeLines(readme, file.path(opt$outdir, "fibrillar_ambient_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
