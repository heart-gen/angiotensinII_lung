## Is basement-membrane deposition a distinct pericyte axis, or the fibrillar axis?
##
## CIRCULARITY WARNING, and the reason this script groups by `pericyte_state`
## rather than `state_program`: after the gate escalated, `state_program` is
## assigned by a relative-enrichment argmax that INCLUDES the BM panel. Testing
## "BM score differs by state_program" would therefore be testing the BM score
## against a label derived from the BM score. `pericyte_state` (the stable Leiden
## clusters) is panel-independent -- it comes from clustering on X_pca_harmony --
## so it is the honest grouping variable. state_program is reported alongside as
## a descriptive crosswalk only, explicitly flagged.
##
## Unit of analysis is the donor x cluster pseudobulk, dataset-standardized, with
## (1|study). Depth enters as a covariate everywhere. The BM panel is far better
## detected than AGTR1 ever was (COL4A2 81%, COL4A1 79%, LAMB1 60%, NID1 58% of
## pericytes) so this is not the sparse-gene regime that produced the retracted
## AGTR1 result -- but the same defenses are applied regardless.
##
## Also tests AGTR1 against BM directly (--denoise), at both the pseudobulk and
## the within-donor level. Figure panel F pairs a positive AGTR1-vs-pseudotime
## correlation with a negative BM-vs-pseudotime one; that is two correlations
## through a shared mediator, not an AGTR1-BM antagonism, and the manuscript may
## not claim the latter without this test. Outputs: bm_vs_agtr1_corr.tsv,
## bm_vs_agtr1_models.tsv, bm_vs_agtr1_within_donor.tsv.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(data.table)
    library(lme4)
    library(lmerTest)
    library(emmeans)
})
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--bm-meta", type = "character", dest = "bm_meta"),
    make_option("--state-meta", type = "character", dest = "state_meta"),
    make_option("--continuum", type = "character", default = NA_character_),
    ## scVI-denoised AGTR1, same table and model that
    ## pericyte_states/_h/03.agtr1_lenses.R uses as its robust lens.
    make_option("--denoise", type = "character", default = NA_character_),
    make_option("--den-model", type = "character", default = "Pericyte-only-trained",
                dest = "den_model"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(x, file) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)
}

bm <- fread(opt$bm_meta)
st <- fread(opt$state_meta)
setnames(st, 1, "index")
setnames(bm, 1, "index")
drop <- intersect(setdiff(names(bm), "index"), names(st))
d <- merge(st, bm[, .SD, .SDcols = setdiff(names(bm), drop)], by = "index")
message(sprintf("Merged %d cells", nrow(d)))

d[, pericyte_state := factor(pericyte_state)]
if ("dataset" %in% names(d)) d[, dataset := as.character(dataset)]

z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

## ------------------------------------------- donor x cluster pseudobulk ----
score_cols <- c("basement_membrane_score", "bm_collagen_iv_score",
                "bm_laminin_score", "bm_linker_score", "fibrillar_ecm_score",
                "fibrillar_core_score", "fibrillar_minor_score",
                "ambient_tracer_score",
                "fibroblast_like_noCOL4A1_score", "fibroblast_like_score")
score_cols <- intersect(score_cols, names(d))

pb <- d[, c(lapply(.SD, mean, na.rm = TRUE),
            .(n_cells = .N,
              study = first(study), dataset = first(dataset),
              mean_log10_counts = mean(log10_total_counts, na.rm = TRUE))),
        by = .(donor_id, pericyte_state),
        .SDcols = score_cols]
pb <- pb[n_cells >= opt$min_cells]
message(sprintf("Donor x cluster units (>=%d cells): %d across %d donors",
                opt$min_cells, nrow(pb), uniqueN(pb$donor_id)))

for (cl in score_cols)
    pb[[paste0(cl, "_z")]] <- z_within_dataset(pb[[cl]], pb$dataset)

## ------------------------------------------------- variance-component gate ----
## If between-study SD exceeds residual SD, the cluster claim cannot stand --
## the exact diagnostic that collapsed the earlier fibrosis claim.
vc_fit <- suppressMessages(lmer(
    basement_membrane_score_z ~ 1 + (1 | study) + (1 | dataset) + (1 | donor_id),
    data = pb))
vc <- as.data.frame(VarCorr(vc_fit))
write_tsv_safe(vc, file.path(opt$outdir, "bm_variance_components_pericyte.tsv"))
sd_study <- vc$sdcor[vc$grp == "study"]
sd_resid <- vc$sdcor[vc$grp == "Residual"]
study_dominated <- length(sd_study) && sd_study > sd_resid
message(sprintf("Variance: study SD %.3f vs residual SD %.3f -> %s",
                ifelse(length(sd_study), sd_study, NA), sd_resid,
                ifelse(study_dominated, "STUDY DOMINATED", "ok")))

## ------------------------------------------- BM across Leiden clusters ----
by_cluster <- function(col) {
    if (!col %in% names(pb)) return(NULL)
    f <- suppressMessages(lmer(
        reformulate(c("pericyte_state", "mean_log10_counts", "(1 | study)",
                      "(1 | donor_id)"), col), data = pb))
    e <- emmeans(f, specs = "pericyte_state")
    list(emm = data.frame(score = col, as.data.frame(e)),
         post = data.frame(score = col, as.data.frame(pairs(e, adjust = "BH"))))
}
res <- Filter(Negate(is.null), lapply(paste0(score_cols, "_z"), by_cluster))
write_tsv_safe(rbindlist(lapply(res, `[[`, "emm"), fill = TRUE),
               file.path(opt$outdir, "bm_by_cluster_emmeans.tsv"))
write_tsv_safe(rbindlist(lapply(res, `[[`, "post"), fill = TRUE),
               file.path(opt$outdir, "bm_by_cluster_posthoc.tsv"))

## ------------------------------------ BM vs fibrillar: are they one axis? ----
corr_tbl <- rbindlist(lapply(
    c("fibrillar_ecm_score", "fibroblast_like_noCOL4A1_score",
      "fibroblast_like_score"),
    function(other) {
        if (!other %in% names(pb)) return(NULL)
        ct <- cor.test(pb$basement_membrane_score, pb[[other]], method = "pearson")
        cs <- suppressWarnings(cor.test(pb$basement_membrane_score, pb[[other]],
                                        method = "spearman"))
        data.table(comparison = paste0("basement_membrane vs ", other),
                   pearson_r = unname(ct$estimate),
                   pearson_ci_lo = ct$conf.int[1], pearson_ci_hi = ct$conf.int[2],
                   pearson_p = ct$p.value, spearman_rho = unname(cs$estimate),
                   spearman_p = cs$p.value, n_units = nrow(pb))
    }), fill = TRUE)
write_tsv_safe(corr_tbl, file.path(opt$outdir, "bm_vs_fibrillar_corr.tsv"))
message("BM vs fibrillar correlations (donor x cluster pseudobulk):")
print(corr_tbl)

## Orthogonalization: does the cluster structure of BM survive adjustment for
## the fibrillar axis? If it vanishes, BM is not a separable axis -- report that.
orth <- suppressMessages(lmer(
    basement_membrane_score_z ~ pericyte_state + fibrillar_ecm_score_z +
        mean_log10_counts + (1 | study) + (1 | donor_id), data = pb))
e_orth <- emmeans(orth, specs = "pericyte_state")
write_tsv_safe(data.frame(as.data.frame(e_orth), model = "adjusted_for_fibrillar"),
               file.path(opt$outdir, "bm_orthogonalized_emmeans.tsv"))
write_tsv_safe(as.data.frame(pairs(e_orth, adjust = "BH")),
               file.path(opt$outdir, "bm_orthogonalized_posthoc.tsv"))

## Residual BM axis, carried forward to the continuum test.
pb[, bm_resid := residuals(suppressMessages(lmer(
    basement_membrane_score_z ~ fibrillar_ecm_score_z + mean_log10_counts +
        (1 | study) + (1 | donor_id), data = pb)))]

## ------------------------------------------- direct AGTR1 vs BM association ----
## Panel F of figure_pericyte_layer shows AGTR1 rising and BM falling along
## pseudotime. Those are two SEPARATE correlations sharing a mediator: they do not
## establish that AGTR1 and BM oppose each other, and until now nothing in this
## module or in pericyte_states/ tested the pair directly. This block does.
##
## The readout is the scVI-DENOISED AGTR1, not the raw lens. Raw AGTR1 and the BM
## panel both scale with capture depth, so a raw-lens association here is
## confounded by precisely the dropout that 03.agtr1_lenses.R showed manufactures
## AGTR1's apparent program bias -- reporting a raw-only result would repeat the
## error the three-lens analysis exists to catch. Raw and detection lenses are
## carried alongside as the dropout-sensitive comparison: raw-negative with
## denoised-null is the signature of shared dropout, not of biology (same logic as
## the ACTA2 control in pericyte_states/_h/05.acta2_control.R).
##
## Two levels are tested because they answer different questions, and panel F only
## answers the first: donor x cluster pseudobulk (do lungs/clusters with more AGTR1
## carry less BM?) and within-donor cell-level (does AGTR1 oppose BM INSIDE a
## lung?). The BM-vs-pseudotime trend is already known to be between-donor only --
## flat within donors, bm_continuum_summary.tsv median rho 0.02, p = 0.84 -- so the
## within-donor test is the one that would license an antagonism claim in the text.
agtr1_readme <- character(0)
agtr1_cols <- intersect(c("AGTR1_expr", "AGTR1_detect"), names(d))
dl <- copy(d)
if (!is.na(opt$denoise) && file.exists(opt$denoise)) {
    den <- fread(opt$denoise, select = c("index", "Model", "AGTR1_scvi"))
    den <- den[Model == opt$den_model]
    if (anyDuplicated(den$index))
        stop("denoised AGTR1 table has duplicate cells for model ", opt$den_model)
    dl <- merge(dl, den[, .(index, AGTR1_scvi)], by = "index", all.x = TRUE)
    message(sprintf("Denoised AGTR1 (%s): %d/%d cells (%.1f%%)", opt$den_model,
                    sum(is.finite(dl$AGTR1_scvi)), nrow(dl),
                    100 * sum(is.finite(dl$AGTR1_scvi)) / nrow(dl)))
    agtr1_cols <- c("AGTR1_scvi", agtr1_cols)
} else {
    warning("no --denoise table supplied; the denoised lens is the DISCRIMINATING ",
            "readout for AGTR1-vs-BM, so the raw-lens results below cannot settle ",
            "the question on their own.", call. = FALSE)
}

if (length(agtr1_cols)) {
    dl[, (agtr1_cols) := lapply(.SD, as.numeric), .SDcols = agtr1_cols]

    ## (i) donor x cluster pseudobulk, on the SAME units as every other test here
    apb <- dl[, lapply(.SD, mean, na.rm = TRUE),
              by = .(donor_id, pericyte_state), .SDcols = agtr1_cols]
    pba <- merge(pb, apb, by = c("donor_id", "pericyte_state"))
    for (cl in agtr1_cols)
        pba[[paste0(cl, "_z")]] <- z_within_dataset(pba[[cl]], pba$dataset)

    ## Unadjusted correlations, formatted to sit beside bm_vs_fibrillar_corr.tsv.
    corr_agtr1 <- rbindlist(lapply(agtr1_cols, function(lens) {
        ok <- is.finite(pba[[lens]]) & is.finite(pba$basement_membrane_score)
        if (sum(ok) < 10) return(NULL)
        x <- pba$basement_membrane_score[ok]; y <- pba[[lens]][ok]
        ct <- cor.test(x, y, method = "pearson")
        cs <- suppressWarnings(cor.test(x, y, method = "spearman"))
        data.table(comparison = paste0("basement_membrane vs ", lens),
                   pearson_r = unname(ct$estimate),
                   pearson_ci_lo = ct$conf.int[1], pearson_ci_hi = ct$conf.int[2],
                   pearson_p = ct$p.value, spearman_rho = unname(cs$estimate),
                   spearman_p = cs$p.value, n_units = sum(ok))
    }), fill = TRUE)
    write_tsv_safe(corr_agtr1, file.path(opt$outdir, "bm_vs_agtr1_corr.tsv"))

    ## Depth-adjusted mixed models. Outcome 1 is the BM score; outcome 2 is the
    ## fibrillar-orthogonalized residual, so a surviving slope there means AGTR1
    ## tracks BM specifically rather than the shared stromal-ECM axis.
    agtr1_fit <- function(lens, outcome, extra = character()) {
        lz <- paste0(lens, "_z")
        if (!lz %in% names(pba) || !outcome %in% names(pba)) return(NULL)
        f <- tryCatch(suppressMessages(lmer(
            reformulate(c(lz, extra, "(1 | study)", "(1 | donor_id)"), outcome),
            data = pba)), error = function(e) NULL)
        if (is.null(f)) return(NULL)
        cf <- summary(f)$coefficients
        if (!lz %in% rownames(cf)) return(NULL)
        ## Satterthwaite `df`/`Pr(>|t|)` come from lmerTest::lmer. Degrade to NA
        ## rather than erroring if a bare lme4 fit ever reaches here.
        col <- function(nm) if (nm %in% colnames(cf)) cf[lz, nm] else NA_real_
        data.table(outcome = outcome, lens = lens,
                   estimate = col("Estimate"), SE = col("Std. Error"),
                   df = col("df"), t_ratio = col("t value"),
                   p_value = col("Pr(>|t|)"),
                   depth_adjusted = "mean_log10_counts" %in% extra,
                   n_units = nrow(pba), n_donors = uniqueN(pba$donor_id))
    }
    mods_agtr1 <- rbindlist(c(
        lapply(agtr1_cols, agtr1_fit, outcome = "basement_membrane_score_z",
               extra = "mean_log10_counts"),
        lapply(agtr1_cols, agtr1_fit, outcome = "bm_resid")), fill = TRUE)
    if (nrow(mods_agtr1)) {
        mods_agtr1[, p_BH := p.adjust(p_value, method = "BH"), by = outcome]
        write_tsv_safe(mods_agtr1, file.path(opt$outdir, "bm_vs_agtr1_models.tsv"))
    }

    ## (ii) within-donor, cell level: per-donor Spearman then a one-sample test on
    ## the donor rhos. Never pools cells across donors -- the same pattern as the
    ## continuum block below, so the two are directly comparable.
    wd_agtr1 <- rbindlist(lapply(agtr1_cols, function(lens) {
        sub <- dl[is.finite(get(lens)) & is.finite(basement_membrane_score)]
        dr <- sub[, .(n = .N, rho = suppressWarnings(
                          cor(get(lens), basement_membrane_score, method = "spearman"))),
                  by = donor_id][n >= 20 & is.finite(rho)]
        if (nrow(dr) < 3) return(NULL)
        w <- suppressWarnings(wilcox.test(dr$rho))
        data.table(lens = lens, n_donors = nrow(dr),
                   median_rho = median(dr$rho),
                   q25 = unname(quantile(dr$rho, 0.25)),
                   q75 = unname(quantile(dr$rho, 0.75)),
                   p_wilcox = w$p.value)
    }), fill = TRUE)
    if (nrow(wd_agtr1)) {
        wd_agtr1[, p_BH := p.adjust(p_wilcox, method = "BH")]
        write_tsv_safe(wd_agtr1, file.path(opt$outdir, "bm_vs_agtr1_within_donor.tsv"))
    }

    message("AGTR1 vs BM -- pseudobulk correlations:");     print(corr_agtr1)
    message("AGTR1 vs BM -- depth-adjusted mixed models:"); print(mods_agtr1)
    message("AGTR1 vs BM -- within-donor Spearman:");       print(wd_agtr1)
    agtr1_readme <- c(
        "",
        "Direct AGTR1 vs BM (denoised lens is the readout; raw/detection are the",
        "dropout-sensitive comparison, NOT the answer):",
        paste(utils::capture.output(print(mods_agtr1)), collapse = "\n"),
        "",
        "Within-donor (cell-level Spearman per donor, >=20 cells, one-sample test):",
        paste(utils::capture.output(print(wd_agtr1)), collapse = "\n"))
} else {
    message("No AGTR1 columns available -- skipping the AGTR1-vs-BM test.")
}

## ------------------------------------------- BM along the injury continuum ----
if (!is.na(opt$continuum) && file.exists(opt$continuum)) {
    cont <- fread(opt$continuum)
    setnames(cont, 1, "index")
    pt_col <- intersect(c("dpt_pseudotime", "pseudotime"), names(cont))
    if (length(pt_col)) {
        ## Score columns carried onto the continuum. The first three metrics
        ## below reproduce the published rho_bm / rho_fib / rho_switch exactly;
        ## the 2026-08-10 additions resolve the fibrillar side into the
        ## fibril-forming (I/III) and regulatory (V/XI) blocks and add the
        ## ambient-tracer trend as a control -- if soup burden itself rises along
        ## pseudotime, the fibrillar rise is confounded and cannot be claimed.
        want <- c("basement_membrane_score", "fibrillar_ecm_score",
                  "fibrillar_core_score", "fibrillar_minor_score",
                  "ambient_tracer_score")
        want <- intersect(want, names(d))
        dc <- merge(d[, .SD, .SDcols = c("index", "donor_id",
                                         "log10_total_counts", want)],
                    cont[, .SD, .SDcols = c("index", pt_col[1])], by = "index")
        setnames(dc, pt_col[1], "pt")

        ## metric -> per-cell vector, so every rho is computed the same way.
        metrics <- list(rho_bm = "basement_membrane_score",
                        rho_fib = "fibrillar_ecm_score",
                        rho_core = "fibrillar_core_score",
                        rho_minor = "fibrillar_minor_score",
                        rho_tracer = "ambient_tracer_score")
        metrics <- metrics[vapply(metrics, `%in%`, logical(1), want)]
        switches <- list(
            rho_switch = c("basement_membrane_score", "fibrillar_ecm_score"),
            rho_switch_core = c("basement_membrane_score", "fibrillar_core_score"))
        switches <- switches[vapply(switches, function(p) all(p %in% want),
                                    logical(1))]

        ## Per-donor Spearman, then a one-sample test on the donor rhos -- never
        ## pools cells across donors (the existing pseudotime-trend pattern).
        spear <- function(x, pt) suppressWarnings(cor(x, pt, method = "spearman"))
        donor_rho <- dc[, c(
            .(n = .N),
            lapply(metrics, function(cl) spear(get(cl), pt)),
            lapply(switches, function(p) spear(get(p[1]) - get(p[2]), pt))),
            by = donor_id][n >= 20]
        write_tsv_safe(donor_rho, file.path(opt$outdir, "bm_continuum_donor_rho.tsv"))

        summ <- rbindlist(lapply(c(names(metrics), names(switches)), function(cc) {
            v <- donor_rho[[cc]]; v <- v[is.finite(v)]
            if (length(v) < 3) return(NULL)
            w <- suppressWarnings(wilcox.test(v))
            data.table(metric = cc, n_donors = length(v), median_rho = median(v),
                       q25 = unname(quantile(v, 0.25)),
                       q75 = unname(quantile(v, 0.75)),
                       p_wilcox = w$p.value)
        }), fill = TRUE)
        ## BH over the three original metrics only, so the published adjusted
        ## p-values do not move because the panel grew; the expansion metrics
        ## are corrected as their own family.
        orig <- c("rho_bm", "rho_fib", "rho_switch")
        summ[, bh_family := fifelse(metric %in% orig, "original", "expansion")]
        summ[, p_BH := p.adjust(p_wilcox, method = "BH"), by = bh_family]
        write_tsv_safe(summ, file.path(opt$outdir, "bm_continuum_summary.tsv"))
        message("Continuum trends (donor-level Spearman vs pseudotime):")
        print(summ)
    }
}

readme <- c(
    "BM pericyte axis -- generated summary",
    sprintf("Donor x cluster units (>=%d cells): %d; donors: %d",
            opt$min_cells, nrow(pb), uniqueN(pb$donor_id)),
    sprintf("Variance components: study SD %.3f, residual SD %.3f -> %s",
            ifelse(length(sd_study), sd_study, NA_real_), sd_resid,
            ifelse(study_dominated, "STUDY DOMINATED - claim not supportable", "ok")),
    "",
    "Grouping variable is pericyte_state (panel-independent Leiden clusters).",
    "state_program is NOT used as an outcome grouping: after the gate escalated it",
    "is derived from the BM score, so testing against it would be circular.",
    "",
    "BM vs fibrillar (donor x cluster pseudobulk):",
    paste(utils::capture.output(print(corr_tbl)), collapse = "\n"),
    agtr1_readme)
writeLines(readme, file.path(opt$outdir, "bm_pericyte_axis_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
