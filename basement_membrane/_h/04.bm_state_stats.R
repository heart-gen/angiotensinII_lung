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

## ------------------------------------------- BM along the injury continuum ----
if (!is.na(opt$continuum) && file.exists(opt$continuum)) {
    cont <- fread(opt$continuum)
    setnames(cont, 1, "index")
    pt_col <- intersect(c("dpt_pseudotime", "pseudotime"), names(cont))
    if (length(pt_col)) {
        dc <- merge(d[, .(index, donor_id, basement_membrane_score,
                          fibrillar_ecm_score, log10_total_counts)],
                    cont[, .SD, .SDcols = c("index", pt_col[1])], by = "index")
        setnames(dc, pt_col[1], "pt")
        ## Per-donor Spearman, then a one-sample test on the donor rhos -- never
        ## pools cells across donors (the existing pseudotime-trend pattern).
        donor_rho <- dc[, .(n = .N,
                            rho_bm = suppressWarnings(
                                cor(basement_membrane_score, pt, method = "spearman")),
                            rho_fib = suppressWarnings(
                                cor(fibrillar_ecm_score, pt, method = "spearman")),
                            rho_switch = suppressWarnings(
                                cor(basement_membrane_score - fibrillar_ecm_score,
                                    pt, method = "spearman"))),
                        by = donor_id][n >= 20]
        write_tsv_safe(donor_rho, file.path(opt$outdir, "bm_continuum_donor_rho.tsv"))
        summ <- rbindlist(lapply(c("rho_bm", "rho_fib", "rho_switch"), function(cc) {
            v <- donor_rho[[cc]]; v <- v[is.finite(v)]
            if (length(v) < 3) return(NULL)
            w <- suppressWarnings(wilcox.test(v))
            data.table(metric = cc, n_donors = length(v), median_rho = median(v),
                       p_wilcox = w$p.value)
        }), fill = TRUE)
        summ[, p_BH := p.adjust(p_wilcox, method = "BH")]
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
    paste(utils::capture.output(print(corr_tbl)), collapse = "\n"))
writeLines(readme, file.path(opt$outdir, "bm_pericyte_axis_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
