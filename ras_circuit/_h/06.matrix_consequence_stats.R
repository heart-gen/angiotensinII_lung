## Figure 5E -- is the AT1R-response programme coupled to basement-membrane versus
## fibrillar matrix?
##
## THE ESTIMAND IS THE CONTRAST (BM - fibrillar), never either matrix score alone
## (memory: basement-membrane-program; the contrast's null is centred near zero,
## the single-score null near +0.45). BM-alone and fibrillar-alone rows are
## written with claimable = FALSE for completeness.
##
## Unit: donor x pericyte_state, the table 05 wrote. Primary:
##   lmer(bm_minus_fib_z ~ at1r_z + mean_log10_counts + (1|donor_id) + (1|study))
## where bm_minus_fib = basement_membrane_score - fibrillar_collagen_score (the
## contrast the AGTR1 count-model claim uses). Arms: + ambient tracer (it rises
## along pseudotime, rho +0.203, so it is a live confound), depth spline, the
## fibrillar_ecm contrast variant. Null: the primary refitted on every
## detection-matched null panel; p_emp against the null's own centre.
## Within-donor: per-donor Spearman of (BM - fib) with the score, partialled on
## depth and tracer; Wilcoxon on donor rhos.
## Context rows (READ, not recomputed): bm_continuum_summary.tsv (rho_switch,
## rho_bm, rho_tracer) and the AGTR1 count-model contrast, for the legend caveats.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(lmerTest)
    library(splines); library(parallel)
})
source("../_h/_stats_common.R")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--units", default = "./at1r_unit_pseudobulk.tsv.gz"),
    make_option("--meta", default = "./at1r_response_metadata.tsv.gz"),
    make_option("--null-pb", dest = "null_pb", default = "./at1r_null_pseudobulk.tsv.gz"),
    make_option("--bm-stats", dest = "bm_stats",
                default = "../../basement_membrane/_m/stats_data"),
    make_option("--outdir", default = "./stats_data"),
    make_option("--max-null", type = "integer", default = 0L, dest = "max_null")
)))
CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

pb <- read_req(opt$units)
pb[, `:=`(bm_minus_fib = basement_membrane_score - fibrillar_collagen_score,
          bm_minus_fibecm = basement_membrane_score - fibrillar_ecm_score)]
for (v in c("bm_minus_fib", "bm_minus_fibecm", "basement_membrane_score",
            "fibrillar_collagen_score", "at1r_response_score", "ambient_tracer_score"))
    pb[[paste0(v, "_z")]] <- z_within_dataset(pb[[v]], pb$dataset)

fit1 <- function(y, rhs, data, spec, claimable = TRUE, pred = "at1r_response_score_z") {
    f <- as.formula(paste(y, "~", pred, "+", rhs, "+ (1|donor_id) + (1|study)"))
    fit <- try(suppressMessages(lmer(f, data = data)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    r <- tidy_row(fit, pred, "pseudobulk", "lmer(+1|donor)+(1|study)", nrow(data),
                  uniqueN(data$donor_id), spec = spec)
    r[, `:=`(outcome = y, covariates = rhs, claimable = claimable)]
}
res <- rbindlist(list(
    fit1("bm_minus_fib_z", "mean_log10_counts", pb, "primary"),
    fit1("bm_minus_fib_z", "mean_log10_counts + ambient_tracer_score_z", pb, "_tracer_adj"),
    fit1("bm_minus_fib_z", "ns(mean_log10_counts, 3)", pb, "_depth_spline"),
    fit1("bm_minus_fibecm_z", "mean_log10_counts", pb, "_fibrillar_ecm_contrast"),
    fit1("basement_membrane_score_z", "mean_log10_counts", pb, "_bm_alone", FALSE),
    fit1("fibrillar_collagen_score_z", "mean_log10_counts", pb, "_fib_alone", FALSE)),
    fill = TRUE)

## null on the primary
nl <- read_req(opt$null_pb)
nl[, `:=`(donor_id = as.character(donor_id), pericyte_state = as.character(pericyte_state))]
null_cols <- grep("^null_", names(nl), value = TRUE)
if (opt$max_null > 0) null_cols <- head(null_cols, opt$max_null)
pbn <- merge(pb[, .(donor_id = as.character(donor_id),
                    pericyte_state = as.character(pericyte_state),
                    study, dataset, mean_log10_counts, bm_minus_fib_z,
                    basement_membrane_score_z)],
             nl[, c("donor_id", "pericyte_state", null_cols), with = FALSE],
             by = c("donor_id", "pericyte_state"))
if (nrow(pbn) != nrow(pb)) stop("null pseudobulk does not cover the units")
for (cl in null_cols) pbn[[paste0(cl, "_z")]] <- z_within_dataset(pbn[[cl]], pbn$dataset)
nulls <- rbindlist(mclapply(null_cols, function(cl) {
    out <- lapply(c("bm_minus_fib_z", "basement_membrane_score_z"), function(y) {
        fit <- try(suppressMessages(lmer(as.formula(sprintf(
            "%s ~ %s_z + mean_log10_counts + (1|donor_id) + (1|study)", y, cl)),
            data = pbn)), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        data.table(panel = cl, outcome = y, estimate = fixef(fit)[[paste0(cl, "_z")]])
    })
    rbindlist(out)
}, mc.cores = CORES))
fwrite(nulls, file.path(opt$outdir, "matrix_vs_at1r_null.tsv"), sep = "\t")
res[, `:=`(null_mean = NA_real_, null_sd = NA_real_, p_emp = NA_real_)]
for (y in c("bm_minus_fib_z", "basement_membrane_score_z")) {
    nv <- nulls[outcome == y, estimate]
    i <- res$outcome == y & res$covariates == "mean_log10_counts"
    res[i, `:=`(null_mean = mean(nv), null_sd = sd(nv), p_emp = emp_p(estimate, nv))]
}
write_tsv_safe(res, file.path(opt$outdir, "matrix_vs_at1r.tsv"))
print(res[, .(spec, outcome, estimate, SE, p_value, null_mean, p_emp, claimable)])

## within-donor Spearman
m <- read_req(opt$meta); setnames(m, 1, "index")
m[, bm_minus_fib := basement_membrane_score - fibrillar_collagen_score]
spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
wd <- m[, if (.N >= 20) {
    rx <- resid(lm(at1r_response_score ~ log10_total_counts + ambient_tracer_score))
    ry <- resid(lm(bm_minus_fib ~ log10_total_counts + ambient_tracer_score))
    .(n = .N, rho = spear(at1r_response_score, bm_minus_fib), prho = spear(rx, ry))
}, by = donor_id]
fwrite(wd, file.path(opt$outdir, "matrix_vs_at1r_donor_rho.tsv"), sep = "\t")
ws <- rbindlist(lapply(c("rho", "prho"), function(k) {
    v <- wd[[k]][is.finite(wd[[k]])]
    data.table(kind = fifelse(k == "rho", "raw", "partial_depth_tracer"),
               n_donors = length(v), median_rho = median(v),
               q25 = quantile(v, 0.25), q75 = quantile(v, 0.75),
               p_wilcox = suppressWarnings(wilcox.test(v))$p.value)
}))
write_tsv_safe(ws, file.path(opt$outdir, "matrix_vs_at1r_donor_rho_summary.tsv"))

## context rows for the legend caveats (read, not recomputed)
ctx <- list()
f1 <- file.path(opt$bm_stats, "bm_continuum_summary.tsv")
if (file.exists(f1)) ctx$cont <- fread(f1)[metric %in% c("rho_bm", "rho_switch", "rho_tracer",
                                                         "rho_fib", "rho_collagen")][
    , source := "basement_membrane/_m/stats_data/bm_continuum_summary.tsv"]
f2 <- file.path(opt$bm_stats, "agtr1_count_models.tsv")
if (file.exists(f2)) ctx$cnt <- fread(f2)[grepl("bm_minus_fibrillar", predictor)][
    , source := "basement_membrane/_m/stats_data/agtr1_count_models.tsv"]
write_tsv_safe(rbindlist(ctx, fill = TRUE), file.path(opt$outdir, "matrix_context_rows.tsv"))

cat("\nReproducibility information:\n"); print(sessionInfo())
