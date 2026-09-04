#!/usr/bin/env Rscript
##
## Is AGTR1 -> matrix being tested against the right reference?
##
## The TGF-beta specificity run found that the null for a panel score against
## `basement_membrane_score_z` is centred at +0.45, not 0. 04.bm_state_stats.R
## runs the identical model with AGTR1 as predictor and reports AGTR1 -> BM as
## significantly POSITIVE (+0.170 expr, BH 0.025; +0.209 detect, BH 0.017). If
## the single-gene null is also above zero, those estimates are measured against
## the wrong reference and the claim's sign is in question, not just its size.
##
## AGTR1 is one gene, not a panel, so the panel result does not transfer. This
## fits the SAME model to 1,000 detection-matched random single genes, under
## both lenses the module uses (mean log-expression and detection rate).
##
## Design is copied from assoc_block() in 04.bm_state_stats.R -- same units, same
## depth covariate, same random effects, same within-dataset z -- because the
## whole point is comparability with the published estimates.
##
## AGTR1_scvi has no null here: building one would need scVI retrained per null
## gene. Reported as a limitation, not approximated.
suppressPackageStartupMessages({
    library(data.table); library(optparse); library(lme4); library(lmerTest)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--bm-meta", type = "character", dest = "bm_meta"),
    make_option("--state-meta", type = "character", dest = "state_meta"),
    make_option("--null-pseudobulk", type = "character", dest = "null_pb"),
    make_option("--observed", type = "character", dest = "observed",
                help = "bm_vs_agtr1_models.tsv, the published estimates"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

bm <- fread(opt$bm_meta); st <- fread(opt$state_meta)
setnames(st, 1, "index"); setnames(bm, 1, "index")
drop <- intersect(setdiff(names(bm), "index"), names(st))
d <- merge(st, bm[, .SD, .SDcols = setdiff(names(bm), drop)], by = "index")
d[, pericyte_state := factor(pericyte_state)]
if ("dataset" %in% names(d)) d[, dataset := as.character(dataset)]

score_cols <- intersect(c("basement_membrane_score", "fibrillar_collagen_score"),
                        names(d))
pb <- d[, c(lapply(.SD, mean, na.rm = TRUE),
            .(n_cells = .N, study = first(study), dataset = first(dataset),
              mean_log10_counts = mean(log10_total_counts, na.rm = TRUE))),
        by = .(donor_id, pericyte_state), .SDcols = score_cols]
pb <- pb[n_cells >= opt$min_cells]
for (cl in score_cols)
    pb[[paste0(cl, "_z")]] <- z_within_dataset(pb[[cl]], pb$dataset)
pb[, bm_minus_fibrillar := basement_membrane_score_z - fibrillar_collagen_score_z]
pb[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]
message(sprintf("Units: %d across %d donors", nrow(pb), uniqueN(pb$donor_id)))

OUTCOMES <- c("basement_membrane_score_z", "fibrillar_collagen_score_z",
              "bm_minus_fibrillar_z")

nl <- fread(opt$null_pb)
nl[, `:=`(donor_id = as.character(donor_id),
          pericyte_state = as.character(pericyte_state))]
pbk <- pb[, .(donor_id = as.character(donor_id),
              pericyte_state = as.character(pericyte_state),
              study, dataset, mean_log10_counts,
              basement_membrane_score_z, fibrillar_collagen_score_z,
              bm_minus_fibrillar_z)]
pbn <- merge(pbk, nl, by = c("donor_id", "pericyte_state"))
if (nrow(pbn) != nrow(pb))
    stop("null covers ", nrow(pbn), " of ", nrow(pb), " units")

null_cols <- grep("^null(expr|det)_", names(pbn), value = TRUE)
message(sprintf("Null predictors: %d, on %d units", length(null_cols), nrow(pbn)))
for (cl in null_cols)
    pbn[[paste0(cl, "_z")]] <- z_within_dataset(pbn[[cl]], pbn$dataset)

fit_one <- function(dt, outcome, pred) {
    f <- reformulate(c(pred, "mean_log10_counts", "(1 | study)", "(1 | donor_id)"),
                     response = outcome)
    fit <- try(suppressMessages(lmerTest::lmer(
        f, data = dt, control = lmerControl(calc.derivs = FALSE))), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    co <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    setnames(co, c("term", "estimate", "SE", "df", "t_ratio", "p_value"))
    co <- co[term == pred]
    g <- tryCatch(max(abs(fit@optinfo$derivs$gradient)), error = function(e) NA_real_)
    co[, `:=`(outcome = outcome, converged = is.na(g) || g < 0.01)]
    co[]
}

res <- rbindlist(lapply(seq_along(null_cols), function(i) {
    cl <- null_cols[i]
    if (i %% 250 == 0) message(sprintf("  fit %d/%d", i, length(null_cols)))
    rbindlist(lapply(OUTCOMES, function(o) {
        r <- fit_one(pbn, o, paste0(cl, "_z"))
        if (is.null(r)) return(NULL)
        r[, `:=`(panel = cl,
                 lens = fifelse(grepl("^nullexpr_", cl), "expr", "detect"))][]
    }), fill = TRUE)
}), fill = TRUE)
fwrite(res, file.path(opt$outdir, "agtr1_null_models.tsv"), sep = "\t")

## ---- compare the published estimates against their own matched null ---------
obs <- fread(opt$observed)
obs <- obs[outcome %in% OUTCOMES &
           predictor %in% c("AGTR1_expr", "AGTR1_detect", "AGTR1_scvi")]
obs[, lens := fcase(predictor == "AGTR1_expr", "expr",
                    predictor == "AGTR1_detect", "detect",
                    default = NA_character_)]

summ <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
    row <- obs[i]
    if (is.na(row$lens)) {
        ## The denoised lens has no matched null -- say so in the table rather
        ## than leaving a gap a reader fills in with an assumption.
        return(data.table(predictor = row$predictor, outcome = row$outcome,
                          beta_observed = row$estimate, p_model = row$p_value,
                          n_null = NA_integer_, null_mean = NA_real_,
                          null_sd = NA_real_, z_vs_null = NA_real_,
                          empirical_p = NA_real_, detectable_effect_80 = NA_real_,
                          verdict = "no matched null (scVI cannot be permuted here)"))
    }
    nd <- res[lens == row$lens & outcome == row$outcome & converged == TRUE]
    mu <- mean(nd$estimate); sdev <- sd(nd$estimate)
    ## p and power referred to the null's OWN centre -- the null is not at zero.
    emp <- (1 + sum(abs(nd$estimate - mu) >= abs(row$estimate - mu))) / (1 + nrow(nd))
    data.table(predictor = row$predictor, outcome = row$outcome,
               beta_observed = row$estimate, p_model = row$p_value,
               n_null = nrow(nd), null_mean = mu, null_sd = sdev,
               null_q025 = quantile(nd$estimate, 0.025),
               null_q975 = quantile(nd$estimate, 0.975),
               z_vs_null = (row$estimate - mu) / sdev,
               empirical_p = emp,
               detectable_effect_80 = quantile(abs(nd$estimate - mu), 0.80),
               verdict = fifelse(emp < 0.05, "outside its matched null",
                                 "INSIDE its matched null"))
}), fill = TRUE)
fwrite(summ, file.path(opt$outdir, "agtr1_null_summary.tsv"), sep = "\t")

message("\n---- AGTR1 vs its detection-matched null ----")
print(summ[, .(predictor, outcome, beta_observed, null_mean, z_vs_null,
               empirical_p, verdict)])
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir, "agtr1_null_sessionInfo.txt"))
