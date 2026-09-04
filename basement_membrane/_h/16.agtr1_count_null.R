#!/usr/bin/env Rscript
##
## Does the COUNT MODEL -- this module's designated arbiter -- survive a
## detection-matched null?
##
## 14.agtr1_null_models.R found every AGTR1-versus-matrix REGRESSION estimate
## inside its detection-matched null, including the headline
## AGTR1 -> BM - fibrillar (+0.197, BH 0.006; emp. p 0.28). But those
## regressions are not the arbiter. The standing rule in this repo is that the
## count model -- AGTR1 integer counts as the RESPONSE of an NB GLMM with a
## library-size offset, no imputation -- arbitrates group contrasts. Its design
## is inverted relative to the regressions and it handles depth twice on
## purpose, so it may behave differently. Whether the biological claim stands or
## falls depends on this file.
##
## Method: substitute each of the N detection-matched null genes for AGTR1 as
## the response, refit the SAME specifications, and ask where the real AGTR1
## estimate sits in the resulting distribution.
##
## Two specs, both taken verbatim from 10.agtr1_count_models.R:
##   pseudobulk NB   : y_sum ~ score_z + mean_log10_counts + (1|study) + (1|donor) + offset(log(total_sum))
##   cell-level NB   : y ~ score_z + log10_total_counts + (1|study) + (1|donor) + offset(log(raw_total_counts))
## The cell-level fit is the one memory cites as the arbiter (+0.080, p 3.1e-7).
suppressPackageStartupMessages({
    library(data.table); library(optparse); library(lme4)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input", type = "character",
                help = "agtr1_count_input.tsv.gz (the arbiter's own input)"),
    make_option("--null-counts", type = "character", dest = "null_counts",
                help = "null_gene_count_input.tsv.gz"),
    make_option("--observed", type = "character",
                help = "agtr1_count_models.tsv"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--level", type = "character", default = "cell",
                help = "'cell' (the arbiter) or 'pseudobulk'"),
    make_option("--max-genes", type = "integer", default = 0L, dest = "max_genes",
                help = "0 = all null genes"),
    ## glmer.nb on 11,680 cells is slow enough that 350 fits run for hours
    ## serially, so the null genes are split across a SLURM array and the
    ## per-chunk fits are concatenated by 17.agtr1_count_null_summarise.R.
    make_option("--chunk", type = "integer", default = 0L,
                help = "1-based chunk index; 0 = run everything in one process"),
    make_option("--n-chunks", type = "integer", default = 1L, dest = "n_chunks")
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

d  <- fread(opt$input)
nc <- fread(opt$null_counts)
setnames(d, 1, "index"); setnames(nc, 1, "index")
nc[, raw_total_counts := NULL]          # already in `d`; keep one copy
dl <- merge(d, nc, by = "index")
message(sprintf("Merged %d cells x %d null genes", nrow(dl),
                sum(grepl("^null_", names(dl)))))
if (nrow(dl) != nrow(d))
    stop("null counts do not cover the arbiter's cells: ", nrow(dl), " of ", nrow(d))

for (s in c("basement_membrane_score", "fibrillar_collagen_score"))
    dl[[paste0(s, "_z")]] <- z_within_dataset(dl[[s]], dl$dataset)
dl[, bm_minus_fibrillar := basement_membrane_score - fibrillar_collagen_score]
dl[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]

PREDICTORS <- c("basement_membrane_score_z", "bm_minus_fibrillar_z")
null_genes <- grep("^null_", names(dl), value = TRUE)
if (opt$max_genes > 0) null_genes <- head(null_genes, opt$max_genes)
if (opt$chunk > 0) {
    grp <- rep(seq_len(opt$n_chunks), length.out = length(null_genes))
    null_genes <- null_genes[grp == opt$chunk]
    message(sprintf("Chunk %d of %d", opt$chunk, opt$n_chunks))
}
message(sprintf("Null genes to fit: %d  x %d predictors = %d fits",
                length(null_genes), length(PREDICTORS),
                length(null_genes) * length(PREDICTORS)))

## Same NB-with-offset specification the arbiter uses. glmer.nb is slow and
## occasionally fails to converge; a failure is recorded and dropped from the
## null rather than silently contributing a meaningless coefficient.
fit_null <- function(y, pred) {
    f <- as.formula(sprintf(
        "%s ~ %s + log10_total_counts + (1|study) + (1|donor_id) + offset(log(raw_total_counts))",
        y, pred))
    fit <- tryCatch(suppressMessages(glmer.nb(f, data = dl)),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients
    if (!pred %in% rownames(co)) return(NULL)
    data.table(gene = sub("^null_", "", y), predictor = pred,
               estimate = co[pred, 1], SE = co[pred, 2],
               z_value = co[pred, 3], p_value = co[pred, 4],
               singular = isSingular(fit))
}

t0 <- Sys.time()
res <- rbindlist(lapply(seq_along(null_genes), function(i) {
    y <- null_genes[i]
    if (i %% 10 == 0)
        message(sprintf("  gene %d/%d (%.1f min elapsed)", i, length(null_genes),
                        as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    rbindlist(lapply(PREDICTORS, function(p) fit_null(y, p)), fill = TRUE)
}), fill = TRUE)
fits_out <- if (opt$chunk > 0)
    sprintf("agtr1_count_null_fits_chunk%02d.tsv", opt$chunk) else
    "agtr1_count_null_fits.tsv"
fwrite(res, file.path(opt$outdir, fits_out), sep = "\t")
message(sprintf("Converged null fits: %d of %d attempted",
                nrow(res), length(null_genes) * length(PREDICTORS)))

if (opt$chunk > 0) {
    message("Chunk complete; summary is built by 17.agtr1_count_null_summarise.R")
    quit(save = "no", status = 0)
}

obs <- fread(opt$observed)
obs <- obs[level == "cell" & spec == "primary" & model == "NB GLMM" &
           predictor %in% PREDICTORS]

summ <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
    row <- obs[i]
    nd <- res[predictor == row$predictor & is.finite(estimate)]
    if (!nrow(nd)) return(NULL)
    mu <- mean(nd$estimate); sdev <- sd(nd$estimate)
    ## Referred to the null's OWN centre: the null need not sit at zero.
    emp <- (1 + sum(abs(nd$estimate - mu) >= abs(row$estimate - mu))) / (1 + nrow(nd))
    data.table(predictor = row$predictor, level = row$level, spec = row$spec,
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
fwrite(summ, file.path(opt$outdir, "agtr1_count_null_summary.tsv"), sep = "\t")

message("\n---- count model vs its detection-matched null ----")
print(summ[, .(predictor, beta_observed, null_mean, null_sd, z_vs_null,
               empirical_p, verdict)])
writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir, "agtr1_count_null_sessionInfo.txt"))
