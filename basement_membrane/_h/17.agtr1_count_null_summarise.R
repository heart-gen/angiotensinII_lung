#!/usr/bin/env Rscript
## Concatenate the array chunks from 16.agtr1_count_null.R and place the
## arbiter's observed estimate against the pooled detection-matched null.
suppressPackageStartupMessages({ library(data.table); library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
    make_option("--dir", type = "character", help = "stats_data holding the chunks"),
    make_option("--observed", type = "character"),
    make_option("--outdir", type = "character"))))

f <- list.files(opt$dir, pattern = "^agtr1_count_null_fits_chunk[0-9]+\\.tsv$",
                full.names = TRUE)
if (!length(f)) stop("no chunk files in ", opt$dir)
res <- rbindlist(lapply(f, fread), fill = TRUE)
message(sprintf("Pooled %d chunks -> %d fits over %d genes (%d non-converged)",
                length(f), nrow(res), uniqueN(res$gene), sum(!res$converged, na.rm = TRUE)))
fwrite(res, file.path(opt$outdir, "agtr1_count_null_fits.tsv"), sep = "\t")

obs <- fread(opt$observed)
obs <- obs[level == "cell" & spec == "primary" & model == "NB GLMM" &
           predictor %in% unique(res$predictor)]
summ <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
    row <- obs[i]
    nd <- res[predictor == row$predictor & is.finite(estimate) & converged == TRUE]
    if (!nrow(nd)) return(NULL)
    mu <- mean(nd$estimate); sdev <- sd(nd$estimate)
    emp <- (1 + sum(abs(nd$estimate - mu) >= abs(row$estimate - mu))) / (1 + nrow(nd))
    data.table(predictor = row$predictor, level = row$level, spec = row$spec,
               beta_observed = row$estimate, p_model = row$p_value,
               n_null = nrow(nd), null_mean = mu, null_sd = sdev,
               null_q025 = quantile(nd$estimate, 0.025),
               null_q975 = quantile(nd$estimate, 0.975),
               z_vs_null = (row$estimate - mu) / sdev, empirical_p = emp,
               detectable_effect_80 = quantile(abs(nd$estimate - mu), 0.80),
               verdict = fifelse(emp < 0.05, "outside its matched null",
                                 "INSIDE its matched null"))
}), fill = TRUE)
fwrite(summ, file.path(opt$outdir, "agtr1_count_null_summary.tsv"), sep = "\t")
message("\n---- count model vs its detection-matched null ----")
print(summ[, .(predictor, beta_observed, null_mean, null_sd, z_vs_null,
               empirical_p, verdict)])
