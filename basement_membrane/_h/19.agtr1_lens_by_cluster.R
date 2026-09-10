## The three observational/imputed AGTR1 lenses, refit at the SAME unit as the
## count-model arbiter, so panel D can put all four on one panel honestly.
##
## WHY THIS SCRIPT EXISTS. `pericyte_states/_h/03.agtr1_lenses.R` fits each lens
## as `lmer(y ~ state_program + (1 | donor_id))` on all 11,680 CELLS -- no study
## term, no depth covariate, and (for the denoised lens) no propagation of
## imputation uncertainty. `10.agtr1_count_models.R` fits the count arbiter on
## 214 donor x cluster pseudobulks with `(1 | study) + (1 | donor_id)` and a
## library-size offset. Drawing those two side by side compares error bars that
## were never on the same footing: on a common log(AGTR1 per 10^4 transcripts)
## scale the cell-level denoised SEs are 0.061-0.076 and the count SEs are
## 0.131-0.204, a 2.0-2.7x gap that is almost entirely the unit of analysis.
## Refit at the pseudobulk unit the same denoised lens gives SEs of 0.110-0.221,
## within 5-43% of the count model's. The gap is the unit, not the lens, and a
## figure that lets a reader conclude otherwise is misleading.
##
## This also replaces `agtr1_lens_by_cluster_emmeans.tsv`, which existed in
## `_m/stats_data` with NO producing script anywhere in the repo (it carried a
## superseded "denoised OLD (invalid)" series from the scVI model that failed its
## validity gate). An orphan table cannot feed a main figure.
##
## GROUPING IS `pericyte_state` (Leiden P0-P5), NOT `state_program`. Two reasons,
## both from 10.agtr1_count_models.R's README:
##   - non-circularity: the clusters come from X_pca_harmony over 2,000 HVGs that
##     exclude AGTR1 (highly_variable = FALSE), so the grouping is independent of
##     the readout. `state_program` is a marker-panel argmax.
##   - power: at program level only 1 of 5 count-model specs separates
##     basement-membrane from vascular-stabilizing; at cluster level all twelve
##     pseudobulk BM-vs-VS contrasts across the three pseudobulk specs agree in
##     sign, 10 of 12 significant.
##
## THE PSEUDOBULK MUST MATCH 10.agtr1_count_models.R EXACTLY -- same cell table,
## same `donor_id x pericyte_state` unit, same `n_cells >= --min-cells` filter --
## or the panel is again comparing different footings. It is therefore built from
## `agtr1_count_input.tsv.gz` (the count model's own input), with AGTR1_expr /
## AGTR1_detect / AGTR1_scvi merged on by barcode. The script refuses to run if
## the resulting unit count does not match what the count table reports.
##
## SCALES, and why the denoised lens is logged here. `AGTR1_scvi` is a RATE
## (denoised AGTR1 per 10^4 transcripts), the count model reports a LOG rate.
## The pseudobulk denoised mean is logged so the two share an axis and an SE
## interpretation; the SE travels by the delta method, SE(log x) = SE(x)/x.
## AGTR1_expr (mean log1p-normalised) and AGTR1_detect (detected fraction) are
## NOT rates and are left on their own scale -- the figure facets on `scale_grp`
## for exactly this reason. Never centre these two against the other two.
##
## Input:  agtr1_count_input.tsv.gz  (09.agtr1_counts_prep.py)
##         pericytes_states_metadata.tsv.gz  (AGTR1_expr, AGTR1_detect)
##         pericytes_airspace_denoising.tsv  (AGTR1_scvi, retrained model)
##         agtr1_count_by_cluster.tsv  (10.agtr1_count_models.R; read only to
##                                      cross-check the unit count)
## Output: agtr1_lens_by_cluster_emmeans.tsv, agtr1_lens_by_cluster_posthoc.tsv
suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(lmerTest)
    library(emmeans)
})

option_list <- list(
    make_option("--input", type = "character", default = "./agtr1_count_input.tsv.gz"),
    make_option("--state-meta", type = "character", dest = "state_meta",
                default = "../../pericyte_states/_m/pericytes_states_metadata.tsv.gz"),
    make_option("--denoise", type = "character",
                default = "../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv"),
    make_option("--den-model", type = "character", dest = "den_model",
                default = "Pericyte-only-trained"),
    make_option("--count-table", type = "character", dest = "count_table",
                default = "./stats_data/agtr1_count_by_cluster.tsv"),
    make_option("--outdir", type = "character", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--seed", type = "integer", default = 13L)
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## The by-cluster fits are on 214 pseudobulk rows, far under emmeans' default
## ceiling, so Satterthwaite d.f. are computed without raising any limit. That is
## a deliberate difference from 03.agtr1_lenses.R, which had to raise
## lmerTest.limit because it fits on 11,680 cells.

## ---- assemble the cell table -------------------------------------------
dl <- fread(opt$input)
stopifnot(all(c("index", "donor_id", "study", "pericyte_state",
                "log10_total_counts", "raw_total_counts") %in% names(dl)))

meta <- fread(opt$state_meta)
setnames(meta, 1, "barcode")
meta <- meta[, .(barcode, AGTR1_expr, AGTR1_detect)]

den <- fread(opt$denoise, select = c("index", "Model", "AGTR1_scvi"))
setnames(den, "index", "barcode")
den <- den[Model == opt$den_model]
stopifnot(!anyDuplicated(den$barcode))
message(sprintf("denoised model: %s (%d cells)", opt$den_model, nrow(den)))

setnames(dl, "index", "barcode")
n0 <- nrow(dl)
dl <- merge(dl, meta, by = "barcode")
dl <- merge(dl, den[, .(barcode, AGTR1_scvi)], by = "barcode")
message(sprintf("cells: %d of %d count-model cells carry all three lenses (%.1f%%)",
                nrow(dl), n0, 100 * nrow(dl) / n0))
## A silent shortfall here would change which cells each lens is fit on, which is
## the whole failure this script exists to remove.
if (nrow(dl) != n0)
    stop(sprintf("lens merge dropped %d of %d cells; the pseudobulk would no ",
                 n0 - nrow(dl), n0),
         "longer match the count model's unit")

## ---- pseudobulk: donor x Leiden cluster, identical to 10.agtr1_count_models.R
pb <- dl[, .(n_cells = .N,
             mean_log10_counts = mean(log10_total_counts),
             AGTR1_expr = mean(AGTR1_expr),
             AGTR1_detect = mean(AGTR1_detect),
             AGTR1_scvi = mean(AGTR1_scvi),
             study = study[1], dataset = dataset[1]),
         by = .(donor_id, pericyte_state)]
pb <- pb[n_cells >= opt$min_cells]
pb[, cl := factor(pericyte_state)]
pb[, donor_id := factor(donor_id)]
message(sprintf("pseudobulk units (>=%d cells): %d from %d donors",
                opt$min_cells, nrow(pb), uniqueN(pb$donor_id)))

## Cross-check against the count model rather than trusting the recipe. If the
## two ever diverge -- a changed --min-cells, a re-run upstream prep -- the panel
## silently stops being unit-matched, so this is a hard stop.
if (file.exists(opt$count_table)) {
    cnt <- fread(opt$count_table)
    n_expect <- sum(unique(cnt[level == "pseudobulk",
                               .(pericyte_state, n_units)])$n_units)
    if (n_expect != nrow(pb))
        stop(sprintf("pseudobulk unit mismatch: %d built here vs %d in %s",
                     nrow(pb), n_expect, opt$count_table))
    message(sprintf("unit count matches %s (%d)", basename(opt$count_table), n_expect))
} else {
    warning(opt$count_table, " absent; unit count NOT cross-checked")
}

## ---- fit each lens ------------------------------------------------------
## `mean_log10_counts` is a free covariate in every lens for the same reason it
## is in the count model: sequencing depth varies ~2-fold across these clusters
## and drives both the raw and the detection readout. Leaving it out is what
## manufactures the raw lens's apparent cluster bias in the first place.
LENSES <- c(AGTR1_expr = "raw AGTR1_expr",
            AGTR1_detect = "AGTR1_detect",
            AGTR1_scvi = "denoised (retrained)")

emm_all <- list(); ph_all <- list()
for (nm in names(LENSES)) {
    d <- copy(pb); d[, y := get(nm)]
    fit <- suppressMessages(lmer(
        y ~ cl + mean_log10_counts + (1 | study) + (1 | donor_id), data = d))
    e <- as.data.table(as.data.frame(emmeans(fit, ~ cl)))
    p <- as.data.table(as.data.frame(pairs(emmeans(fit, ~ cl), adjust = "BH")))
    setnames(e, c("emmean", "SE"), c("emmean", "SE"), skip_absent = TRUE)

    ## The denoised lens is a rate; put it on the log scale the count model uses.
    ## Delta method for both the SE and the CI half-width.
    if (nm == "AGTR1_scvi") {
        if (any(e$emmean <= 0))
            stop("non-positive denoised pseudobulk mean; cannot take a log")
        e[, `:=`(SE = SE / emmean,
                 lower.CL = log(lower.CL), upper.CL = log(upper.CL),
                 emmean = log(emmean))]
        ## A contrast of rates becomes a contrast of logs to first order by
        ## dividing by the mean rate over the levels being compared.
        p[, `:=`(SE = SE / mean(exp(e$emmean)),
                 estimate = estimate / mean(exp(e$emmean)))]
        e[, scale_note := "log AGTR1 per 10^4 transcripts (delta method)"]
        p[, scale_note := "log AGTR1 per 10^4 transcripts (delta method)"]
    } else {
        e[, scale_note := if (nm == "AGTR1_detect") "detected fraction"
                          else "mean log1p-normalised expression"]
        p[, scale_note := if (nm == "AGTR1_detect") "detected fraction"
                          else "mean log1p-normalised expression"]
    }
    ## `n_fitted_*` are the fit-wide totals. Per-CLUSTER support is merged on
    ## below as n_units/n_donors, matching what those two names mean in
    ## agtr1_count_by_cluster.tsv -- panel D joins the two tables, so a column
    ## meaning "per cluster" in one and "whole fit" in the other is a trap.
    e[, `:=`(lens = LENSES[[nm]], level = "pseudobulk",
             n_fitted_units = nrow(pb), n_fitted_donors = uniqueN(pb$donor_id),
             singular = isSingular(fit))]
    p[, `:=`(lens = LENSES[[nm]], level = "pseudobulk",
             n_fitted_units = nrow(pb), n_fitted_donors = uniqueN(pb$donor_id),
             singular = isSingular(fit))]
    emm_all[[nm]] <- e; ph_all[[nm]] <- p
    cat("\n== ", nm, " by cluster (pseudobulk emmeans) ==\n", sep = "")
    print(e[, .(cl, emmean = round(emmean, 4), SE = round(SE, 4))])
}

emm <- rbindlist(emm_all, fill = TRUE)
ph  <- rbindlist(ph_all,  fill = TRUE)

## Per-cluster support, so an underpowered point cannot be read as a firm one.
## The thresholds mirror 10.agtr1_count_models.R's `underpowered` flag.
sup <- pb[, .(n_units = .N, n_donors = uniqueN(donor_id), n_cells = sum(n_cells)),
          by = .(cl)]
sup[, underpowered := n_units < 20L]
emm <- merge(emm, sup, by = "cl", sort = FALSE)
cat("\n== per-cluster support ==\n"); print(sup)

fwrite(emm, file.path(opt$outdir, "agtr1_lens_by_cluster_emmeans.tsv"), sep = "\t")
fwrite(ph,  file.path(opt$outdir, "agtr1_lens_by_cluster_posthoc.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()),
           file.path(opt$outdir, "agtr1_lens_by_cluster_sessionInfo.txt"))
message("done")
