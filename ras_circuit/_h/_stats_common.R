## Shared statistics helpers for the ras_circuit module.
##
## Copied VERBATIM from the modules that established them, so that ras_circuit
## estimates are on the same scale and carry the same sign convention:
##   flip_contrast()     agt_axis/_h/01.ras_landscape_stats.R:39
##   z_within_dataset()  basement_membrane/_h/10.agtr1_count_models.R:56
##   tidy_row()          basement_membrane/_h/10.agtr1_count_models.R:99
## If one of those changes, change this file too.
##
## Source with: source("../_h/_stats_common.R") (scripts run from ras_circuit/_m).

suppressPackageStartupMessages({
    library(data.table)
    library(lme4)
})

## emmeans' trt.vs.ctrl returns (other - REF); tables here report (REF - other),
## so the WHOLE contrast is reversed -- estimate, t/z ratio and label together.
flip_contrast <- function(ct, ref) {
    ct$estimate <- -ct$estimate
    if ("t.ratio" %in% names(ct)) ct$t.ratio <- -ct$t.ratio
    if ("z.ratio" %in% names(ct)) ct$z.ratio <- -ct$z.ratio
    ct$contrast <- paste0(ref, " - ",
                          sub(" - .*$", "", as.character(ct$contrast)))
    ct
}

## Identical to basement_membrane/_h/{04,10}: centre and scale within dataset.
z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

## One coefficient row from an lme4/lmerTest fit. `term` defaults to `pred`.
tidy_row <- function(fit, pred, level, model, n, n_don, note = "", spec = "primary",
                     term = NULL) {
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients
    term <- if (is.null(term)) pred else term
    if (!term %in% rownames(co)) return(NULL)
    co <- co[term, , drop = FALSE]; rownames(co) <- pred
    pcol <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(co))
    data.table(level = level, model = model, spec = spec, predictor = pred,
               estimate = co[pred, "Estimate"], SE = co[pred, "Std. Error"],
               z_value = co[pred, 3],
               p_value = if (length(pcol)) co[pred, pcol[1]] else NA_real_,
               n_units = n, n_donors = n_don,
               singular = isSingular(fit), note = note)
}

write_tsv_safe <- function(x, file) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    fwrite(as.data.table(x), file, sep = "\t", na = "NA", quote = FALSE)
}

## Refuse to run on a missing input rather than silently skipping a block
## (same contract as figures/_h/assemble_mechanism_figures.R::read_req).
read_req <- function(path, ...) {
    if (!file.exists(path))
        stop("missing input: ", path, "\n  Run the producing step first.", call. = FALSE)
    fread(path, ...)
}

## Empirical p of an observed estimate against its matched null, referred to the
## null's OWN centre (memory: gene-set-score-null-not-zero). Two-sided on the
## centred scale; the floor is 1/(N+1).
emp_p <- function(obs, null) {
    null <- null[is.finite(null)]
    if (!length(null) || !is.finite(obs)) return(NA_real_)
    c0 <- mean(null)
    (1 + sum(abs(null - c0) >= abs(obs - c0))) / (length(null) + 1)
}
