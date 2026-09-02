## Does the matrix program track AGTR1 abundance, with dropout inside the model?
##
## WHY THIS LENS EXISTS. Every other AGTR1 readout in this project is a derived
## quantity: `AGTR1_expr`/`AGTR1_detect` come from the soupX ambient-corrected
## matrix, and `AGTR1_scvi` is a model output that has to earn its keep through
## the validity gate in localization/airspace_analysis/_h/02.qc_airspace.py. This
## script asks the same association question WITHOUT imputing anything. AGTR1
## integer counts are the RESPONSE of a negative-binomial GLMM with a library
## size offset, so sampling zeros are part of the likelihood rather than
## something to be removed beforehand, and nothing is borrowed across genes --
## which is the circularity that any denoiser risks when the predictor is itself
## a multi-gene score.
##
## ROLE REVERSAL, and what it does and does not change. In 04.bm_state_stats.R
## the matrix score is the outcome and AGTR1 the predictor. Here the matrix score
## is the predictor and AGTR1 the response, because the measurement error that
## motivated the whole three-lens design is on the AGTR1 side and a GLM can only
## model it there. The null tested is the same -- no association between the
## matrix program and AGTR1 abundance -- and a standardized slope is still
## reported, but this is an association test, not a direction-of-effect test.
## Positive beta = a one-SD higher matrix score goes with higher AGTR1 per
## transcript sampled.
##
## THREE MODELS, deliberately:
##   (1) pseudobulk NB GLMM (primary) -- the donor x cluster unit used everywhere
##       else in the module, so the estimate is comparable to bm_vs_agtr1_models.tsv;
##   (2) pseudobulk binomial GLMM on detected/undetected cells -- makes no
##       distributional assumption about counts beyond independence, so it is the
##       robustness check on (1);
##   (3) cell-level NB GLMM -- the highest-resolution version, reported second
##       because pseudoreplication defenses matter more than resolution here.
##
## Depth enters twice on purpose: `offset(log(raw_total_counts))` fixes the AGTR1
## exposure, and `mean_log10_counts` is a free covariate because the matrix
## SCORE is itself depth-dependent (BM rises with depth). The offset handles one
## side of the association, the covariate the other.
##
## Input:  agtr1_count_input.tsv.gz (09.agtr1_counts_prep.py)
## Output: agtr1_count_models.tsv, agtr1_count_README.txt
suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(splines)
    library(emmeans)
})

option_list <- list(
    make_option("--input", type = "character", default = "./agtr1_count_input.tsv.gz"),
    make_option("--outdir", type = "character", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--seed", type = "integer", default = 13L)
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

## Identical to 04.bm_state_stats.R so the two are on one scale
z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

dl <- fread(opt$input)
message(sprintf("cells: %d | donors: %d | AGTR1 detected in %.1f%%",
                nrow(dl), uniqueN(dl$donor_id), 100 * mean(dl$AGTR1_count > 0)))

## ---- pseudobulk: donor x Leiden cluster ---------------------------------
## `pericyte_state` and not `state_program`: the latter is assigned by a
## relative-enrichment argmax that includes the BM panel, so it is not
## independent of the predictor (same reasoning as 04.bm_state_stats.R).
pb <- dl[, .(n_cells = .N,
             agtr1_sum = sum(AGTR1_count),
             n_pos = sum(AGTR1_count > 0),
             total_sum = sum(raw_total_counts),
             mean_log10_counts = mean(log10_total_counts),
             basement_membrane_score = mean(basement_membrane_score),
             fibrillar_collagen_score = mean(fibrillar_collagen_score),
             study = study[1], dataset = dataset[1]),
         by = .(donor_id, pericyte_state)]
pb <- pb[n_cells >= opt$min_cells]
message(sprintf("pseudobulk units (>=%d cells): %d from %d donors",
                opt$min_cells, nrow(pb), uniqueN(pb$donor_id)))

for (s in c("basement_membrane_score", "fibrillar_collagen_score"))
    pb[[paste0(s, "_z")]] <- z_within_dataset(pb[[s]], pb$dataset)
pb[, bm_minus_fibrillar := basement_membrane_score - fibrillar_collagen_score]
pb[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]
pb[, bm_resid_collagen := residuals(suppressMessages(lmer(
    basement_membrane_score_z ~ fibrillar_collagen_score_z +
        mean_log10_counts + (1 | study) + (1 | donor_id), data = pb)))]

PREDICTORS <- c("basement_membrane_score_z", "fibrillar_collagen_score_z",
                "bm_minus_fibrillar_z", "bm_resid_collagen")

## ---- model fitting ------------------------------------------------------
tidy_row <- function(fit, pred, level, model, n, n_don, note = "", spec = "primary",
                     term = NULL) {
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients
    term <- if (is.null(term)) pred else term
    if (!term %in% rownames(co)) return(NULL)
    co <- co[term, , drop = FALSE]; rownames(co) <- pred
    data.table(level = level, model = model, spec = spec, predictor = pred,
               estimate = co[pred, "Estimate"], SE = co[pred, "Std. Error"],
               z_value = co[pred, 3], p_value = co[pred, 4],
               n_units = n, n_donors = n_don,
               singular = isSingular(fit), note = note)
}

## NB GLMM, with a Poisson + observation-level random effect fallback. The
## fallback is recorded in `model` so it can never be mistaken for the NB fit.
fit_nb <- function(form_nb, form_pois, data, pred, level, n, n_don) {
    fit <- tryCatch(suppressMessages(glmer.nb(form_nb, data = data)),
                    error = function(e) NULL, warning = function(w) NULL)
    if (!is.null(fit)) return(tidy_row(fit, pred, level, "NB GLMM", n, n_don))
    message("  glmer.nb failed for ", pred, " at ", level,
            "; falling back to Poisson + OLRE")
    data <- copy(data)[, olre := factor(seq_len(.N))]
    fit <- tryCatch(suppressMessages(glmer(form_pois, data = data,
                                           family = poisson)),
                    error = function(e) NULL)
    tidy_row(fit, pred, level, "Poisson+OLRE", n, n_don,
             note = "glmer.nb did not converge")
}

res <- list()
n_pb <- nrow(pb); n_don_pb <- uniqueN(pb$donor_id)
for (p in PREDICTORS) {
    ## (1) primary: counts with a library-size offset
    f_nb <- as.formula(sprintf(
        "agtr1_sum ~ %s + mean_log10_counts + (1|study) + (1|donor_id) + offset(log(total_sum))", p))
    f_po <- as.formula(sprintf(
        "agtr1_sum ~ %s + mean_log10_counts + (1|study) + (1|donor_id) + (1|olre) + offset(log(total_sum))", p))
    res[[length(res) + 1]] <- fit_nb(f_nb, f_po, pb, p, "pseudobulk", n_pb, n_don_pb)

    ## (2) robustness: detected / not detected, no count distribution assumed
    f_bin <- as.formula(sprintf(
        "cbind(n_pos, n_cells - n_pos) ~ %s + mean_log10_counts + (1|study) + (1|donor_id)", p))
    fit <- tryCatch(suppressMessages(glmer(f_bin, data = pb, family = binomial)),
                    error = function(e) NULL)
    res[[length(res) + 1]] <- tidy_row(fit, p, "pseudobulk", "binomial GLMM",
                                       n_pb, n_don_pb)
}

## (3) cell level. Scores are z-standardized within dataset on the cell table so
## the slope is on the same scale as the pseudobulk one.
for (s in c("basement_membrane_score", "fibrillar_collagen_score"))
    dl[[paste0(s, "_z")]] <- z_within_dataset(dl[[s]], dl$dataset)
dl[, bm_minus_fibrillar := basement_membrane_score - fibrillar_collagen_score]
dl[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]
n_cl <- nrow(dl); n_don_cl <- uniqueN(dl$donor_id)
for (p in intersect(PREDICTORS, names(dl))) {
    f_nb <- as.formula(sprintf(
        "AGTR1_count ~ %s + log10_total_counts + (1|study) + (1|donor_id) + offset(log(raw_total_counts))", p))
    f_po <- as.formula(sprintf(
        "AGTR1_count ~ %s + log10_total_counts + (1|study) + (1|donor_id) + (1|olre) + offset(log(raw_total_counts))", p))
    res[[length(res) + 1]] <- fit_nb(f_nb, f_po, dl, p, "cell", n_cl, n_don_cl)
}

## ---- sensitivities ------------------------------------------------------
## The primary fits could be produced by three things other than biology, and
## each has a specific antidote:
##   ambient   -- AGTR1 in the soup would inflate counts in cells that also
##                carry more matrix transcripts. raw/X is NOT ambient-corrected,
##                so the tracer score enters as a covariate.
##   depth     -- a linear term in log10 depth may not absorb a nonlinear
##                relationship; a natural spline is the stricter control.
##   donor     -- a between-donor association can be driven by donor-level
##                confounding. The Mundlak decomposition splits the predictor
##                into a donor mean and a within-donor deviation; the deviation
##                coefficient is the within-donor estimate, which no
##                time-invariant donor characteristic can produce.
SENS <- c("basement_membrane_score_z", "bm_minus_fibrillar_z")

add_mundlak <- function(d, p) {
    d <- copy(d)
    d[, dmean := mean(get(p)), by = donor_id]
    d[, dev := get(p) - dmean]
    d
}

for (p in SENS) {
    ## cell-level NB (the specification that converged natively)
    if (p %in% names(dl)) {
        f <- as.formula(sprintf("AGTR1_count ~ %s + log10_total_counts + ambient_tracer_score + (1|study) + (1|donor_id) + offset(log(raw_total_counts))", p))
        fit <- tryCatch(suppressMessages(glmer.nb(f, data = dl)), error = function(e) NULL,
                        warning = function(w) NULL)
        res[[length(res)+1]] <- tidy_row(fit, p, "cell", "NB GLMM", n_cl, n_don_cl,
                                         spec = "+ambient")
        f <- as.formula(sprintf("AGTR1_count ~ %s + ns(log10_total_counts, 3) + (1|study) + (1|donor_id) + offset(log(raw_total_counts))", p))
        fit <- tryCatch(suppressMessages(glmer.nb(f, data = dl)), error = function(e) NULL,
                        warning = function(w) NULL)
        res[[length(res)+1]] <- tidy_row(fit, p, "cell", "NB GLMM", n_cl, n_don_cl,
                                         spec = "depth spline")
        dm <- add_mundlak(dl, p)
        f <- as.formula("AGTR1_count ~ dev + dmean + log10_total_counts + (1|study) + (1|donor_id) + offset(log(raw_total_counts))")
        fit <- tryCatch(suppressMessages(glmer.nb(f, data = dm)), error = function(e) NULL,
                        warning = function(w) NULL)
        res[[length(res)+1]] <- tidy_row(fit, p, "cell", "NB GLMM", n_cl, n_don_cl,
                                         spec = "within-donor (Mundlak)", term = "dev")
    }
    ## pseudobulk binomial (the specification that makes no count assumption)
    f <- as.formula(sprintf("cbind(n_pos, n_cells - n_pos) ~ %s + mean_log10_counts + (1|study) + (1|donor_id)", p))
    dmp <- add_mundlak(pb, p)
    fit <- tryCatch(suppressMessages(glmer(as.formula("cbind(n_pos, n_cells - n_pos) ~ dev + dmean + mean_log10_counts + (1|study) + (1|donor_id)"),
                                           data = dmp, family = binomial)),
                    error = function(e) NULL)
    res[[length(res)+1]] <- tidy_row(fit, p, "pseudobulk", "binomial GLMM", n_pb,
                                     n_don_pb, spec = "within-donor (Mundlak)",
                                     term = "dev")
    f <- as.formula(sprintf("cbind(n_pos, n_cells - n_pos) ~ %s + ns(mean_log10_counts, 3) + (1|study) + (1|donor_id)", p))
    fit <- tryCatch(suppressMessages(glmer(f, data = pb, family = binomial)),
                    error = function(e) NULL)
    res[[length(res)+1]] <- tidy_row(fit, p, "pseudobulk", "binomial GLMM", n_pb,
                                     n_don_pb, spec = "depth spline")
}

## ---- AGTR1 abundance across pericyte clusters ---------------------------
## This is the count-model replacement for a by-cluster comparison of raw
## detection, and it is not the same question. Detection rate and library size
## track each other across these clusters (46.1% / 2510 counts in cluster 0 down
## to 20.9% / 1323 in cluster 4), so a detection comparison is substantially a
## depth comparison; the offset removes that.
##
## NOT CIRCULAR: the Leiden clusters come from X_pca_harmony, built on 2,000 HVGs
## that do NOT include AGTR1 (highly_variable = FALSE, dispersions_norm = 0.081,
## highly_variable_intersection = FALSE in pericyte_states.h5ad). AGTR1 therefore
## played no part in defining the groups being compared -- a stronger guarantee
## than the BM score has against state_program.
##
## TWO ESTIMANDS, both reported, because they answer different questions:
##   with_offset  -- AGTR1 per transcript sampled (a concentration). Library size
##                   varies ~2-fold across these clusters, so a cluster with
##                   globally lower output can show a higher AGTR1 FRACTION while
##                   holding the same absolute number of AGTR1 molecules.
##   depth_covar  -- no offset, log10 depth as a free covariate: closer to a
##                   per-cell abundance. Divergence between the two is
##                   interpretable, not an error.
## `cl` is pericyte_state as a FACTOR. The column is stored as an integer, and
## passing it to the formula bare fits cluster as a LINEAR covariate -- the model
## then estimates a slope across cluster labels, which is meaningless, and
## emmeans returns one marginal mean at the covariate average rather than six.
## That is what the first run of this block did; never use the bare column here.
## Fitting strategy. glmer.nb converges for the continuous-predictor models above
## but NOT once cluster enters as a 5-df factor: every by-cluster NB fit failed,
## and the Poisson + OLRE fallbacks it dropped to were themselves non-converged
## (max|grad| up to 0.34 against a 0.002 tolerance, with emmeans SEs collapsing to
## ~9e-4 -- an artefact, not precision). Three changes fix that:
##   1. a warning no longer discards a fit. The old tryCatch(warning = NULL) threw
##      away usable fits that merely emitted a convergence note; fits are now
##      graded on max|grad| instead of on whether lme4 was chatty.
##   2. if glmer.nb's joint theta search fails, theta is estimated once from the
##      fixed-effect-only MASS::glm.nb and held fixed while lme4 fits the random
##      effects. This is still a negative binomial, not a Poisson.
##   3. Poisson + OLRE remains, but only as the last resort, and whichever fit is
##      used carries its max|grad| and a converged flag into the output so a
##      non-converged estimate cannot be read as a clean one.
MAX_GRAD_OK <- 0.01

fit_keep_warnings <- function(expr) {
    w <- character()
    val <- withCallingHandlers(
        tryCatch(suppressMessages(expr), error = function(e) NULL),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
        })
    list(fit = val, warns = w)
}

max_grad <- function(fit) {
    g <- tryCatch(fit@optinfo$derivs$gradient, error = function(e) NULL)
    if (is.null(g) || !length(g)) return(NA_real_)
    max(abs(g))
}

## `group_var` generalizes this block. It is called twice:
##   pericyte_state -- the Leiden clusters. AGTR1 is not among the 2,000 HVGs
##                     that define them, so the grouping is AGTR1-independent.
##   state_program  -- the marker-panel argmax label. This is NOT independent of
##                     the BM score (which is one of the panels), so it must
##                     never be used for a BM-vs-AGTR1 test. It IS independent of
##                     AGTR1, which is the only guarantee an AGTR1-by-program
##                     comparison needs, and that comparison is what
##                     figureS_acta2_control panel A asks.
by_group <- function(data, resp_expr, depth_var, offset_expr, level, spec,
                     group_var = "pericyte_state") {
    data <- copy(data)[, cl := factor(get(group_var))]
    rhs <- sprintf("cl + %s", depth_var)
    f_full  <- as.formula(sprintf("%s ~ %s + (1|study) + (1|donor_id)%s",
                                  resp_expr, rhs, offset_expr))
    f_fixed <- as.formula(sprintf("%s ~ %s%s", resp_expr, rhs, offset_expr))
    f_olre  <- as.formula(sprintf("%s ~ %s + (1|study) + (1|donor_id) + (1|olre)%s",
                                  resp_expr, rhs, offset_expr))
    ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

    cands <- list()
    if (grepl("cbind", resp_expr)) {
        r <- fit_keep_warnings(glmer(f_full, data = data, family = binomial,
                                     control = ctrl))
        if (!is.null(r$fit)) cands[["binomial GLMM"]] <- r$fit
    } else {
        r <- fit_keep_warnings(glmer.nb(f_full, data = data, control = ctrl))
        if (!is.null(r$fit)) cands[["NB GLMM"]] <- r$fit

        th <- fit_keep_warnings(MASS::glm.nb(f_fixed, data = data))$fit
        if (!is.null(th)) {
            r <- fit_keep_warnings(glmer(f_full, data = data, control = ctrl,
                                         family = MASS::negative.binomial(theta = th$theta)))
            if (!is.null(r$fit)) cands[["NB GLMM (fixed theta)"]] <- r$fit
        }

        data[, olre := factor(seq_len(.N))]
        r <- fit_keep_warnings(glmer(f_olre, data = data, family = poisson,
                                     control = ctrl))
        if (!is.null(r$fit)) cands[["Poisson+OLRE"]] <- r$fit
    }
    if (!length(cands)) {
        message("  by_group fit failed: ", level, " / ", spec); return(NULL)
    }

    ## Prefer the first strategy that converges; otherwise take the least-bad and
    ## say so. Order of cands is the preference order.
    grads <- vapply(cands, max_grad, numeric(1))
    good  <- which(is.finite(grads) & grads < MAX_GRAD_OK)
    idx   <- if (length(good)) good[1] else which.min(ifelse(is.finite(grads), grads, Inf))
    fit   <- cands[[idx]]; mdl <- names(cands)[idx]; mg <- grads[idx]
    converged <- isTRUE(is.finite(mg) && mg < MAX_GRAD_OK)
    if (!converged) {
        message("  by_group NOT CONVERGED: ", level, " / ", spec,
                " (best = ", mdl, ", max|grad| = ", signif(mg, 3), ")")
    } else if (idx > 1) {
        message("  by_group: ", level, " / ", spec, " used ", mdl,
                " (max|grad| = ", signif(mg, 3), ")")
    }

    ## offset = 0 so the marginal means are log rates per unit exposure rather
    ## than at the mean library size of the data
    emm <- suppressMessages(emmeans(fit, ~ cl, offset = 0))
    em <- as.data.table(summary(emm))
    setnames(em, names(em)[2:3], c("estimate", "SE"))
    ## A per-10k rate is only defined when a library-size offset fixes the
    ## exposure; without one, exp(estimate) is on an arbitrary scale.
    ## The rate scale follows WHICH offset was used: a transcript-library offset
    ## gives AGTR1 per 10k transcripts, a cell-count offset gives AGTR1 per cell.
    ## Reporting both as "per 10k" would silently mix two different quantities.
    is_bin <- grepl("cbind", resp_expr)
    ru <- if (is_bin) NA_character_
          else if (grepl("n_cells", offset_expr)) "AGTR1 per cell"
          else if (nzchar(offset_expr)) "AGTR1 per 10k transcripts"
          else NA_character_
    ## The Leiden labels are integers; state_program labels are strings. Keep the
    ## column named after the grouping actually used so the two output tables
    ## cannot be confused for one another.
    grp <- as.character(em$cl)
    em[, (group_var) := if (group_var == "pericyte_state") as.integer(grp) else grp]
    em[, `:=`(grouping = group_var,
              level = level, spec = spec, model = mdl,
              max_grad = mg, converged = converged, rate_unit = ru,
              rate = if (is.na(ru)) NA_real_
                     else if (ru == "AGTR1 per cell") exp(estimate)
                     else exp(estimate) * 1e4)]
    em[, cl := NULL]
    ct <- as.data.table(summary(contrast(emm, "pairwise", adjust = "BH")))
    ct[, `:=`(grouping = group_var, level = level, spec = spec, model = mdl,
              max_grad = mg, converged = converged)]
    list(emmeans = em, posthoc = ct, singular = isSingular(fit))
}

## Power is very unequal and two clusters cannot support a claim: cluster 5 is
## 4 donors / 113 cells / 46 AGTR1 counts, cluster 4 is 13 donors. Every contrast
## involving cluster 5 rests on 4 paired donors. Emitted as a column so the
## caller cannot forget.
cl_n <- pb[, .(n_units = .N, n_donors = uniqueN(donor_id), n_cells = sum(n_cells),
               agtr1_counts = sum(agtr1_sum)), by = pericyte_state]
cl_n[, underpowered := n_donors < 20]
setorder(cl_n, pericyte_state)
print(cl_n)

bc_em <- list(); bc_ph <- list()
specs <- list(
    list(d = pb, r = "agtr1_sum", dep = "mean_log10_counts",
         off = " + offset(log(total_sum))", lev = "pseudobulk", sp = "with_offset"),
    ## offset(log(n_cells)), NOT a bare intercept. A pseudobulk unit's AGTR1 total
    ## scales with how many cells it pools, and unit size is wildly unequal across
    ## clusters (90 cells/unit in cluster 0 vs 16 in cluster 3). With no exposure
    ## term the fit recovers log(90/16) = 1.72 -- which is almost exactly the 1.62
    ## "cluster 0 >> cluster 3" estimate the first run produced. That was unit size,
    ## not biology. Pooling cells is the exposure here; depth stays a free covariate,
    ## which is what makes this the per-cell-abundance estimand.
    list(d = pb, r = "agtr1_sum", dep = "mean_log10_counts",
         off = " + offset(log(n_cells))", lev = "pseudobulk", sp = "depth_covar"),
    list(d = pb, r = "cbind(n_pos, n_cells - n_pos)", dep = "mean_log10_counts",
         off = "", lev = "pseudobulk", sp = "detection"),
    list(d = pb, r = "agtr1_sum", dep = "ns(mean_log10_counts, 3)",
         off = " + offset(log(total_sum))", lev = "pseudobulk", sp = "depth spline"),
    list(d = dl, r = "AGTR1_count", dep = "log10_total_counts",
         off = " + offset(log(raw_total_counts))", lev = "cell", sp = "with_offset"),
    list(d = dl, r = "AGTR1_count", dep = "log10_total_counts + ambient_tracer_score",
         off = " + offset(log(raw_total_counts))", lev = "cell", sp = "+ambient")
)
for (s in specs) {
    r <- by_group(s$d, s$r, s$dep, s$off, s$lev, s$sp,
                  group_var = "pericyte_state")
    if (is.null(r)) next
    bc_em[[length(bc_em) + 1]] <- r$emmeans
    bc_ph[[length(bc_ph) + 1]] <- r$posthoc
}
if (length(bc_em)) {
    bc_em <- merge(rbindlist(bc_em, fill = TRUE), cl_n, by = "pericyte_state",
                   all.x = TRUE)
    setorder(bc_em, spec, level, pericyte_state)
    fwrite(bc_em, file.path(opt$outdir, "agtr1_count_by_cluster.tsv"), sep = "\t")
    bc_ph <- rbindlist(bc_ph, fill = TRUE)
    ## Flag any contrast that leans on an underpowered cluster
    low <- cl_n[underpowered == TRUE, as.character(pericyte_state)]
    bc_ph[, involves_underpowered := grepl(paste(low, collapse = "|"), contrast)]
    setorder(bc_ph, spec, level, contrast)
    fwrite(bc_ph, file.path(opt$outdir, "agtr1_count_by_cluster_posthoc.tsv"), sep = "\t")
    message("wrote agtr1_count_by_cluster{,_posthoc}.tsv")
    print(bc_em[spec == "with_offset" & level == "pseudobulk",
                .(pericyte_state, estimate = round(estimate, 3),
                  rate = round(rate, 3), rate_unit, n_donors, underpowered)])
}

## ---- AGTR1 abundance across state PROGRAMS ------------------------------
## Same machinery, different grouping, for one specific question:
## figureS_acta2_control asks whether AGTR1 is reducible to ACTA2+ contractile
## identity. That question is posed over `state_program`, so the arbiter for it
## has to be estimated over `state_program` too -- reading it off the by-cluster
## table would mean comparing a claim to a different grouping's numbers.
##
## WHY THIS IS ADMISSIBLE HERE AND NOT ABOVE. `state_program` is a
## relative-enrichment argmax over the marker panels, one of which is the BM
## panel; it is therefore NOT independent of the BM score, and using it for a
## BM-vs-AGTR1 test would be circular (04.bm_state_stats.R says the same). AGTR1
## is in no panel, so the grouping is independent of the RESPONSE, which is the
## only independence an AGTR1-by-program comparison requires. The circularity is
## on the predictor side, and there is no matrix predictor in this fit.
##
## Programs with too few donors cannot support a contrast; the flag travels with
## the table rather than living in a comment.
pb_prog <- dl[, .(n_cells = .N,
                  agtr1_sum = sum(AGTR1_count),
                  n_pos = sum(AGTR1_count > 0),
                  total_sum = sum(raw_total_counts),
                  mean_log10_counts = mean(log10_total_counts),
                  study = study[1], dataset = dataset[1]),
              by = .(donor_id, state_program)]
pb_prog <- pb_prog[n_cells >= opt$min_cells]
message(sprintf("by-program pseudobulk units (>=%d cells): %d from %d donors",
                opt$min_cells, nrow(pb_prog), uniqueN(pb_prog$donor_id)))

pg_n <- pb_prog[, .(n_units = .N, n_donors = uniqueN(donor_id),
                    n_cells = sum(n_cells), agtr1_counts = sum(agtr1_sum)),
                by = state_program]
pg_n[, underpowered := n_donors < 20]
setorder(pg_n, state_program)
print(pg_n)

pg_em <- list(); pg_ph <- list()
prog_specs <- list(
    list(d = pb_prog, r = "agtr1_sum", dep = "mean_log10_counts",
         off = " + offset(log(total_sum))", lev = "pseudobulk", sp = "with_offset"),
    list(d = pb_prog, r = "agtr1_sum", dep = "mean_log10_counts",
         off = " + offset(log(n_cells))", lev = "pseudobulk", sp = "depth_covar"),
    list(d = pb_prog, r = "cbind(n_pos, n_cells - n_pos)", dep = "mean_log10_counts",
         off = "", lev = "pseudobulk", sp = "detection"),
    list(d = dl, r = "AGTR1_count", dep = "log10_total_counts",
         off = " + offset(log(raw_total_counts))", lev = "cell", sp = "with_offset"),
    list(d = dl, r = "AGTR1_count", dep = "log10_total_counts + ambient_tracer_score",
         off = " + offset(log(raw_total_counts))", lev = "cell", sp = "+ambient")
)
for (s in prog_specs) {
    r <- by_group(s$d, s$r, s$dep, s$off, s$lev, s$sp, group_var = "state_program")
    if (is.null(r)) next
    pg_em[[length(pg_em) + 1]] <- r$emmeans
    pg_ph[[length(pg_ph) + 1]] <- r$posthoc
}
if (length(pg_em)) {
    pg_em <- merge(rbindlist(pg_em, fill = TRUE), pg_n, by = "state_program",
                   all.x = TRUE)
    setorder(pg_em, spec, level, state_program)
    fwrite(pg_em, file.path(opt$outdir, "agtr1_count_by_program.tsv"), sep = "\t")
    pg_ph <- rbindlist(pg_ph, fill = TRUE)
    low <- pg_n[underpowered == TRUE, as.character(state_program)]
    pg_ph[, involves_underpowered := if (length(low))
              grepl(paste(low, collapse = "|"), contrast) else FALSE]
    setorder(pg_ph, spec, level, contrast)
    fwrite(pg_ph, file.path(opt$outdir, "agtr1_count_by_program_posthoc.tsv"),
           sep = "\t")
    message("wrote agtr1_count_by_program{,_posthoc}.tsv")
    print(pg_em[spec == "with_offset" & level == "pseudobulk",
                .(state_program, estimate = round(estimate, 3),
                  rate = round(rate, 3), rate_unit, n_donors, underpowered)])
}

## Is cluster confounded with disease or study? Reported, not adjusted for --
## adjusting would change the estimand.
if ("lung_condition" %in% names(dl)) {
    tab <- dl[, .N, by = .(pericyte_state, lung_condition)]
    fwrite(dcast(tab, pericyte_state ~ lung_condition, value.var = "N", fill = 0),
           file.path(opt$outdir, "agtr1_count_cluster_composition.tsv"), sep = "\t")
}

out <- rbindlist(res, fill = TRUE)
## One BH family per (level, model, spec): these are the same correlated matrix
## contrasts asked several ways, not independent questions.
out[, p_BH := p.adjust(p_value, "BH"), by = .(level, model, spec)]
setorder(out, spec, level, model, predictor)
fwrite(out, file.path(opt$outdir, "agtr1_count_models.tsv"), sep = "\t")
print(out[, .(spec, level, model, predictor, estimate = round(estimate, 4),
              SE = round(SE, 4), p_value = signif(p_value, 3),
              p_BH = signif(p_BH, 3))])

writeLines(c(
    "AGTR1 count-model lens -- no imputation",
    strrep("=", 60), "",
    "AGTR1 integer counts (raw/X, not the soupX float matrix) are the RESPONSE of",
    "a negative-binomial GLMM with offset(log(library size)) and (1|study) +",
    "(1|donor_id). Dropout is therefore inside the likelihood, and nothing is",
    "borrowed across genes, so this lens cannot manufacture or erase an",
    "association the way a denoiser can.",
    "",
    "ROLE REVERSAL: the matrix score is the predictor here and the outcome in",
    "bm_vs_agtr1_models.tsv. The null is the same; the sign means 'one SD higher",
    "matrix score goes with higher AGTR1 per transcript sampled'. This is an",
    "association test, not a direction-of-effect test.",
    "",
    "Depth appears twice by design: the offset fixes the AGTR1 exposure, and",
    "mean_log10_counts is free because the matrix score is itself depth-dependent.",
    "",
    sprintf("Pseudobulk units: %d (>= %d cells) from %d donors.", n_pb,
            opt$min_cells, n_don_pb),
    sprintf("Cells: %d from %d donors; AGTR1 detected in %.1f%%.", n_cl, n_don_cl,
            100 * mean(dl$AGTR1_count > 0)),
    "",
    "Rows marked model = 'Poisson+OLRE' are fallbacks where glmer.nb did not",
    "converge; they are not NB fits and should be reported as such.",
    "BH is applied within (level, model).",
    "",
    "GROUP COMPARISONS -- two groupings, two files, different admissibility:",
    "  agtr1_count_by_cluster{,_posthoc}.tsv  grouping = pericyte_state (Leiden).",
    "    AGTR1 is not among the 2,000 HVGs that define these clusters, so the",
    "    grouping is independent of both predictor and response. This is the",
    "    arbiter for any AGTR1-across-clusters claim.",
    "  agtr1_count_by_program{,_posthoc}.tsv  grouping = state_program.",
    "    state_program is a marker-panel argmax that INCLUDES the BM panel, so it",
    "    is not independent of the matrix score and must never arbitrate a",
    "    BM-vs-AGTR1 test. AGTR1 is in no panel, so it is admissible for the",
    "    AGTR1-across-programs question only (figureS_acta2_control panel A).",
    "Both carry `converged` and `max_grad`; an unconverged row is not a result."
), file.path(opt$outdir, "agtr1_count_README.txt"))
message("wrote ", file.path(opt$outdir, "agtr1_count_models.tsv"))
