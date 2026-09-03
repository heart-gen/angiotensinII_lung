#!/usr/bin/env Rscript
##
## Is the TGF-beta/BM association TGF-beta signalling, or a generic
## immediate-early / dissociation-stress program?
##
## The reported claim -- TGF-beta response tracks LOWER basement membrane
## (beta = -0.196, BH 0.0074) -- rests on a 17-gene panel whose variance sits
## almost entirely in six genes that are not TGF-beta-specific. This script
## splits that panel into a SMAD-proximal arm and an immediate-early/mechano arm
## and asks which one carries the association.
##
## THE PARAMETER ANSWERING THE QUESTION is the coefficient on
## tgfb_smad_score_z in the head-to-head model, adjusted for tgfb_ieg_score_z.
##
## The arms differ 3.5x in detection by construction, so a null SMAD arm is
## uninterpretable without knowing what a panel of that sparsity can produce.
## 11.tgfb_null_panels.py supplies that: K detection-matched random panels per
## arm, scored identically, fitted here through the identical model.
##
## Design, decision rule and verdict table were fixed in
## _h/TGFB_SPECIFICITY_PLAN.md BEFORE any of this ran.
##
## Outputs (to --outdir/stats_data):
##   tgfb_specificity_models.tsv   observed arms, marginal and head-to-head
##   tgfb_specificity_null.tsv     per-null-panel betas (the empirical null)
##   tgfb_specificity_summary.tsv  empirical p, power, and the verdict
##   tgfb_specificity_logo.tsv     leave-one-gene-out on the full panel
##   tgfb_specificity_varcomp.tsv  between-study variance per arm (Test C)
suppressPackageStartupMessages({
    library(data.table); library(optparse)
    library(lme4); library(lmerTest)
})

## CLI mirrors 04.bm_state_stats.R exactly -- same inputs, same --outdir
## convention, same --min-cells default -- because the betas here only mean
## something in comparison to the ones that script produced.
opt <- parse_args(OptionParser(option_list = list(
    make_option("--bm-meta", type = "character", dest = "bm_meta"),
    make_option("--state-meta", type = "character", dest = "state_meta"),
    make_option("--null-pseudobulk", type = "character", default = NA_character_,
                dest = "null_pb",
                help = "tgfb_null_pseudobulk.tsv.gz from 11.tgfb_null_panels.py"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--ref-beta", type = "double", default = 0.196, dest = "ref_beta",
                help = "reported |beta| the power statement is measured against")
)))

stopifnot(!is.null(opt$bm_meta), !is.null(opt$state_meta), !is.null(opt$outdir))
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---- pseudobulk, built exactly as 04.bm_state_stats.R builds it -------------
## Same unit, same floor, same within-dataset z. If these drift apart the new
## betas are not comparable to the -0.196 they exist to interrogate.
z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

bm <- fread(opt$bm_meta)
st <- fread(opt$state_meta)
setnames(st, 1, "index"); setnames(bm, 1, "index")
drop <- intersect(setdiff(names(bm), "index"), names(st))
d <- merge(st, bm[, .SD, .SDcols = setdiff(names(bm), drop)], by = "index")
message(sprintf("Merged %d cells", nrow(d)))
d[, pericyte_state := factor(pericyte_state)]
if ("dataset" %in% names(d)) d[, dataset := as.character(dataset)]

ARM_COLS <- c("tgfb_response_score", "tgfb_smad_score", "tgfb_ieg_score")
missing_arms <- setdiff(ARM_COLS, names(d))
if (length(missing_arms))
    stop("metadata lacks the specificity arms: ", paste(missing_arms, collapse = ", "),
         ". Re-run 00.bm_score.py after adding tgfb_smad/tgfb_ieg to bm_panels.PANELS.")

logo_cols <- grep("^logo_", names(d), value = TRUE)
score_cols <- c("basement_membrane_score", "fibrillar_collagen_score",
                ARM_COLS, logo_cols)
score_cols <- intersect(score_cols, names(d))

pb <- d[, c(lapply(.SD, mean, na.rm = TRUE),
            .(n_cells = .N, study = first(study), dataset = first(dataset),
              mean_log10_counts = mean(log10_total_counts, na.rm = TRUE))),
        by = .(donor_id, pericyte_state), .SDcols = score_cols]
pb <- pb[n_cells >= opt$min_cells]
message(sprintf("Donor x cluster units (>=%d cells): %d across %d donors",
                opt$min_cells, nrow(pb), uniqueN(pb$donor_id)))

for (cl in score_cols)
    pb[[paste0(cl, "_z")]] <- z_within_dataset(pb[[cl]], pb$dataset)
pb[, bm_minus_fibrillar :=
       basement_membrane_score_z - fibrillar_collagen_score_z]
pb[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]

OUTCOMES <- c("basement_membrane_score_z", "bm_minus_fibrillar_z")

## ---- one fitter, used for the observed arms and for every null panel --------
## Identical formula throughout; that is the whole point of the null.
fit_one <- function(dt, outcome, preds) {
    f <- reformulate(c(preds, "mean_log10_counts",
                       "(1 | study)", "(1 | donor_id)"), response = outcome)
    fit <- try(suppressMessages(lmerTest::lmer(f, data = dt,
                                               control = lmerControl(
                                                   calc.derivs = FALSE))),
               silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    co <- as.data.table(summary(fit)$coefficients, keep.rownames = "term")
    setnames(co, c("term", "estimate", "SE", "df", "t_ratio", "p_value"))
    co <- co[term %in% preds]
    ## Convergence is graded, not assumed: an unconverged fit's coefficient is
    ## not evidence, and with 4000 null fits a silent failure would otherwise
    ## be averaged into the null distribution.
    gmax <- tryCatch(max(abs(fit@optinfo$derivs$gradient)), error = function(e) NA_real_)
    co[, `:=`(outcome = outcome, converged = is.na(gmax) || gmax < 0.01,
              max_grad = gmax)]
    co[]
}

## ---- observed arms ----------------------------------------------------------
spec <- list(
    list(tag = "full_panel",   preds = "tgfb_response_score_z"),
    list(tag = "smad_alone",   preds = "tgfb_smad_score_z"),
    list(tag = "ieg_alone",    preds = "tgfb_ieg_score_z"),
    list(tag = "head_to_head", preds = c("tgfb_smad_score_z", "tgfb_ieg_score_z"))
)
obs <- rbindlist(lapply(spec, function(s)
    rbindlist(lapply(OUTCOMES, function(o) {
        r <- fit_one(pb, o, s$preds)
        if (is.null(r)) return(NULL)
        r[, model := s$tag][]
    }), fill = TRUE)), fill = TRUE)
obs[, n_units := nrow(pb)][, n_donors := uniqueN(pb$donor_id)]
obs[, p_BH := p.adjust(p_value, "BH"), by = model]
fwrite(obs, file.path(outdir, "tgfb_specificity_models.tsv"), sep = "\t")

## ---- Test C: where does each arm's variance live? ---------------------------
## If the IEG arm is a warm-dissociation artifact it should partition markedly
## more variance to `study` (i.e. to protocol) than the SMAD arm does.
vc <- rbindlist(lapply(c("tgfb_response_score_z", "tgfb_smad_score_z",
                         "tgfb_ieg_score_z"), function(v) {
    f <- reformulate(c("1", "(1 | study)", "(1 | dataset)", "(1 | donor_id)"),
                     response = v)
    fit <- try(suppressMessages(lmer(f, data = pb)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    vv <- as.data.table(VarCorr(fit))
    tot <- sum(vv$vcov)
    data.table(score = v, component = vv$grp, variance = vv$vcov,
               pct = 100 * vv$vcov / tot)
}), fill = TRUE)
fwrite(vc, file.path(outdir, "tgfb_specificity_varcomp.tsv"), sep = "\t")

## ---- leave-one-gene-out on the full panel -----------------------------------
if (length(logo_cols)) {
    logo <- rbindlist(lapply(logo_cols, function(cl) {
        r <- fit_one(pb, "basement_membrane_score_z", paste0(cl, "_z"))
        if (is.null(r)) return(NULL)
        r[, dropped_gene := sub("^logo_", "", cl)][]
    }), fill = TRUE)
    setorder(logo, estimate)
    fwrite(logo, file.path(outdir, "tgfb_specificity_logo.tsv"), sep = "\t")
    message(sprintf("Leave-one-out: beta ranges %.4f to %.4f across %d drops",
                    min(logo$estimate), max(logo$estimate), nrow(logo)))
}

## ---- detection-matched empirical null ---------------------------------------
if (!is.na(opt$null_pb) && file.exists(opt$null_pb)) {
    nl <- fread(opt$null_pb)
    null_cols <- grep("^null_", names(nl), value = TRUE)
    pbn <- merge(pb[, .(donor_id, pericyte_state, study, dataset,
                        mean_log10_counts, basement_membrane_score_z,
                        bm_minus_fibrillar_z)],
                 nl, by = c("donor_id", "pericyte_state"))
    if (nrow(pbn) != nrow(pb))
        stop("null pseudobulk does not cover the analysis units: ", nrow(pbn),
             " of ", nrow(pb), ". The null must be fitted on the same units as ",
             "the observed arms or the comparison is meaningless.")
    message(sprintf("Null panels: %d, fitted on %d units", length(null_cols),
                    nrow(pbn)))

    for (cl in null_cols)
        pbn[[paste0(cl, "_z")]] <- z_within_dataset(pbn[[cl]], pbn$dataset)

    nullres <- rbindlist(lapply(seq_along(null_cols), function(i) {
        cl <- null_cols[i]
        if (i %% 200 == 0) message(sprintf("  null fit %d/%d", i, length(null_cols)))
        rbindlist(lapply(OUTCOMES, function(o) {
            r <- fit_one(pbn, o, paste0(cl, "_z"))
            if (is.null(r)) return(NULL)
            r[, `:=`(panel = cl, arm = sub("^null_([a-z]+)_.*$", "\\1", cl))][]
        }), fill = TRUE)
    }), fill = TRUE)
    fwrite(nullres, file.path(outdir, "tgfb_specificity_null.tsv"), sep = "\t")

    ## Empirical p and power, per arm and outcome.
    obs_alone <- obs[model %in% c("smad_alone", "ieg_alone")]
    obs_alone[, arm := fifelse(model == "smad_alone", "smad", "ieg")]

    summ <- rbindlist(lapply(seq_len(nrow(obs_alone)), function(i) {
        row <- obs_alone[i]
        nd <- nullres[arm == row$arm & outcome == row$outcome & converged == TRUE]
        if (!nrow(nd)) return(NULL)
        ## Two-sided empirical p: how often does a detection-matched panel of
        ## random genes reach an effect at least this extreme?
        emp_p <- (1 + sum(abs(nd$estimate) >= abs(row$estimate))) / (1 + nrow(nd))
        ## Power: could a panel this sparse have produced the reported effect?
        pwr <- mean(abs(nd$estimate) >= opt$ref_beta)
        data.table(
            arm = row$arm, outcome = row$outcome,
            beta_observed = row$estimate, p_model = row$p_value,
            n_null = nrow(nd),
            null_mean = mean(nd$estimate), null_sd = sd(nd$estimate),
            null_q025 = quantile(nd$estimate, 0.025),
            null_q975 = quantile(nd$estimate, 0.975),
            empirical_p = emp_p,
            frac_null_reaching_ref = pwr, ref_beta = opt$ref_beta)
    }), fill = TRUE)

    ## The verdict table from TGFB_SPECIFICITY_PLAN.md section 7, applied
    ## mechanically so the rule cannot drift after seeing the numbers.
    summ[, outside_null := empirical_p < 0.05]
    summ[, well_powered := frac_null_reaching_ref >= 0.80]
    summ[, arm_verdict := fifelse(
        outside_null & beta_observed < 0, "carries the association",
        fifelse(!outside_null & well_powered, "null, adequately powered",
                "inside null, underpowered -- uninformative"))]
    fwrite(summ, file.path(outdir, "tgfb_specificity_summary.tsv"), sep = "\t")

    message("\n---- verdict inputs ----")
    print(summ[, .(arm, outcome, beta_observed, empirical_p,
                   frac_null_reaching_ref, arm_verdict)])
} else {
    warning("no --null-pseudobulk supplied; the arms are reported WITHOUT their ",
            "detection-matched null, and a null SMAD arm cannot be interpreted. ",
            "Run 11.tgfb_null_panels.py first.", call. = FALSE)
}

h2h <- obs[model == "head_to_head" & outcome == "basement_membrane_score_z"]
if (nrow(h2h)) {
    message("\n---- head-to-head (the parameter that answers the question) ----")
    print(h2h[, .(term, estimate, SE, t_ratio, p_value, p_BH)])
}

writeLines(capture.output(sessionInfo()),
           file.path(outdir, "tgfb_specificity_sessionInfo.txt"))
