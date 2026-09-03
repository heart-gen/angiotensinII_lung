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

## Leave-one-out scores are NOT in the cell table: 11.tgfb_null_panels.py
## re-scores them and writes them, already aggregated, into the null pseudobulk
## alongside the null panels. They are picked up in the null section below.
score_cols <- c("basement_membrane_score", "fibrillar_collagen_score", ARM_COLS)
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

## ---- sensitivity 4 of the plan: complexity floor ----------------------------
## Sparse panel scores are least reliable in shallow units, and the SMAD arm is
## the sparse one -- so if its estimate is an artifact of low-depth noise it
## should change when the shallow half is dropped. Listed in
## TGFB_SPECIFICITY_PLAN.md section 8; run here so the pre-specified set is
## complete rather than quietly trimmed.
depth_cut <- median(pb$mean_log10_counts, na.rm = TRUE)
pb_deep <- pb[mean_log10_counts >= depth_cut]
message(sprintf("Complexity floor: %d of %d units at/above median depth %.3f",
                nrow(pb_deep), nrow(pb), depth_cut))
deep <- rbindlist(lapply(spec, function(s)
    rbindlist(lapply(OUTCOMES, function(o) {
        r <- fit_one(pb_deep, o, s$preds)
        if (is.null(r)) return(NULL)
        r[, model := s$tag][]
    }), fill = TRUE)), fill = TRUE)
if (nrow(deep)) {
    deep[, `:=`(arm_set = "deep_half", n_units = nrow(pb_deep),
                depth_cut = depth_cut)]
    fwrite(deep, file.path(outdir, "tgfb_specificity_depth_sensitivity.tsv"),
           sep = "\t")
}

## ---- Test C: where does each arm's variance live? ---------------------------
## If the IEG arm is a warm-dissociation artifact it should partition markedly
## more variance to `study` (i.e. to protocol) than the SMAD arm does.
##
## THIS MUST RUN ON THE RAW SCORES, NOT THE _z ONES. `dataset` nests strictly
## inside `study` here (33 dataset/study pairs, one study per dataset), so
## z_within_dataset() centres away the between-study variance by construction
## and the z version of this test reports ~0% for every score no matter what is
## true. Both are emitted, `scale` distinguishes them, and only the raw rows
## answer the question -- the z rows are kept solely to show the test is vacuous
## on that scale, so nobody re-runs it there and reads the zeros as a result.
vc_scores <- c("tgfb_response_score", "tgfb_smad_score", "tgfb_ieg_score")
vc <- rbindlist(lapply(c("raw", "z"), function(scale) {
    rbindlist(lapply(vc_scores, function(v) {
        vv_col <- if (scale == "z") paste0(v, "_z") else v
        if (!vv_col %in% names(pb)) return(NULL)
        f <- reformulate(c("1", "(1 | study)", "(1 | dataset)", "(1 | donor_id)"),
                         response = vv_col)
        fit <- try(suppressMessages(lmer(f, data = pb)), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        vv <- as.data.table(VarCorr(fit))
        tot <- sum(vv$vcov)
        data.table(scale = scale, score = v, component = vv$grp,
                   variance = vv$vcov, pct = 100 * vv$vcov / tot)
    }), fill = TRUE)
}), fill = TRUE)
fwrite(vc, file.path(outdir, "tgfb_specificity_varcomp.tsv"), sep = "\t")
if (nrow(vc)) {
    message("\n---- Test C: between-study variance, RAW scale ----")
    print(vc[scale == "raw" & component == "study", .(score, pct)])
}

## ---- detection-matched empirical null ---------------------------------------
if (!is.na(opt$null_pb) && file.exists(opt$null_pb)) {
    nl <- fread(opt$null_pb)
    null_cols <- grep("^null_", names(nl), value = TRUE)
    ## Join keys must agree in type. `pb$pericyte_state` is a factor (set on the
    ## cell table above) while the null table round-trips through TSV and comes
    ## back integer, so data.table refuses the merge outright. Cast both to
    ## character rather than relying on either side's storage type.
    nl[, pericyte_state := as.character(pericyte_state)]
    nl[, donor_id := as.character(donor_id)]
    pb_keys <- pb[, .(donor_id = as.character(donor_id),
                      pericyte_state = as.character(pericyte_state),
                      study, dataset, mean_log10_counts,
                      basement_membrane_score_z, bm_minus_fibrillar_z)]
    pbn <- merge(pb_keys, nl, by = c("donor_id", "pericyte_state"))
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

    ## ---- leave-one-gene-out on the full 17-gene panel -----------------------
    ## Sourced from the same file as the null panels. An earlier version looked
    ## for these columns in the cell table, where they never exist, so the whole
    ## sensitivity skipped without a word -- the exact silent-skip this plan was
    ## written to prevent. Absence is now an error, not a shrug.
    logo_cols <- grep("^logo_", names(nl), value = TRUE)
    if (!length(logo_cols)) {
        warning("null pseudobulk carries no logo_* columns; the leave-one-gene-out ",
                "sensitivity did NOT run. Re-run 11.tgfb_null_panels.py -- it "
                , "writes them alongside the null panels.", call. = FALSE)
    } else {
        for (cl in logo_cols)
            pbn[[paste0(cl, "_z")]] <- z_within_dataset(pbn[[cl]], pbn$dataset)
        logo <- rbindlist(lapply(logo_cols, function(cl) {
            rbindlist(lapply(OUTCOMES, function(o) {
                r <- fit_one(pbn, o, paste0(cl, "_z"))
                if (is.null(r)) return(NULL)
                r[, dropped_gene := sub("^logo_", "", cl)][]
            }), fill = TRUE)
        }), fill = TRUE)
        if (nrow(logo)) {
            setorder(logo, outcome, estimate)
            fwrite(logo, file.path(outdir, "tgfb_specificity_logo.tsv"), sep = "\t")
            lb <- logo[outcome == "basement_membrane_score_z"]
            message(sprintf(
                "Leave-one-out (BM): beta %.4f to %.4f across %d drops; full panel %.4f",
                min(lb$estimate), max(lb$estimate), nrow(lb),
                obs[model == "full_panel" &
                    outcome == "basement_membrane_score_z"]$estimate[1]))
        }
    }

    ## Empirical p and power, per arm and outcome.
    obs_alone <- obs[model %in% c("smad_alone", "ieg_alone")]
    obs_alone[, arm := fifelse(model == "smad_alone", "smad", "ieg")]

    summ <- rbindlist(lapply(seq_len(nrow(obs_alone)), function(i) {
        row <- obs_alone[i]
        nd <- nullres[arm == row$arm & outcome == row$outcome & converged == TRUE]
        if (!nrow(nd)) return(NULL)
        ## THE NULL IS NOT CENTRED ON ZERO, and everything here depends on that.
        ## A random panel of well-detected genes predicts the BM score at about
        ## +0.45, because any sc.tl.score_genes score shares a "general
        ## expression level" component with the BM score that the depth
        ## covariate does not fully absorb. So an empirical p built from
        ## |null| >= |observed| is invalid: it asks whether random panels have
        ## large effects (they do) instead of whether THIS panel is unusual.
        ## Both the p-value and the power statement are therefore referred to
        ## the null's OWN centre.
        mu <- mean(nd$estimate); sdev <- sd(nd$estimate)
        emp_p <- (1 + sum(abs(nd$estimate - mu) >= abs(row$estimate - mu))) /
                 (1 + nrow(nd))
        p_lower <- (1 + sum(nd$estimate <= row$estimate)) / (1 + nrow(nd))
        ## Power: the deviation from the null centre that 80% of matched panels
        ## fall within -- i.e. the smallest shift this design can resolve. If it
        ## is below ref_beta, an arm sitting inside its null is a real null and
        ## not a power failure.
        detectable <- quantile(abs(nd$estimate - mu), 0.80)
        data.table(
            arm = row$arm, outcome = row$outcome,
            beta_observed = row$estimate, p_model = row$p_value,
            n_null = nrow(nd),
            null_mean = mu, null_sd = sdev,
            null_q025 = quantile(nd$estimate, 0.025),
            null_q975 = quantile(nd$estimate, 0.975),
            z_vs_null = (row$estimate - mu) / sdev,
            empirical_p = emp_p, empirical_p_lower = p_lower,
            detectable_effect_80 = detectable, ref_beta = opt$ref_beta)
    }), fill = TRUE)

    ## The verdict table from TGFB_SPECIFICITY_PLAN.md section 7, applied
    ## mechanically so the rule cannot drift after seeing the numbers.
    summ[, outside_null := empirical_p < 0.05]
    ## "Well powered" means the design can resolve a shift of ref_beta from the
    ## null centre, so an arm inside its null is informative rather than mute.
    summ[, well_powered := detectable_effect_80 <= ref_beta]
    summ[, arm_verdict := fifelse(
        outside_null & beta_observed < null_mean, "carries the association",
        fifelse(!outside_null & well_powered, "null, adequately powered",
                "inside null, underpowered -- uninformative"))]
    fwrite(summ, file.path(outdir, "tgfb_specificity_summary.tsv"), sep = "\t")

    message("\n---- verdict inputs ----")
    print(summ[, .(arm, outcome, beta_observed, null_mean, z_vs_null,
                   empirical_p, detectable_effect_80, arm_verdict)])
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
