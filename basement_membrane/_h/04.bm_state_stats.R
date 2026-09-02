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
##
## Also tests AGTR1 against BM directly (--denoise), at both the pseudobulk and
## the within-donor level. Figure panel F pairs a positive AGTR1-vs-pseudotime
## correlation with a negative BM-vs-pseudotime one; that is two correlations
## through a shared mediator, not an AGTR1-BM antagonism, and the manuscript may
## not claim the latter without this test. Outputs: bm_vs_agtr1_corr.tsv,
## bm_vs_agtr1_models.tsv, bm_vs_agtr1_within_donor.tsv.
##
## 2026-09-01 (collaborator request). Three changes, all reusing the machinery
## that block already had:
##   (1) The outcome set is no longer BM alone. Both matrix categories are tested
##       -- basement membrane, fibrillar collagen, and their DIFFERENCE -- because
##       the question is which of the two matrices an exposure is associated with,
##       and a per-category association cannot answer that on its own. The
##       difference is the estimand that does: a multiplicative capture constant
##       applies to both scores and cancels in it.
##   (2) TGF-beta signalling enters as a second predictor, run through the exact
##       same two-level design as AGTR1 so the two are directly comparable. The
##       response panel (bm_panels.TGFB_RESPONSE) shares no gene with either
##       matrix panel or with any pericyte_states program -- asserted in
##       bm_panels.py -- so the association is not arithmetic.
##       Headline row: outcome `bm_minus_fibrillar_z`, predictor
##       `tgfb_response_score`. Positive = TGF-beta exposure goes with a BM-shifted
##       matrix; negative = with a fibril-shifted one.
##   (3) `bm_v1_score` (the frozen 13-gene panel) is scored alongside the current
##       20-gene panel and compared in bm_v1_vs_v2_concordance.tsv, so the effect
##       of the 2026-09-01 panel expansion on every cluster-level estimate is
##       auditable. It is provenance, never a result.

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
    ## scVI-denoised AGTR1, same table and model that
    ## pericyte_states/_h/03.agtr1_lenses.R uses as its robust lens.
    make_option("--denoise", type = "character", default = NA_character_),
    make_option("--den-model", type = "character", default = "Pericyte-only-trained",
                dest = "den_model"),
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
score_cols <- c("basement_membrane_score", "bm_v1_score",
                "bm_collagen_iv_score",
                "bm_laminin_score", "bm_linker_score", "fibrillar_ecm_score",
                "fibrillar_collagen_score",
                "fibrillar_core_score", "fibrillar_minor_score",
                "ambient_tracer_score", "tgfb_response_score",
                "tgfb_response_noECM_score", "tgfb_receptor_score",
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

## --------------------- direct matrix-vs-predictor associations (AGTR1, TGF-b) ----
## Panel F of figure_pericyte_layer shows AGTR1 rising and BM falling along
## pseudotime. Those are two SEPARATE correlations sharing a mediator: they do not
## establish that AGTR1 and BM oppose each other, and until now nothing in this
## module or in pericyte_states/ tested the pair directly. This block does, for
## both matrix categories and for two predictors.
##
## The AGTR1 readout is the scVI-DENOISED lens, not the raw one. Raw AGTR1 and the
## matrix panels both scale with capture depth, so a raw-lens association here is
## confounded by precisely the dropout that 03.agtr1_lenses.R showed manufactures
## AGTR1's apparent program bias -- reporting a raw-only result would repeat the
## error the three-lens analysis exists to catch. Raw and detection lenses are
## carried alongside as the dropout-sensitive comparison: raw-negative with
## denoised-null is the signature of shared dropout, not of biology (same logic as
## the ACTA2 control in pericyte_states/_h/05.acta2_control.R).
##
## Two levels are tested because they answer different questions, and panel F only
## answers the first: donor x cluster pseudobulk (do lungs/clusters with more of
## the exposure carry more of this matrix?) and within-donor cell-level (does the
## exposure track matrix INSIDE a lung?). The BM-vs-pseudotime trend is already
## known to be between-donor only -- flat within donors, bm_continuum_summary.tsv
## median rho 0.02, p = 0.84 -- so the within-donor test is the one that would
## license a mechanistic claim in the text.
##
## OUTCOMES. Testing BM alone cannot say which matrix an exposure goes with, so
## the outcome set is both categories plus their difference. Only the difference
## is invariant to a per-unit multiplicative capture constant, which makes it the
## estimand for "shifts the matrix toward BM" as opposed to "more matrix overall".
## `bm_resid` (BM orthogonalized against the frozen fibrillar_ecm panel) is kept
## as the fourth outcome so a surviving slope there means the exposure tracks BM
## specifically rather than the shared stromal-ECM axis.

matrix_readme <- character(0)

## The difference is taken on the RAW panel scores and standardized afterwards,
## not as a difference of two separately standardized scores. Differencing two
## within-dataset z-scores rescales each side by its own SD, which destroys the
## very cancellation the difference is for: the two panels do not have equal
## dispersion, so the "shared capture constant cancels" argument only holds on the
## raw scale. Depth still enters every model as a covariate -- the cancellation is
## an additional defense, not the only one.
if (all(c("basement_membrane_score", "fibrillar_collagen_score") %in% names(pb))) {
    pb[, bm_minus_fibrillar := basement_membrane_score - fibrillar_collagen_score]
    pb[, bm_minus_fibrillar_z := z_within_dataset(bm_minus_fibrillar, dataset)]
}
## Second orthogonalization, matching the category the collaborator asked about.
## `bm_resid` residualizes BM on the FROZEN mixed fibrillar_ecm panel (collagens +
## proteoglycans + glycoproteins), which answers "is this cell fibroblast-like";
## `bm_resid_collagen` residualizes on the fibrillar COLLAGEN panel, which answers
## "does this cell build interstitial fibrils". Reporting only the first would
## orthogonalize against a different contrast than the one being plotted.
if ("fibrillar_collagen_score_z" %in% names(pb))
    pb[, bm_resid_collagen := residuals(suppressMessages(lmer(
        basement_membrane_score_z ~ fibrillar_collagen_score_z +
            mean_log10_counts + (1 | study) + (1 | donor_id), data = pb)))]
OUTCOMES <- intersect(c("basement_membrane_score_z", "fibrillar_collagen_score_z",
                        "bm_minus_fibrillar_z", "bm_resid", "bm_resid_collagen"),
                      names(pb))
message("Matrix outcomes: ", paste(OUTCOMES, collapse = ", "))

dl <- copy(d)
if (all(c("basement_membrane_score", "fibrillar_collagen_score") %in% names(dl)))
    dl[, bm_minus_fibrillar := basement_membrane_score - fibrillar_collagen_score]
CELL_OUTCOMES <- intersect(c("basement_membrane_score", "fibrillar_collagen_score",
                             "bm_minus_fibrillar"), names(dl))

## ---- predictors -------------------------------------------------------------
agtr1_cols <- intersect(c("AGTR1_expr", "AGTR1_detect"), names(dl))
if (!is.na(opt$denoise) && file.exists(opt$denoise)) {
    den <- fread(opt$denoise, select = c("index", "Model", "AGTR1_scvi"))
    den <- den[Model == opt$den_model]
    if (anyDuplicated(den$index))
        stop("denoised AGTR1 table has duplicate cells for model ", opt$den_model)
    dl <- merge(dl, den[, .(index, AGTR1_scvi)], by = "index", all.x = TRUE)
    message(sprintf("Denoised AGTR1 (%s): %d/%d cells (%.1f%%)", opt$den_model,
                    sum(is.finite(dl$AGTR1_scvi)), nrow(dl),
                    100 * sum(is.finite(dl$AGTR1_scvi)) / nrow(dl)))
    agtr1_cols <- c("AGTR1_scvi", agtr1_cols)
} else {
    warning("no --denoise table supplied; the denoised lens is a concordant ",
            "sensitivity for AGTR1-vs-matrix, so the raw-lens results below are ",
            "unpaired. The arbiter is 10.agtr1_count_models.R either way.",
            call. = FALSE)
}
tgfb_cols <- intersect(c("tgfb_response_score", "tgfb_response_noECM_score",
                         "tgfb_receptor_score"), names(dl))

## ---- one generic association block, run once per predictor family -----------
## Identical design for AGTR1 and TGF-beta so the two are directly comparable:
## same units, same covariate, same random effects, same within-donor estimator.
assoc_block <- function(pred_cols, tag) {
    pred_cols <- intersect(pred_cols, names(dl))
    if (!length(pred_cols) || !length(OUTCOMES)) {
        message("No predictors/outcomes available for '", tag, "' -- skipping.")
        return(NULL)
    }
    dl[, (pred_cols) := lapply(.SD, as.numeric), .SDcols = pred_cols]

    ## (i) donor x cluster pseudobulk, on the SAME units as every other test here.
    ##
    ## Only predictors NOT already aggregated into `pb` are merged in. The TGF-beta
    ## panels are themselves panel scores and are therefore already columns of
    ## `pb`; merging a second copy suffixed both to .x/.y, so every `pba[[pred]]`
    ## lookup silently returned NULL and the whole TGF-beta pseudobulk arm wrote
    ## an empty table while the within-donor arm looked fine. Assert instead.
    new_preds <- setdiff(pred_cols, names(pb))
    pba <- if (length(new_preds)) {
        apb <- dl[, lapply(.SD, mean, na.rm = TRUE),
                  by = .(donor_id, pericyte_state), .SDcols = new_preds]
        merge(pb, apb, by = c("donor_id", "pericyte_state"))
    } else copy(pb)
    missing_pred <- setdiff(pred_cols, names(pba))
    if (length(missing_pred))
        stop("[", tag, "] predictors missing from the pseudobulk table after the ",
             "merge: ", paste(missing_pred, collapse = ", "),
             " -- refusing to write an empty result table.")
    for (cl in pred_cols)
        pba[[paste0(cl, "_z")]] <- z_within_dataset(pba[[cl]], pba$dataset)

    ## Unadjusted correlations against the raw matrix scores, formatted to sit
    ## beside bm_vs_fibrillar_corr.tsv.
    raw_matrix <- intersect(c("basement_membrane_score", "fibrillar_collagen_score"),
                            names(pba))
    corr <- rbindlist(lapply(pred_cols, function(pc) rbindlist(lapply(
        raw_matrix, function(mc) {
            ok <- is.finite(pba[[pc]]) & is.finite(pba[[mc]])
            if (sum(ok) < 10) return(NULL)
            x <- pba[[mc]][ok]; y <- pba[[pc]][ok]
            ct <- cor.test(x, y, method = "pearson")
            cs <- suppressWarnings(cor.test(x, y, method = "spearman"))
            data.table(comparison = paste0(sub("_score$", "", mc), " vs ", pc),
                       matrix_score = mc, predictor = pc,
                       pearson_r = unname(ct$estimate),
                       pearson_ci_lo = ct$conf.int[1],
                       pearson_ci_hi = ct$conf.int[2],
                       pearson_p = ct$p.value,
                       spearman_rho = unname(cs$estimate),
                       spearman_p = cs$p.value, n_units = sum(ok))
        }), fill = TRUE)), fill = TRUE)
    write_tsv_safe(corr, file.path(opt$outdir, sprintf("bm_vs_%s_corr.tsv", tag)))

    ## Depth-adjusted mixed models. `bm_resid` is already residualized on depth,
    ## so it does not take the covariate a second time.
    fit_one <- function(pc, outcome) {
        pz <- paste0(pc, "_z")
        if (!pz %in% names(pba) || !outcome %in% names(pba)) return(NULL)
        ## The residual outcomes are already residualized on depth; adding the
        ## covariate again would double-adjust.
        extra <- if (grepl("^bm_resid", outcome)) character() else
            "mean_log10_counts"
        f <- tryCatch(suppressMessages(lmer(
            reformulate(c(pz, extra, "(1 | study)", "(1 | donor_id)"), outcome),
            data = pba)), error = function(e) NULL)
        if (is.null(f)) return(NULL)
        cf <- summary(f)$coefficients
        if (!pz %in% rownames(cf)) return(NULL)
        ## Satterthwaite `df`/`Pr(>|t|)` come from lmerTest::lmer. Degrade to NA
        ## rather than erroring if a bare lme4 fit ever reaches here.
        col <- function(nm) if (nm %in% colnames(cf)) cf[pz, nm] else NA_real_
        data.table(outcome = outcome, predictor = pc, family = tag,
                   estimate = col("Estimate"), SE = col("Std. Error"),
                   df = col("df"), t_ratio = col("t value"),
                   p_value = col("Pr(>|t|)"),
                   depth_adjusted = "mean_log10_counts" %in% extra,
                   n_units = nrow(pba), n_donors = uniqueN(pba$donor_id))
    }
    mods <- rbindlist(unlist(lapply(OUTCOMES, function(o)
        lapply(pred_cols, fit_one, outcome = o)), recursive = FALSE), fill = TRUE)
    if (!nrow(mods))
        stop("[", tag, "] every model failed although ", length(pred_cols),
             " predictors and ", length(OUTCOMES), " outcomes were available; ",
             "refusing to write an empty result table.")
    if (nrow(mods)) {
        ## BH within outcome: the predictors are the family being screened, and
        ## correcting across outcomes would penalise the difference endpoint for
        ## its own two components.
        mods[, p_BH := p.adjust(p_value, method = "BH"), by = outcome]
        write_tsv_safe(mods, file.path(opt$outdir,
                                       sprintf("bm_vs_%s_models.tsv", tag)))
    }

    ## (ii) within-donor, cell level: per-donor Spearman then a one-sample test on
    ## the donor rhos. Never pools cells across donors -- the same pattern as the
    ## continuum block below, so the two are directly comparable.
    ##
    ## The PARTIAL rho, controlling per-cell sequencing depth, is reported beside
    ## the raw one and is the readout. Every sc.tl.score_genes score rises with
    ## capture depth, so two scores correlated within a donor is the expected
    ## result under pure technical variation -- this is the same confound that
    ## produced the retracted AGTR1-marks-the-stabilizing-pole claim, and the
    ## pseudobulk arm defends against it with a depth covariate while an
    ## unadjusted cell-level rho would not. Same rank-residual estimator as
    ## pericyte_states/_h/02.continuum_dpt.py:partial_spearman.
    pspear <- function(x, y, z) {
        ok <- is.finite(x) & is.finite(y) & is.finite(z)
        if (sum(ok) < 5) return(NA_real_)
        rx <- rank(x[ok]); ry <- rank(y[ok]); rz <- rank(z[ok])
        if (stats::sd(rz) == 0) return(suppressWarnings(cor(rx, ry)))
        suppressWarnings(cor(residuals(stats::lm(rx ~ rz)),
                             residuals(stats::lm(ry ~ rz))))
    }
    depth_col <- if ("log10_total_counts" %in% names(dl)) "log10_total_counts" else NA
    wd <- rbindlist(unlist(lapply(CELL_OUTCOMES, function(oc)
        lapply(pred_cols, function(pc) {
            sub <- dl[is.finite(get(pc)) & is.finite(get(oc))]
            if (!nrow(sub)) return(NULL)
            dr <- sub[, .(n = .N,
                          rho = suppressWarnings(
                              cor(get(pc), get(oc), method = "spearman")),
                          partial_rho = if (is.na(depth_col)) NA_real_ else
                              pspear(get(pc), get(oc), get(depth_col))),
                      by = donor_id][n >= 20 & is.finite(rho)]
            if (nrow(dr) < 3) return(NULL)
            w <- suppressWarnings(wilcox.test(dr$rho))
            pv <- dr$partial_rho[is.finite(dr$partial_rho)]
            wp <- if (length(pv) >= 3) suppressWarnings(wilcox.test(pv)) else NULL
            data.table(outcome = oc, predictor = pc, family = tag,
                       n_donors = nrow(dr), median_rho = median(dr$rho),
                       q25 = unname(quantile(dr$rho, 0.25)),
                       q75 = unname(quantile(dr$rho, 0.75)),
                       p_wilcox = w$p.value,
                       n_donors_partial = length(pv),
                       median_partial_rho = if (length(pv)) median(pv) else NA_real_,
                       p_wilcox_partial = if (is.null(wp)) NA_real_ else wp$p.value,
                       depth_control = if (is.na(depth_col)) "NONE" else depth_col)
        })), recursive = FALSE), fill = TRUE)
    if (nrow(wd)) {
        wd[, p_BH := p.adjust(p_wilcox, method = "BH"), by = outcome]
        wd[, p_BH_partial := p.adjust(p_wilcox_partial, method = "BH"), by = outcome]
        write_tsv_safe(wd, file.path(opt$outdir,
                                     sprintf("bm_vs_%s_within_donor.tsv", tag)))
    }

    message(sprintf("[%s] pseudobulk correlations:", tag));  print(corr)
    message(sprintf("[%s] depth-adjusted mixed models:", tag)); print(mods)
    message(sprintf("[%s] within-donor Spearman:", tag));    print(wd)
    list(corr = corr, models = mods, within = wd)
}

res_agtr1 <- assoc_block(agtr1_cols, "agtr1")
res_tgfb  <- assoc_block(tgfb_cols, "tgfb")

## The one row the collaborator's hypothesis turns on: does TGF-beta exposure go
## with a BM-shifted or a fibril-shifted matrix? Surfaced explicitly so it cannot
## be lost among the other 11 rows of the model table.
if (!is.null(res_tgfb) && nrow(res_tgfb$models)) {
    hd <- res_tgfb$models[outcome == "bm_minus_fibrillar_z" &
                          predictor == "tgfb_response_score"]
    if (nrow(hd))
        message(sprintf(
            "TGF-beta directional headline: BM - fibrillar ~ TGF-beta response = %.4f (SE %.4f, p = %.3g, BH = %.3g) -- %s",
            hd$estimate, hd$SE, hd$p_value, hd$p_BH,
            if (!is.finite(hd$p_BH) || hd$p_BH >= 0.05) "NOT significant" else
                if (hd$estimate > 0) "BM-shifted" else "fibril-shifted"))
    lo <- res_tgfb$models[outcome == "bm_minus_fibrillar_z" &
                          predictor == "tgfb_response_noECM_score"]
    if (nrow(hd) && nrow(lo))
        message(sprintf(
            "  TGFBI leave-out sensitivity: %.4f (p = %.3g) vs %.4f with TGFBI -- %s",
            lo$estimate, lo$p_value, hd$estimate,
            if (is.finite(lo$estimate) && is.finite(hd$estimate) &&
                sign(lo$estimate) == sign(hd$estimate)) "same direction" else
                "DIRECTION CHANGES; the association is carried by an ECM gene"))
}

matrix_readme <- c(
    "",
    "Matrix-vs-predictor associations. Outcomes are BOTH matrix categories, their",
    "difference (taken on the raw scores, then standardized) and two",
    "orthogonalizations. Depth is a covariate in every pseudobulk model and is",
    "partialled out of every within-donor rho: all score_genes scores rise with",
    "capture depth, so an unadjusted correlation between two of them is the",
    "expected result under pure technical variation.",
    "",
    "WHICH LENS ARBITRATES (changed 2026-09-02). The denoised lens was originally",
    "designated the AGTR1 readout. It is NOT the arbiter any more: it is a model",
    "output that can in principle manufacture or erase an association by borrowing",
    "across genes, and the first version of it was trained on the soupX float",
    "matrix and failed its validity gate (donor rho 0.014, p 0.91). It is now",
    "reported as a CONCORDANT SENSITIVITY. Between-group contrasts are arbitrated",
    "by 10.agtr1_count_models.R -- AGTR1 integer counts from raw/X as the response",
    "of an NB GLMM with offset(log(library size)), so dropout sits inside the",
    "likelihood and nothing is imputed. raw/detection remain the dropout-sensitive",
    "comparison. Note the matrix scores themselves have no denoised counterpart, so",
    "EVERY row here still has one depth-affected side, which is the other reason",
    "the offset-bearing count model carries the claim.",
    paste(utils::capture.output(print(
        rbindlist(Filter(Negate(is.null),
                         list(res_agtr1$models, res_tgfb$models)), fill = TRUE))),
          collapse = "\n"),
    "",
    "Within-donor (cell-level Spearman per donor, >=20 cells, one-sample test):",
    paste(utils::capture.output(print(
        rbindlist(Filter(Negate(is.null),
                         list(res_agtr1$within, res_tgfb$within)), fill = TRUE))),
          collapse = "\n"))

## --------------------------- AUDIT: 13-gene vs 20-gene BM panel ----
## The 2026-09-01 expansion moves basement_membrane_score, and therefore every
## BM estimate in this file. This table makes that movement measurable instead of
## inferred: score correlation at both levels, and the by-cluster marginal means
## side by side. `bm_v1_score` is provenance and must never carry a claim.
if (all(c("bm_v1_score", "basement_membrane_score") %in% names(pb))) {
    conc_row <- function(x, y, level, n)
        data.table(level = level,
                   pearson_r = suppressWarnings(cor(x, y)),
                   spearman_rho = suppressWarnings(cor(x, y, method = "spearman")),
                   n = n)
    conc <- rbind(
        conc_row(d$bm_v1_score, d$basement_membrane_score, "cell", nrow(d)),
        conc_row(pb$bm_v1_score, pb$basement_membrane_score, "donor_x_cluster", nrow(pb)))
    emm_both <- rbindlist(Filter(Negate(is.null), lapply(
        c("basement_membrane_score_z", "bm_v1_score_z"), function(cl) {
            r <- by_cluster(cl); if (is.null(r)) NULL else as.data.table(r$emm)
        })), fill = TRUE)
    wide <- dcast(emm_both, pericyte_state ~ score, value.var = "emmean")
    if (all(c("basement_membrane_score_z", "bm_v1_score_z") %in% names(wide)))
        wide[, delta := basement_membrane_score_z - bm_v1_score_z]
    write_tsv_safe(conc, file.path(opt$outdir, "bm_v1_vs_v2_concordance.tsv"))
    write_tsv_safe(wide, file.path(opt$outdir, "bm_v1_vs_v2_by_cluster.tsv"))
    message("BM panel version concordance (13-gene vs 20-gene):"); print(conc)
    print(wide)
}

## ------------------------------------------- BM along the injury continuum ----
if (!is.na(opt$continuum) && file.exists(opt$continuum)) {
    cont <- fread(opt$continuum)
    setnames(cont, 1, "index")
    pt_col <- intersect(c("dpt_pseudotime", "pseudotime"), names(cont))
    if (length(pt_col)) {
        ## Score columns carried onto the continuum. The first three metrics
        ## below reproduce the published rho_bm / rho_fib / rho_switch exactly;
        ## the 2026-08-10 additions resolve the fibrillar side into the
        ## fibril-forming (I/III) and regulatory (V/XI) blocks and add the
        ## ambient-tracer trend as a control -- if soup burden itself rises along
        ## pseudotime, the fibrillar rise is confounded and cannot be claimed.
        want <- c("basement_membrane_score", "fibrillar_ecm_score",
                  "fibrillar_collagen_score",
                  "fibrillar_core_score", "fibrillar_minor_score",
                  "ambient_tracer_score", "tgfb_response_score")
        want <- intersect(want, names(d))
        dc <- merge(d[, .SD, .SDcols = c("index", "donor_id",
                                         "log10_total_counts", want)],
                    cont[, .SD, .SDcols = c("index", pt_col[1])], by = "index")
        setnames(dc, pt_col[1], "pt")

        ## metric -> per-cell vector, so every rho is computed the same way.
        metrics <- list(rho_bm = "basement_membrane_score",
                        rho_fib = "fibrillar_ecm_score",
                        rho_collagen = "fibrillar_collagen_score",
                        rho_core = "fibrillar_core_score",
                        rho_minor = "fibrillar_minor_score",
                        rho_tracer = "ambient_tracer_score",
                        rho_tgfb = "tgfb_response_score")
        metrics <- metrics[vapply(metrics, `%in%`, logical(1), want)]
        switches <- list(
            rho_switch = c("basement_membrane_score", "fibrillar_ecm_score"),
            rho_switch_core = c("basement_membrane_score", "fibrillar_core_score"),
            rho_switch_collagen = c("basement_membrane_score",
                                    "fibrillar_collagen_score"))
        switches <- switches[vapply(switches, function(p) all(p %in% want),
                                    logical(1))]

        ## Per-donor Spearman, then a one-sample test on the donor rhos -- never
        ## pools cells across donors (the existing pseudotime-trend pattern).
        spear <- function(x, pt) suppressWarnings(cor(x, pt, method = "spearman"))
        donor_rho <- dc[, c(
            .(n = .N),
            lapply(metrics, function(cl) spear(get(cl), pt)),
            lapply(switches, function(p) spear(get(p[1]) - get(p[2]), pt))),
            by = donor_id][n >= 20]
        write_tsv_safe(donor_rho, file.path(opt$outdir, "bm_continuum_donor_rho.tsv"))

        summ <- rbindlist(lapply(c(names(metrics), names(switches)), function(cc) {
            v <- donor_rho[[cc]]; v <- v[is.finite(v)]
            if (length(v) < 3) return(NULL)
            w <- suppressWarnings(wilcox.test(v))
            data.table(metric = cc, n_donors = length(v), median_rho = median(v),
                       q25 = unname(quantile(v, 0.25)),
                       q75 = unname(quantile(v, 0.75)),
                       p_wilcox = w$p.value)
        }), fill = TRUE)
        ## BH over the three original metrics only, so the published adjusted
        ## p-values do not move because the panel grew; the expansion metrics
        ## are corrected as their own family.
        orig <- c("rho_bm", "rho_fib", "rho_switch")
        summ[, bh_family := fifelse(metric %in% orig, "original", "expansion")]
        summ[, p_BH := p.adjust(p_wilcox, method = "BH"), by = bh_family]
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
    paste(utils::capture.output(print(corr_tbl)), collapse = "\n"),
    matrix_readme)
writeLines(readme, file.path(opt$outdir, "bm_pericyte_axis_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
