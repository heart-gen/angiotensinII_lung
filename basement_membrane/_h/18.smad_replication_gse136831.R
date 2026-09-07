## =============================================================================
## Independent evaluation of the SMAD -> (BM - fibrillar) association in
## GSE136831 (Adams/Kaminski).  Defect P1-20(a).
##
## WHAT IS BEING EVALUATED.  `04.bm_state_stats.R` + the pre-registered arm
## decomposition (`TGFB_SPECIFICITY_PLAN.md`) report, in HLCA pericytes:
##
##     smad arm -> bm_minus_fibrillar_z :  beta = -0.187,
##                 z vs matched null = -2.18, empirical p = 0.013
##
## and, in the same run, that the SMAD arm is an adequately powered NULL on the
## BM score alone (+0.026, p = 0.70) while the IEG arm carries that endpoint.
## The BM - fibrillar result sits on the plan's SECONDARY endpoint and was
## GENERATED, not confirmed, by the analysis that reports it.  GSE136831 is the
## only independent dataset in this repository with a pericyte compartment
## defined, so it is where the claim can be moved from "generated here" to "held
## up once".
##
## WHAT THIS CAN AND CANNOT DO, stated before any number.
##
##  * THE DONOR SET IS THIN AND DISEASE-SKEWED.  At a 5-cell floor GSE136831 has
##    22 pericyte units: 1 Control, 6 COPD, 15 IPF.  That is fine for an
##    ASSOCIATION between two panel scores measured on the same cells, which is
##    what the claim is, and it is NOT enough for anything disease-stratified.
##    No disease contrast is reported here; `disease` enters only as a nuisance
##    term.  See `06.bm_copd_stats.R` for why a pericyte disease test is not
##    estimable in this dataset at all.
##
##  * ONLY THE BM - FIBRILLAR ENDPOINT IS EVALUABLE, and the reason is the
##    non-zero null established in [gene-set-score-null-not-zero].  A score-on-
##    score coefficient does not have 0 as its null: in HLCA the matched null for
##    smad -> BM alone is +0.047 (sd 0.070), so beta-vs-0 is the WRONG test
##    there.  For smad -> bm_minus_fibrillar the matched null mean is +0.026
##    (sd 0.098) -- a contrast of two panels differences most of the shared
##    non-zero component away -- so beta-vs-0 is approximately the right test on
##    THAT endpoint and only that one.
##
##    The BM-alone endpoint would need a matched null of random gene panels, and
##    this pseudobulk cannot supply one: `gse136831_bm_pseudobulk.tsv.gz` was
##    built with `--genes bm_panel_genes.tsv`, so it contains the panel genes and
##    nothing to draw a background from.  Building one means another pass over
##    the 8.4 GB h5ad with a large random gene list.  The BM-alone arm is
##    therefore reported as NOT EVALUABLE rather than tested against 0.
##
##  * THE IEG ARM IS THE DISCRIMINATING CONTROL, not a second test.  In HLCA it
##    carries 59.2% between-study variance against the SMAD arm's 0.0%, which is
##    a warm-dissociation signature.  GSE136831 is a SINGLE study processed one
##    way, so that particular confound cannot express itself here -- which makes
##    this dataset informative about the SMAD arm specifically.  If IEG carried
##    BM - fibrillar here too, the specificity claim would weaken.
##
##  * SIGNS AND SIGNIFICANCE, NOT MAGNITUDES.  HLCA scores are z-scored within
##    dataset; these are z-scored within compartment inside one study.  The two
##    betas are on different scales and are never differenced or averaged.
## =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--panels", type = "character"),
    make_option("--outdir", type = "character", default = "stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--compartment", type = "character", default = "Pericyte")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) fwrite(as.data.table(x), file.path(opt$outdir, f), sep = "\t")

pan <- fread(opt$panels)
pb  <- fread(opt$pseudobulk)

## Panel scores as the mean over the panel genes PRESENT in the pseudobulk. The
## present/total count travels with every row: a panel that lost half its genes
## to the gene filter is a different panel, and that must be visible.
panel_score <- function(dt, genes) {
    cols <- paste0(genes, "__expr")
    cols <- intersect(cols, names(dt))
    list(score = if (length(cols)) rowMeans(as.matrix(dt[, ..cols]), na.rm = TRUE)
                 else rep(NA_real_, nrow(dt)),
         n_present = length(cols), n_total = length(genes))
}

PANELS <- c("basement_membrane", "fibrillar_collagen", "tgfb_smad", "tgfb_ieg")
cov <- rbindlist(lapply(PANELS, function(p) {
    g <- pan[panel == p, unique(gene)]
    r <- panel_score(pb, g)
    pb[[paste0(p, "_score")]] <<- r$score
    data.table(panel = p, n_genes_panel = r$n_total, n_genes_present = r$n_present,
               genes_missing = paste(setdiff(g, sub("__expr$", "", intersect(paste0(g, "__expr"), names(pb)))),
                                     collapse = ";"))
}))
wt(cov, "smad_replication_panel_coverage.tsv")
cat("== panel coverage in the GSE136831 pseudobulk ==\n"); print(cov)

d <- pb[compartment == opt$compartment & n_cells >= opt$min_cells]
cat(sprintf("\n%s units at >=%d cells: %d\n", opt$compartment, opt$min_cells, nrow(d)))
print(table(d$disease))

## Guard rather than assume. Two levels of `disease` are needed for the nuisance
## term to be estimable, and a singleton arm contributes no residual information.
keep_dx <- names(which(table(d$disease) >= 2L))
if (length(keep_dx) < length(unique(d$disease)))
    cat("dropping disease arms with <2 units:",
        setdiff(unique(d$disease), keep_dx), "\n")
d <- d[disease %in% keep_dx]
d[, disease := factor(disease)]

zc <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
d[, `:=`(bm_z        = zc(basement_membrane_score),
         fib_z       = zc(fibrillar_collagen_score),
         smad_z      = zc(tgfb_smad_score),
         ieg_z       = zc(tgfb_ieg_score))]
d[, bm_minus_fibrillar_z := zc(basement_membrane_score - fibrillar_collagen_score)]

## HLCA reference values, hard-coded from tgfb_specificity_summary.tsv so the two
## sit in one table. They are NOT recomputed here and NOT differenced against the
## GSE136831 betas -- different scales; see the header.
HLCA <- data.table(
    arm      = c("smad", "ieg", "smad", "ieg"),
    outcome  = c("bm_minus_fibrillar_z", "bm_minus_fibrillar_z",
                 "basement_membrane_score_z", "basement_membrane_score_z"),
    hlca_beta = c(-0.1869, -0.0138, -0.0025, -0.1867),
    hlca_null_mean = c(0.0263, 0.2240, 0.0472, 0.4495),
    hlca_empirical_p = c(0.0130, 0.000999, 0.5245, 0.000999))

fit_arm <- function(arm_col, arm_name, outcome, evaluable, why) {
    if (!evaluable)
        return(data.table(arm = arm_name, outcome = outcome, evaluable = FALSE,
                          n_units = nrow(d), estimate = NA_real_, SE = NA_real_,
                          t = NA_real_, p_value = NA_real_, mde_80pct = NA_real_,
                          verdict = why))
    rhs <- c(arm_col, "mean_log10_counts", if (nlevels(d$disease) > 1) "disease")
    f <- lm(reformulate(rhs, outcome), data = d)
    co <- summary(f)$coefficients[arm_col, ]
    ## The MDE converts "we found nothing" into the quantitative statement "this
    ## design could not have detected anything smaller than X". With 21 units it is
    ## the number that decides whether a null here means anything at all.
    data.table(arm = arm_name, outcome = outcome, evaluable = TRUE,
               n_units = nrow(f$model), estimate = unname(co[1]), SE = unname(co[2]),
               t = unname(co[3]), p_value = unname(co[4]),
               mde_80pct = 2.802 * unname(co[2]),
               verdict = NA_character_)
}

res <- rbindlist(list(
    fit_arm("smad_z", "smad", "bm_minus_fibrillar_z", TRUE, NA_character_),
    fit_arm("ieg_z",  "ieg",  "bm_minus_fibrillar_z", TRUE, NA_character_),
    fit_arm("smad_z", "smad", "bm_z", FALSE,
            "NOT EVALUABLE: the BM-alone endpoint has a matched null of +0.047, not 0, and this pseudobulk carries no background genes to rebuild one"),
    fit_arm("ieg_z",  "ieg",  "bm_z", FALSE,
            "NOT EVALUABLE: same reason; HLCA matched null +0.450")))

res <- merge(res, HLCA[, .(arm, outcome = fifelse(outcome == "basement_membrane_score_z", "bm_z", outcome),
                           hlca_beta, hlca_null_mean, hlca_empirical_p)],
             by = c("arm", "outcome"), all.x = TRUE, sort = FALSE)
res[, `:=`(
    same_sign_as_hlca = fifelse(is.na(estimate) | is.na(hlca_beta), NA,
                                sign(estimate) == sign(hlca_beta)),
    scale_note = "HLCA is z-within-dataset across 25 studies; this is z-within-compartment inside ONE study. Signs and significance only -- never difference the two betas.")]
wt(res, "smad_replication_gse136831.tsv")

cat("\n== SMAD arm replication in GSE136831 (signs and significance only) ==\n")
print(res[, .(arm, outcome, evaluable, n_units, estimate, p_value,
              hlca_beta, same_sign_as_hlca)])

## State the conclusion in the log rather than leaving it to be assembled.
sm <- res[arm == "smad" & outcome == "bm_minus_fibrillar_z"]
ie <- res[arm == "ieg"  & outcome == "bm_minus_fibrillar_z"]
cat("\n== VERDICT ==\n")
if (isTRUE(sm$same_sign_as_hlca) && sm$p_value < 0.05) {
    cat("The SMAD -> (BM - fibrillar) association HOLDS UP in GSE136831: same sign,\n",
        "p = ", signif(sm$p_value, 3), " on ", sm$n_units, " pericyte units.\n", sep = "")
    if (isTRUE(ie$p_value < 0.05))
        cat("BUT the IEG control is also significant here, so the arm SPECIFICITY does\n",
            "not replicate even though the direction does.\n", sep = "")
    else
        cat("The IEG control is null here (p = ", signif(ie$p_value, 3),
            "), so the arm specificity replicates too.\n", sep = "")
} else if (isTRUE(sm$same_sign_as_hlca)) {
    cat("Same sign as HLCA but not significant (p = ", signif(sm$p_value, 3),
        ") on ", sm$n_units, " units.\n", sep = "")
    ## Do not let "underpowered" be an assertion. Compare the MDE to the effect
    ## being replicated: if the MDE exceeds it, this dataset could not have
    ## detected the HLCA effect even if it were exactly true, and the null carries
    ## no information in either direction.
    if (is.finite(sm$mde_80pct) && sm$mde_80pct > abs(sm$hlca_beta))
        cat("THIS TEST WAS NOT INFORMATIVE. Its 80%-power MDE is ",
            signif(sm$mde_80pct, 3), ", larger than the HLCA effect it is meant to\n",
            "evaluate (", sm$hlca_beta, "). GSE136831 could not have detected an effect\n",
            "of that size even if it were exactly true, so this is neither support nor\n",
            "a failed replication -- it is an absent test. P1-20(a) is NOT resolved by\n",
            "this dataset, and the claim stays a supplementary observation.\n", sep = "")
    else
        cat("DIRECTIONALLY CONSISTENT and underpowered (MDE ", signif(sm$mde_80pct, 3),
            " vs HLCA ", sm$hlca_beta, "). Not the same as a failed replication, and\n",
            "not the same as support. The claim stays a supplementary observation.\n", sep = "")
} else {
    cat("OPPOSITE sign to HLCA (", signif(sm$estimate, 3), " vs -0.187), p = ",
        signif(sm$p_value, 3), " on ", sm$n_units, " units. The secondary-endpoint\n",
        "claim does NOT hold up and must not be promoted to a manuscript claim.\n", sep = "")
}
cat("\nEither way: 1 Control / 6 COPD / 15 IPF is a disease-skewed donor set, and\n",
    "no disease-stratified statement may be made from this file.\n", sep = "")
if (any(cov$n_genes_present < cov$n_genes_panel))
    cat("\nPANEL SHRINKAGE: ",
        paste(cov[n_genes_present < n_genes_panel,
                  sprintf("%s lost %s", panel, genes_missing)], collapse = "; "),
        ".\nA panel that lost genes is a different panel -- weigh the arm accordingly.\n",
        sep = "")

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
