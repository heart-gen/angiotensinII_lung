## Supplementary Table S6 -- ACTA2 benchmark (Part A) and CoGAPS validation (Part B).
##
## Part A re-emits the existing figures/mechanism/tableS_acta2_control.tsv through
## the manifest, joined to the pairwise contrast P values that table omitted.
##
## Part B assembles the CoGAPS rank-selection and validation outputs. Two labelling
## fixes are applied here because the source files cannot express them:
##   * "selected vs sensitivity" (nP=8 main, nP=9 sensitivity) existed only as
##     prose and as step_*.sh arguments, never as a column;
##   * the reconstruction-error column is `frac_unexplained`, NOT meanChiSq --
##     distributed CoGAPS does not populate meanChiSq, and mislabelling it would
##     invite a reviewer to compare it against published chi-square values.
##
## Per the table specification, the CoGAPS DISEASE associations
## (validation_np*/pattern_disease_ols.tsv.gz) are deliberately EXCLUDED: they
## belong with the disease Results section, in Table S13.
##
## Outputs: tsv/tableS06A.tsv, tableS06B1.tsv, tableS06B2.tsv, tableS06B3.tsv

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

## =========================================================================
## S6A -- ACTA2 benchmark
## =========================================================================
acta <- read_src(P("figures", "mechanism", "tableS_acta2_control.tsv"))
ph   <- read_src(P("pericyte_states", "_m", "stats_data",
                   "acta2_by_program_posthoc.tsv"))
if (!is.null(acta)) {
    if (!is.null(ph)) {
        ph2 <- ph[, .(panel = "A (contrasts)", readout = lens, program = contrast,
                      emmean = estimate, SE = SE,
                      pearson_p = p.value)]
        acta <- rbindlist(list(acta, ph2), fill = TRUE)
    }
    if ("pearson_p" %in% names(acta)) acta[, p_formatted := fmt_p(pearson_p)]
    write_part(acta, "06A",
        "ACTA2 benchmark: AGTR1 is not reducible to ACTA2+ contractile identity",
        supports = "Figure S4",
        sources = c("figures/mechanism/tableS_acta2_control.tsv",
                    "pericyte_states/_m/stats_data/acta2_by_program_posthoc.tsv"),
        notes = paste("Panel A rows are centred donor-aware program marginal means;",
                      "panel B rows are donor x program pseudobulk correlations",
                      "(n = 154). The 'A (contrasts)' rows are the pairwise program",
                      "contrasts, which the original source-data file omitted.",
                      "Both raw AGTR1-ACTA2 correlations vanish under denoising --",
                      "the signature of shared dropout, not shared biology."))
}

## =========================================================================
## S6B -- CoGAPS validation
## =========================================================================
CG <- function(...) P("pericyte_cogaps", "_m", ...)

## ---- B1: rank selection -------------------------------------------------
np <- read_src(CG("cogaps_nP_selection.tsv"))
if (!is.null(np)) {
    ## nP=8 is the main factorization and nP=9 the sensitivity run. This was only
    ## recoverable from prose and from which validation_np* directories exist.
    val_dirs <- basename(Sys.glob(CG("validation_np*")))
    val_np <- as.integer(sub("^validation_np", "", val_dirs))
    np[, selection_status := fcase(
        np == 8L, "SELECTED (main)",
        np == 9L, "sensitivity",
        np %in% val_np, "validated (secondary)",
        default = "swept only")]
    np[, validated := np %in% val_np]
    setnames(np, "frac_unexplained", "reconstruction_error_frac_unexplained")
    setnames(np, c("mean_r", "min_r"),
             c("cross_seed_r_mean", "cross_seed_r_weakest_pattern"))
    write_part(np, "06B1",
        "CoGAPS rank selection: cross-seed pattern reproducibility and reconstruction error",
        supports = "Figure S8",
        sources = "pericyte_cogaps/_m/cogaps_nP_selection.tsv",
        notes = paste("`cross_seed_r_weakest_pattern` is the minimum across",
                      "patterns, i.e. the least reproducible pattern at that rank.",
                      "`n_seeds` = 3 MATCHED replicate seeds; the canonical seed",
                      "(13) is the reference each replicate is matched to, not a",
                      "fourth comparison. Reconstruction error is reported as",
                      "fraction-unexplained: distributed CoGAPS does not populate",
                      "meanChiSq, so this is NOT a chi-square statistic.",
                      "Validation was run for nP = 5, 7, 8, 9 only."))
}

## ---- B2: per-pattern cross-seed stability -------------------------------
ss <- read_src(CG("cogaps_seed_stability_summary.tsv"))
if (!is.null(ss)) {
    setnames(ss, c("mean_r", "min_r"), c("cross_seed_r_mean", "cross_seed_r_min"))
    write_part(ss, "06B2",
        "CoGAPS per-pattern cross-seed reproducibility at the selected rank",
        supports = "Figure S8",
        sources = "pericyte_cogaps/_m/cogaps_seed_stability_summary.tsv",
        notes = "Best-match correlation of each reference pattern to each replicate seed.")
}

## ---- B3: pattern-to-program correspondence at each validated rank -------
bits <- list()
for (d in Sys.glob(CG("validation_np*"))) {
    npv <- sub("^validation_np", "", basename(d))
    ov <- read_src(file.path(d, "pattern_panel_overlap.tsv.gz"))
    if (!is.null(ov)) {
        ov[, `:=`(nP = as.integer(npv), block = "PatternMarker x panel overlap")]
        setnames(ov, c("phyper", "padj"), c("hypergeometric_p", "hypergeometric_p_BH"),
                 skip_absent = TRUE)
        bits[[paste0("ov", npv)]] <- ov
    }
    sc <- read_src(file.path(d, "pattern_score_spearman.tsv.gz"))
    if (!is.null(sc)) {
        sc[, `:=`(nP = as.integer(npv), block = "pattern vs program score (Spearman)")]
        bits[[paste0("sc", npv)]] <- sc
    }
    ag <- read_src(file.path(d, "pattern_AGTR1_spearman.tsv.gz"))
    if (!is.null(ag)) {
        ag[, `:=`(nP = as.integer(npv), block = "pattern vs AGTR1 (Spearman)")]
        setnames(ag, "rho_AGTR1", "rho", skip_absent = TRUE)
        ag[, score := "AGTR1"]
        bits[[paste0("ag", npv)]] <- ag
    }
}
if (length(bits)) {
    b3 <- rbindlist(bits, fill = TRUE)
    setcolorder(b3, intersect(c("nP", "block", "pattern", "program", "score", "rho",
                                "n_markers", "panel_size", "overlap", "jaccard",
                                "hypergeometric_p", "hypergeometric_p_BH"),
                              names(b3)))
    setorder(b3, nP, block, pattern)
    write_part(b3, "06B3",
        "CoGAPS pattern-to-program correspondence: score correlations and PatternMarker overlaps",
        supports = "Figure S8",
        sources = c("pericyte_cogaps/_m/validation_np{5,7,8,9}/pattern_panel_overlap.tsv.gz",
                    ".../pattern_score_spearman.tsv.gz",
                    ".../pattern_AGTR1_spearman.tsv.gz"),
        notes = paste("Hypergeometric P tests overlap between a pattern's",
                      "PatternMarker set and a curated panel; BH is across",
                      "pattern x panel pairs within a rank, as adjusted by the",
                      "source script. Disease associations of the CoGAPS patterns",
                      "are deliberately NOT included here -- they belong with the",
                      "disease Results section (Table S13)."))
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
