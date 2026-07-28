## Supplementary Tables S4 (stable pericyte cluster characteristics) and
## S5 (AGTR1 dropout, airspace-affinity, and the three-lens program models).
##
## S4 is mostly annotations/state_summary.tsv, plus two columns that existed
## nowhere: per-cluster bootstrap Jaccard (00.state_discovery.py computed it but
## persisted only the median across clusters) and mean scVI-denoised AGTR1 per
## cluster (03.agtr1_lenses.R does the barcode merge but only ever aggregates by
## PROGRAM, never by cluster).
##
## S5 assembles the dropout model, the airspace models, and the three-lens
## program models, and adds the BH adjustment the airspace outputs never carried.
##
## Outputs: tsv/tableS04.tsv, tableS05A.tsv, tableS05B.tsv, tableS05C*.tsv

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

## =========================================================================
## S4 -- stable cluster characteristics
## =========================================================================
ss <- read_src(P("pericyte_states", "_m", "annotations", "state_summary.tsv"))

## ---- denoised AGTR1 per cluster (new) -----------------------------------
## Same merge as pericyte_states/_h/03.agtr1_lenses.R:42-48 -- metadata barcode
## joined to the pericyte-only scVI model -- but grouped by pericyte_state.
DEN_MODEL <- "Pericyte-only-trained"
den_by_cluster <- NULL
meta_f <- P("pericyte_states", "_m", "pericytes_states_metadata.tsv.gz")
den_f  <- P("localization", "airspace_analysis", "_m", "airspace",
            "pericytes_airspace_denoising.tsv")
if (file.exists(meta_f) && file.exists(den_f)) {
    meta <- fread(meta_f); setnames(meta, 1, "barcode")
    den  <- fread(den_f, select = c("index", "Model", "AGTR1_scvi"))
    setnames(den, "index", "barcode")
    den  <- den[Model == DEN_MODEL]
    stopifnot(!anyDuplicated(den$barcode))
    mm <- merge(meta[, .(barcode, pericyte_state, AGTR1_expr, AGTR1_detect)],
                den[, .(barcode, AGTR1_scvi)], by = "barcode")
    cat(sprintf("denoised merge: %d of %d cells (%.1f%%)\n",
                nrow(mm), nrow(meta), 100 * nrow(mm) / nrow(meta)))
    den_by_cluster <- mm[, .(
        n_cells_denoised   = .N,
        AGTR1_scvi_mean    = mean(AGTR1_scvi, na.rm = TRUE),
        AGTR1_scvi_median  = as.numeric(median(AGTR1_scvi, na.rm = TRUE)),
        AGTR1_raw_mean_chk = mean(AGTR1_expr, na.rm = TRUE),
        AGTR1_detect_chk   = mean(AGTR1_detect, na.rm = TRUE)
    ), by = .(pericyte_state)]
}

## ---- per-cluster bootstrap Jaccard (from the --stability-only re-run) ----
jac_f <- P("pericyte_states", "_m", "stability", "cluster_bootstrap_jaccard.tsv")
jac <- read_src(jac_f)
s4_status <- "complete"
if (is.null(jac)) {
    s4_status <- "pending_upstream"
    cat("NOTE: cluster_bootstrap_jaccard.tsv absent -- run\n",
        "  sbatch -D pericyte_states/_m pericyte_states/_h/step_0c.sh\n")
}

if (!is.null(ss)) {
    s4 <- copy(ss)
    setnames(s4, "pericyte_state", "cluster_id")
    ## Cluster ids are integers in state_summary.tsv but strings in the metadata
    ## and the Jaccard table; join on character throughout.
    s4[, cluster_id := as.character(cluster_id)]
    if (!is.null(den_by_cluster)) {
        den_by_cluster[, cluster_id := as.character(pericyte_state)][
            , pericyte_state := NULL]
        s4 <- merge(s4, den_by_cluster, by = "cluster_id", all.x = TRUE)
    }
    if (!is.null(jac)) {
        j <- jac[!is.na(pericyte_state),
                 .(cluster_id = as.character(pericyte_state), n_boot,
                   jaccard_mean, jaccard_median, jaccard_min,
                   jaccard_lo, jaccard_hi)]
        s4 <- merge(s4, j, by.x = "cluster_id", by.y = "cluster_id", all.x = TRUE)
    }
    setnames(s4, "AGTR1_pos_frac", "AGTR1_detect_frac", skip_absent = TRUE)
    setnames(s4, "AGTR1_mean", "AGTR1_raw_mean", skip_absent = TRUE)
    setcolorder(s4, intersect(c("cluster_id", "state_program", "n_cells", "n_donors",
                                "n_studies", "n_datasets", "jaccard_mean",
                                "jaccard_median", "jaccard_lo", "jaccard_hi",
                                "AGTR1_detect_frac", "AGTR1_raw_mean",
                                "AGTR1_scvi_mean"), names(s4)))
    write_part(s4, "04",
        "Stable pericyte cluster characteristics and program assignments",
        supports = "Figure 2; Figure S7",
        sources = c("pericyte_states/_m/annotations/state_summary.tsv",
                    "pericyte_states/_m/stability/cluster_bootstrap_jaccard.tsv",
                    "localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv"),
        status = s4_status,
        notes = paste("Clusters are bootstrap-stable Leiden clusters on the",
                      "study-integrated embedding, not marker-defined.",
                      "`*_relenrich` is enrichment of a program score relative to",
                      "the other clusters. Denoised AGTR1 is the",
                      DEN_MODEL, "scVI model.",
                      ## The published median (0.966) is not representative: it is
                      ## a median over clusters, and reporting it alone would let a
                      ## reader assume every cluster is that reproducible.
                      if (!is.null(jac) && "jaccard_mean" %in% names(s4))
                          sprintf(paste("Per-cluster mean bootstrap Jaccard ranges",
                                        "%.2f-%.2f across the %d clusters; the",
                                        "solution-level median of 0.966 conceals",
                                        "this spread, so per-cluster values should",
                                        "be cited rather than the median alone."),
                                  min(s4$jaccard_mean, na.rm = TRUE),
                                  max(s4$jaccard_mean, na.rm = TRUE),
                                  sum(!is.na(s4$jaccard_mean))) else "",
                      if (s4_status != "complete")
                          "PENDING: per-cluster Jaccard requires step_0c.sh." else ""))
}

## =========================================================================
## S5A -- AGTR1 dropout model
## =========================================================================
dr <- read_src(P("localization", "pericyte_analysis", "_m", "qc_results",
                 "dropout_expectation_results.tsv"))
if (!is.null(dr)) {
    setnames(dr, c("obs_zeros", "exp_zeros", "ratio_obs_to_exp", "pval"),
             c("observed_undetected", "expected_undetected",
               "ratio_observed_to_expected", "p_value"), skip_absent = TRUE)
    dr[, p_formatted := fmt_p(p_value)]
    write_part(dr, "05A",
        "AGTR1 dropout model: observed vs depth-matched expected undetected cells",
        supports = "Figure 1; Results (localization)",
        sources = "localization/pericyte_analysis/_m/qc_results/dropout_expectation_results.tsv",
        notes = paste("Expected zeros from depth- and expression-matched control",
                      "genes (FFT method, n_matched genes). The model is fit for",
                      "AGTR1 only -- no AGTR2 or ACTA2 row exists.",
                      "A ratio near 1 with a non-significant P means AGTR1 is NOT",
                      "excessively dropped out relative to matched transcripts."))
}

## =========================================================================
## S5B -- airspace-affinity models
## =========================================================================
ab <- list()
pc_f <- P("localization", "airspace_analysis", "_m", "airspace",
          "airspace_effect_AGTR1_per_cell.csv")
dn_f <- P("localization", "airspace_analysis", "_m", "airspace",
          "airspace_effect_AGTR1.csv")
pcm <- read_src(pc_f); dnm <- read_src(dn_f)
if (!is.null(pcm)) ab$cell <- pcm[, .(
    level = "cell", predictor = term, estimate, se, ci_lower, ci_upper,
    p_value = pval, n_cells, n_donors, formula, fallback)]
if (!is.null(dnm)) ab$donor <- dnm[, .(
    level = "donor", predictor = term, estimate, se,
    ci_lower = estimate - 1.96 * se, ci_upper = estimate + 1.96 * se,
    p_value = pval, n_cells = NA_integer_, n_donors = 56L,
    formula = "airspace_score ~ frac_AGTR1_pos + age + C(sex)",
    fallback = NA)]
if (length(ab)) {
    s5b <- rbindlist(ab, fill = TRUE)
    ## Neither source file carried an adjusted P. The two models are the same
    ## hypothesis at two units of analysis, so BH across them is the honest family.
    s5b[, p_BH := bh(p_value)]
    s5b[, p_formatted := fmt_p(p_value)]
    write_part(s5b, "05B",
        "Airspace-affinity models for AGTR1 (cell level and donor level)",
        supports = "Figure 1; Figure S1",
        sources = c(sub(paste0(ROOT, "/"), "", pc_f), sub(paste0(ROOT, "/"), "", dn_f)),
        notes = paste("BH ADDED HERE: neither source file carried an adjusted P.",
                      "The donor-level CI is derived as estimate +/- 1.96*SE; the",
                      "source CSV stored only estimate, se and p. Full model",
                      "summaries with all covariate rows are in",
                      "airspace_model_summary.txt and",
                      "airspace_lmm_per_cell_summary.txt."))
}

## =========================================================================
## S5C -- program models under three lenses
## =========================================================================
SD <- function(f) P("pericyte_states", "_m", "stats_data", f)
lens_emm <- read_src(SD("agtr1_lenses_by_program_emmeans.tsv"))
lens_ph  <- read_src(SD("agtr1_lenses_by_program_posthoc.tsv"))
lens_top <- read_src(SD("agtr1_lenses_top_program.tsv"))

LENS_LAB <- c(AGTR1_expr = "raw expression", AGTR1_detect = "binary detection",
              AGTR1_scvi = "scVI-denoised")
if (!is.null(lens_emm)) {
    lens_emm[, lens_label := LENS_LAB[lens]]
    write_part(lens_emm, "05C1",
        "AGTR1 by pericyte program under three measurement lenses: marginal means",
        supports = "Figure 2; Figure S5",
        sources = "pericyte_states/_m/stats_data/agtr1_lenses_by_program_emmeans.tsv",
        notes = paste("Donor-aware mixed model per lens. The program ranking",
                      "REVERSES between the raw and denoised lenses, which is why",
                      "AGTR1 is reported as a compartment label and not a state",
                      "marker."))
}
if (!is.null(lens_ph)) {
    lens_ph[, lens_label := LENS_LAB[lens]]
    lens_ph[, p_formatted := fmt_p(p.value)]
    write_part(lens_ph, "05C2",
        "AGTR1 by pericyte program under three lenses: pairwise contrasts",
        supports = "Figure 2; Figure S5",
        sources = "pericyte_states/_m/stats_data/agtr1_lenses_by_program_posthoc.tsv",
        notes = "P values are BH-adjusted within each lens by the source script.")
}
if (!is.null(lens_top)) {
    lens_top[, lens_label := LENS_LAB[lens]]
    write_part(lens_top, "05C3",
        "Top-ranked pericyte program under each AGTR1 measurement lens",
        supports = "Figure 2",
        sources = "pericyte_states/_m/stats_data/agtr1_lenses_top_program.tsv",
        notes = "The three-row summary of the lens reversal.")
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
