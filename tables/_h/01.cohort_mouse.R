## Supplementary Tables S1 (cohort/dataset characteristics) and S2 (mouse mural
## composition and raw Agtr1a detection).
##
## S1 has no precursor anywhere in the repo -- no cohort_summary.tsv exists -- so
## it is computed here in one pass over the per-cell pericyte metadata, with the
## non-HLCA cohorts (mouse, IPF, GTEx, LungMAP) added as their own rows.
##
## S2 is assembly, EXCEPT Part C's depth dependence: the median UMI per pericyte
## and the Agtr1a-vs-depth correlation are quoted in writings/SPECIES_SUMMARY.md
## and MECHANISM_ANALYSES.md but were never written to any file. They are
## recomputed here from species_comparability_pericyte_cells.tsv.
##
## Outputs: tsv/tableS01.tsv, tableS02A.tsv, tableS02B.tsv, tableS02C_*.tsv

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

## =========================================================================
## S1 -- cohort and dataset characteristics
## =========================================================================
meta <- read_src(P("pericyte_states", "_m", "pericytes_states_metadata.tsv.gz"))

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc) ~ "Healthy",
        lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc,
              ignore.case = TRUE) ~ "Fibrotic_ILD",
        TRUE ~ "Other")
}

if (!is.null(meta)) {
    meta[, age := suppressWarnings(as.numeric(age_or_mean_of_age_range))]
    meta[, disease_group := map_disease_group(lung_condition)]

    ## Per-donor first, so donor-level attributes are not weighted by cell count.
    donor <- meta[, .(n_cells = .N, dataset = first(dataset), study = first(study),
                      disease_group = first(disease_group),
                      sex = first(sex), age = mean(age, na.rm = TRUE),
                      smoking = first(smoking_status), bmi = first(BMI)),
                  by = donor_id]

    pct <- function(x) round(100 * mean(x), 1)
    ## min/max over an all-NA age vector returns Inf/-Inf with a warning, which
    ## then prints as a literal "Inf" age in the table. Several HLCA datasets
    ## deposit no age at all, so guard rather than emit a nonsense bound.
    safe <- function(f, x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else round(f(x), 1) }
    s1_hlca <- donor[, .(
        n_donors      = .N,
        n_studies     = uniqueN(study),
        n_cells       = sum(n_cells),
        cells_median  = as.numeric(median(n_cells)),
        cells_min     = min(n_cells),
        cells_max     = max(n_cells),
        n_healthy     = sum(disease_group == "Healthy"),
        n_copd        = sum(disease_group == "COPD"),
        n_fibrotic    = sum(disease_group == "Fibrotic_ILD"),
        n_other       = sum(disease_group == "Other"),
        age_median    = safe(median, age),
        age_min       = safe(min, age),
        age_max       = safe(max, age),
        pct_age_known = pct(!is.na(age)),
        pct_female    = pct(tolower(as.character(sex)) %in% c("female", "f")),
        pct_smoking_known = pct(!is.na(smoking) & nzchar(as.character(smoking)) &
                                    !as.character(smoking) %in% c("nan", "NA", "unknown")),
        pct_bmi_known = pct(!is.na(bmi))
    ), by = .(dataset)][order(-n_donors)]
    s1_hlca[, cohort := "HLCA (pericytes)"][, species := "human"]
    ## The cohorts do not carry equivalent metadata, and blank cells in a
    ## supplementary table read as an omission unless the table says otherwise.
    ## These two columns state, per cohort, what n_cells counts and which
    ## donor attributes the source actually deposits.
    s1_hlca[, cell_count_basis := "pericytes"]
    s1_hlca[, donor_metadata_available := paste(
        "age; sex; smoking status; BMI; ethnicity; lung condition (disease group)")]

    ## Donors at each analysis threshold. These two denominators were previously
    ## produced by different modules (>=10 in disease_association, >=20 in
    ## pericyte_states) and never reported side by side, so which donors an
    ## analysis actually used was not auditable.
    thr <- donor[, .(cohort = "HLCA (pericytes)", dataset = "ALL DATASETS",
                     n_donors = .N,
                     n_donors_ge10 = sum(n_cells >= 10),
                     n_donors_ge20 = sum(n_cells >= 20),
                     n_donors_age_sex = sum(!is.na(age) & !is.na(sex)),
                     n_cells = sum(n_cells), n_studies = uniqueN(study))]

    s1 <- rbindlist(list(s1_hlca, thr), fill = TRUE)
} else s1 <- NULL

## ---- non-HLCA cohorts ---------------------------------------------------
extra <- list()

sp <- read_src(P("cross_species", "_m", "stats_data",
                 "species_comparability_summary.tsv"))
if (!is.null(sp)) {
    g <- function(k) { v <- sp[metric == k, value]; if (length(v)) as.numeric(v[1]) else NA_real_ }
    extra$mouse <- data.table(
        cohort = "Mouse lung (CELLxGENE)", species = "mouse", dataset = "4 datasets",
        n_donors = g("n_pericyte_donors"), n_cells = g("n_mural_cells"),
        n_studies = g("n_pericyte_datasets"),
        cell_count_basis = "mural cells",
        donor_metadata_available =
            "dataset of origin only -- no age, sex or disease",
        notes = paste("mural cells; pericytes =", g("n_pericytes"),
                      "| donor_id is 'pooled' in 2 of 4 datasets"))
}

ipf <- read_src(P("disease_association", "ipf_analysis", "_h", "sample_demo.csv"))
if (!is.null(ipf)) {
    nm <- names(ipf)
    sexc <- grep("^Sex$", nm, value = TRUE, ignore.case = TRUE)
    agec <- grep("^Age", nm, value = TRUE, ignore.case = TRUE)
    grp  <- grep("Disease", nm, value = TRUE, ignore.case = TRUE)
    extra$ipf <- data.table(
        cohort = "IPF cohort (GEO)", species = "human", dataset = "sample_demo",
        n_donors = uniqueN(ipf[[1]]),
        pct_female = if (length(sexc)) round(100 * mean(tolower(ipf[[sexc[1]]]) %in%
                                                        c("female", "f")), 1) else NA_real_,
        age_median = if (length(agec)) round(median(suppressWarnings(
            as.numeric(ipf[[agec[1]]])), na.rm = TRUE), 1) else NA_real_,
        cell_count_basis = "bulk RNA-seq libraries (no single cells)",
        donor_metadata_available =
            "age; sex; race; ever-smoker; disease group -- no BMI",
        notes = if (length(grp)) paste("groups:", paste(sort(unique(
            as.character(ipf[[grp[1]]]))), collapse = "/")) else NA_character_)
}

gtex <- read_src(P("inputs", "gtex", "_m", "gtex_v11_lung_sample_data.tsv"))
if (!is.null(gtex)) {
    extra$gtex <- data.table(
        cohort = "GTEx v11 lung", species = "human", dataset = "GTEx",
        n_donors = uniqueN(gtex$SUBJID), n_cells = NA_real_,
        pct_female = round(100 * mean(unique(gtex[, .(SUBJID, SEX)])$SEX == 2), 1),
        cell_count_basis = "bulk lung samples (no single cells)",
        donor_metadata_available = paste(
            "sex; age BRACKET (decade, not exact age); race; death circumstance",
            "-- no smoking or BMI"),
        notes = paste(nrow(gtex), "lung samples; bulk, sample-level metadata only.",
                      "age_median is left NA because GTEx releases AGE as a",
                      "decade bracket, which has no median comparable to the",
                      "exact ages in the HLCA rows."))
}

## ---- LungMAP ------------------------------------------------------------
## LungMAP (GSE161382) is deliberately NOT filled in to the depth of the HLCA
## rows, because the source does not support it. The deposited metadata carries
## donor, age, sex and batch and nothing else -- no smoking history, no BMI, no
## ethnicity, and no disease field. Those columns stay NA and are named in
## `donor_metadata_available` so the blanks read as a property of GSE161382
## rather than as missing work here.
##
## Two further asymmetries are made explicit rather than papered over:
##   * Age is deposited as strings on two scales -- "31wk" (post-conceptional
##     weeks) and "3yr"/"31yr". This is a developmental atlas spanning preterm
##     to young adult, so its age distribution is not comparable to the HLCA
##     adult distribution even after conversion to years.
##   * Per-donor pericyte counts come from the DEPOSITED annotation (803 cells),
##     which is the only pericyte set resolvable per donor from plain-text
##     inputs. The replication analysis instead used transferred HLCA labels
##     (1,051 predicted, 762 retained after confidence/purity filtering); both
##     counts are recorded so neither can be mistaken for the other.
lm_meta <- read_src(P("inputs", "lungmap", "_m", "GSE161382_metadata.txt.gz"))
lm_f    <- read_src(P("localization", "lungmap_replication", "_m",
                      "annotated_clusters_filtering_summary.tsv"))
if (!is.null(lm_meta) || !is.null(lm_f)) {
    n_all_dep <- n_all_ret <- n_peri_pred <- n_peri_ret <- NA_real_
    if (!is.null(lm_f)) {
        n_all_dep   <- sum(suppressWarnings(as.numeric(lm_f$Original)), na.rm = TRUE)
        n_all_ret   <- sum(suppressWarnings(as.numeric(lm_f$Filtered)), na.rm = TRUE)
        pr          <- lm_f[Cell_Type == "Pericytes"]
        if (nrow(pr)) { n_peri_pred <- as.numeric(pr$Original[1])
                        n_peri_ret  <- as.numeric(pr$Filtered[1]) }
    }

    lm_row <- data.table(cohort = "LungMAP", species = "human",
                         dataset = "GSE161382", n_studies = 1L,
                         cell_count_basis = "pericytes (deposited annotation)")

    if (!is.null(lm_meta)) {
        yrs <- function(s) {
            n <- suppressWarnings(as.numeric(sub("(wk|yr)$", "", s)))
            ifelse(grepl("wk$", s), n / 52, n)
        }
        dm <- lm_meta[, .(n_cells_all = .N,
                          n_peri = sum(celltype == "pericytes"),
                          age = first(age), sex = first(sex)), by = donor]
        dm[, age_years := yrs(age)]
        lm_row[, `:=`(
            n_donors     = nrow(dm),
            n_cells      = sum(dm$n_peri),
            cells_median = as.numeric(median(dm$n_peri)),
            cells_min    = min(dm$n_peri), cells_max = max(dm$n_peri),
            age_median   = round(median(dm$age_years), 2),
            age_min      = round(min(dm$age_years), 2),
            age_max      = round(max(dm$age_years), 2),
            pct_age_known     = 100,
            pct_female        = round(100 * mean(dm$sex == "F"), 1),
            pct_smoking_known = 0, pct_bmi_known = 0,
            donor_metadata_available = paste0(
                "age (", paste(sort(unique(dm$age)), collapse = "/"),
                "); sex; batch -- NO smoking, BMI, ethnicity or disease field"))]
    } else {
        lm_row[, donor_metadata_available :=
                   "not read (GSE161382_metadata.txt.gz missing)"]
    }

    lm_row[, notes := paste0(
        "Developmental atlas: ages span preterm to young adult and are deposited ",
        "on two scales (weeks/years), converted to years here; the age columns ",
        "are NOT comparable to the HLCA adult distribution. ",
        "Disease columns are NA because GSE161382 deposits no disease field ",
        "(all donors are non-diseased by study design). ",
        "n_cells and cells_* are pericytes under the DEPOSITED annotation; the ",
        "replication analysis used transferred HLCA labels (",
        n_peri_pred, " predicted Pericytes, ", n_peri_ret,
        " retained after confidence/purity filtering, out of ", n_all_dep,
        " deposited / ", n_all_ret, " retained cells of all types).")]
    extra$lungmap <- lm_row
}

if (!is.null(s1)) {
    s1 <- rbindlist(c(list(s1), extra), fill = TRUE)
    setcolorder(s1, c("cohort", "species", "dataset", "cell_count_basis",
                      "n_donors", "n_studies", "n_cells"))
    ## Push the two explanatory columns to the end, next to each other: a reader
    ## reaches them only after finding a blank numeric cell.
    tail_cols <- intersect(c("donor_metadata_available", "notes"), names(s1))
    setcolorder(s1, c(setdiff(names(s1), tail_cols), tail_cols))
    write_part(s1, "01",
        "Cohort and dataset characteristics",
        supports = "Methods; all analyses",
        sources = c("pericyte_states/_m/pericytes_states_metadata.tsv.gz",
                    "cross_species/_m/stats_data/species_comparability_summary.tsv",
                    "disease_association/ipf_analysis/_h/sample_demo.csv",
                    "inputs/gtex/_m/gtex_v11_lung_sample_data.tsv",
                    "inputs/lungmap/_m/GSE161382_metadata.txt.gz",
                    "localization/lungmap_replication/_m/annotated_clusters_filtering_summary.tsv"),
        notes = paste("Percentages are donor-level, not cell-weighted.",
                      "The ALL DATASETS row carries the analysis denominators:",
                      ">=10 pericytes is the primary donor filter, >=20 the",
                      "sensitivity filter. The cohorts do NOT carry equivalent",
                      "metadata: `cell_count_basis` states what n_cells counts and",
                      "`donor_metadata_available` names the donor attributes each",
                      "source deposits, so a blank cell can be read as absent at",
                      "source rather than as an omission. LungMAP in particular",
                      "deposits only age, sex and batch, on a developmental age",
                      "scale that is not comparable to the HLCA adult donors."))
}

## =========================================================================
## S2 -- mouse mural composition and raw Agtr1a detection
## =========================================================================
SP <- function(f) P("cross_species", "_m", "stats_data", paste0("species_comparability_", f))

## ---- Part A: composition and comparability ------------------------------
a_bits <- list()
if (!is.null(sp)) a_bits$summary <- sp[, .(block = "headline metric", key = metric,
                                           value = as.character(value))]
ct_ds <- read_src(SP("celltype_by_dataset.tsv"))
if (!is.null(ct_ds)) {
    m <- melt(ct_ds, id.vars = names(ct_ds)[1], variable.name = "key",
              value.name = "value")
    setnames(m, names(m)[1], "cell_type")
    a_bits$ct <- m[, .(block = "cells by dataset", key = paste0(cell_type, " | ", key),
                       value = as.character(value))]
}
st_ct <- read_src(SP("state_by_celltype.tsv"))
if (!is.null(st_ct)) {
    m <- melt(st_ct, id.vars = names(st_ct)[1], variable.name = "key",
              value.name = "value")
    setnames(m, names(m)[1], "pericyte_state")
    a_bits$st <- m[, .(block = "human-program label by mouse cell type",
                       key = paste0(pericyte_state, " | ", key), value = as.character(value))]
}
by_donor <- read_src(SP("pericyte_by_donor.tsv"))
if (!is.null(by_donor))
    a_bits$donor <- by_donor[, .(block = "pericytes per donor",
                                 key = paste0(dataset_id, " | ", donor_id),
                                 value = paste0(n_pericytes, " pericytes, ",
                                                n_agtr1a_pos, " Agtr1a+"))]
if (length(a_bits))
    write_part(rbindlist(a_bits, fill = TRUE), "02A",
        "Mouse lung mural-cell composition and dataset comparability",
        supports = "Figure S2",
        sources = c("cross_species/_m/stats_data/species_comparability_summary.tsv",
                    "..._celltype_by_dataset.tsv", "..._state_by_celltype.tsv",
                    "..._pericyte_by_donor.tsv"),
        notes = paste("The mouse mural set is overwhelmingly smooth muscle;",
                      "human program labels transferred onto it are SMC labels and",
                      "support no state-level claim. Outputs of",
                      "cross_species/_h/03.conservation_stats.R are superseded and",
                      "deliberately excluded."))

## ---- Part B: within-dataset detection tests -----------------------------
tests <- read_src(SP("agtr1a_tests.tsv"))
within <- read_src(SP("agtr1a_within_dataset.tsv"))
if (!is.null(tests)) {
    tests[, p_BH := bh(p_value)]
    tests[, p_formatted := fmt_p(p_value)]
    write_part(tests, "02B",
        "Within-dataset Agtr1a detection tests (pericytes vs smooth muscle)",
        supports = "Figure S2",
        sources = "cross_species/_m/stats_data/species_comparability_agtr1a_tests.tsv",
        notes = paste("Raw counts. Fisher exact within each dataset that contains",
                      "both cell types, plus a Mantel-Haenszel test stratified by",
                      "dataset. BH is across the two Fisher tests."))
}

## ---- Part C1: depth dependence (RECOMPUTED -- no prior source file) ------
pc <- read_src(SP("pericyte_cells.tsv"))
if (!is.null(pc)) {
    ## Two distinct depth relationships, kept apart because they answer different
    ## questions and give different numbers:
    ##   detection ~ depth  is the manuscript's claim ("detection tracks depth"),
    ##                      a point-biserial-style Spearman on the 0/1 indicator.
    ##   counts ~ depth     is the supporting magnitude relationship.
    ## Quoting one for the other would misstate the result, so both are labelled.
    depth_stats <- function(d, label) {
        ok  <- is.finite(d$total_counts_cell) & is.finite(d$Agtr1a_counts)
        det <- as.numeric(d$Agtr1a_counts[ok] > 0)
        umi <- d$total_counts_cell[ok]
        one <- function(y, measure) {
            ct <- if (sum(ok) >= 4 && uniqueN(y) > 1)
                suppressWarnings(cor.test(umi, y, method = "spearman")) else NULL
            data.table(
                dataset_id = label, measure = measure, n_pericytes = sum(ok),
                median_umi_per_cell = as.numeric(median(umi)),
                min_umi = min(umi), max_umi = max(umi),
                n_agtr1a_pos = sum(det), detect_frac = mean(det),
                median_agtr1a_umi_in_pos = as.numeric(median(
                    d$Agtr1a_counts[ok][det == 1])),
                spearman_rho = if (is.null(ct)) NA_real_ else unname(ct$estimate),
                spearman_p   = if (is.null(ct)) NA_real_ else ct$p.value)
        }
        rbind(one(det, "Agtr1a detection (0/1) ~ total UMI"),
              one(d$Agtr1a_counts[ok], "Agtr1a raw counts ~ total UMI"))
    }

    depth <- rbindlist(c(
        lapply(split(pc, pc$dataset_id), function(d) depth_stats(d, d$dataset_id[1])),
        list(depth_stats(pc, "ALL DATASETS POOLED"))), fill = TRUE)
    setorder(depth, measure, dataset_id)
    depth[, p_formatted := fmt_p(spearman_p)]
    write_part(depth, "02C1",
        "Sequencing-depth dependence of raw Agtr1a detection in mouse pericytes",
        supports = "Figure S2; Results (cross-species)",
        sources = "cross_species/_m/stats_data/species_comparability_pericyte_cells.tsv",
        notes = paste("RECOMPUTED for this table: these statistics were quoted in",
                      "writings/SPECIES_SUMMARY.md and MECHANISM_ANALYSES.md but",
                      "were never written to any output file. The reported",
                      "'detection tracks depth' result is the pooled",
                      "detection-vs-UMI row (rho = 0.432, P = 0.0048, n = 41);",
                      "the counts-vs-UMI rows are a different, stronger",
                      "relationship and must not be quoted in its place."))
}

## ---- Part C2: receptor specificity --------------------------------------
lay <- read_src(SP("agtr1a_by_layer.tsv"))
if (!is.null(lay)) {
    setorder(lay, gene, layer, cell_type)
    write_part(lay, "02C2",
        "Receptor specificity of the mouse mural signal (Agtr1a/Agtr1b/Agtr2, raw vs scVI)",
        supports = "Figure S2",
        sources = "cross_species/_m/stats_data/species_comparability_agtr1a_by_layer.tsv",
        notes = paste("The scVI-corrected layer is dense: every cell is 'detected',",
                      "which is why the compartment claim rests on the raw counts."))
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
