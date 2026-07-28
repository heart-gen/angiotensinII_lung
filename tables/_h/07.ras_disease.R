## Supplementary Tables S12 (local renin-angiotensin machinery), S13 (disease
## associations of discrete composition and continuous injury-stromal programs)
## and S14 (net niche index, smoking availability, leave-one-study-out).
##
## Threshold convention throughout S13/S14: >=10 pericytes per donor is PRIMARY
## and lives in the unsuffixed source files; >=20 is the sensitivity analysis in
## the `_mincells20` files. Both are stacked with a `min_cells` column so the two
## denominators can never be conflated the way they were when the modules
## disagreed silently.
##
## Composition and injury-fraction posthoc P values are ALREADY BH-adjusted by
## emmeans::pairs(adjust = "BH") in the source scripts and are relabelled, not
## re-adjusted.
##
## Outputs: tsv/tableS12A.tsv .. S12C, tableS13A.tsv .. S13F, tableS14A.tsv .. S14C

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

AG  <- function(...) P("agt_axis", "_m", ...)
PSS <- function(f) P("pericyte_states", "_m", "stats_data", f)
NI  <- function(f) P("niche_index", "_m", "stats_data", f)
SE  <- function(f) P("sensitivity", "_m", "stats_data", f)
DA  <- function(f) P("disease_association", "_m", "mixed_model_forest", f)

## =========================================================================
## S12 -- local renin-angiotensin machinery
## =========================================================================
prof  <- read_src(AG("stats_data", "ras_celltype_profile.tsv"))
panel <- read_src(P("agt_axis", "_h", "ras_panel.tsv"))
rpbk  <- read_src(AG("ras_pseudobulk_celltype.tsv.gz"))

CATEGORY <- c(ras_substrate = "substrate",
              ras_protease = "angiotensin II-generating or degradative enzyme",
              ras_receptor = "receptor",
              comparator_ligand = "comparator ligand")

if (!is.null(prof)) {
    s12a <- copy(prof)
    if (!is.null(panel)) {
        s12a <- merge(s12a, panel, by = "gene", all.x = TRUE)
        s12a[, functional_category := CATEGORY[panel]]
    }
    ## Full within-gene ranking: ras_top_celltypes.tsv keeps only the top 3.
    s12a[, within_gene_rank := frank(-emmean, ties.method = "min"), by = gene]
    if (!is.null(rpbk)) {
        ## Per-cell-type denominators under the model's own filter
        ## (01.ras_landscape_stats.R: n_cells >= 5, then min_donors >= 5).
        pb <- rpbk[n_cells >= 5]
        den <- pb[, .(n_donors = uniqueN(donor_id), n_profiles = .N,
                      n_cells_total = sum(n_cells)), by = ccc_group][n_donors >= 5]
        s12a <- merge(s12a, den, by = "ccc_group", all.x = TRUE)
        cat(sprintf("RAS cohort: %d profiles, %d donors, %d cell types\n",
                    sum(den$n_profiles),
                    uniqueN(pb[ccc_group %in% den$ccc_group, donor_id]), nrow(den)))
    }
    setorder(s12a, gene, within_gene_rank)
    write_part(s12a, "12A",
        "Cell-type expression landscape of the local lung renin-angiotensin machinery",
        supports = "Figure S10",
        sources = c("agt_axis/_m/stats_data/ras_celltype_profile.tsv",
                    "agt_axis/_h/ras_panel.tsv",
                    "agt_axis/_m/ras_pseudobulk_celltype.tsv.gz"),
        notes = paste("`detect` is the donor-level detection fraction and",
                      "`emmean` the depth-adjusted marginal expression with donor",
                      "and study random effects. Within-gene ranks are computed",
                      "here for all 22 populations; the stored",
                      "ras_top_celltypes.tsv keeps only the top 3."))
}

## ---- S12B: circuit-role classification ----------------------------------
circ <- read_src(AG("stats_data", "ras_circuit_completeness.tsv"))
if (!is.null(circ)) {
    s12b <- copy(circ)
    ## `n_steps_present` counts all FIVE step columns, so it cannot support the
    ## "no cell type carried more than one of the three requirements" claim.
    ## The three core requirements are substrate, an AngII-generating protease,
    ## and the AT1 receptor.
    s12b[, n_core_roles := as.integer(has_substrate) + as.integer(has_protease) +
             as.integer(has_receptor)]
    ## Which gene carried each positive role: the chymase step is max(CMA1, CTSG)
    ## and the source discarded the winner.
    if (!is.null(prof)) {
        wide <- dcast(prof[, .(ccc_group, gene, detect)], ccc_group ~ gene,
                      value.var = "detect")
        gg <- function(g) if (g %in% names(wide))
            wide[match(s12b$ccc_group, wide$ccc_group), get(g)] else NA_real_
        cma <- gg("CMA1"); ctsg <- gg("CTSG")
        s12b[, chymase_gene := fifelse(is.na(cma) | is.na(ctsg), NA_character_,
                                       fifelse(cma >= ctsg, "CMA1", "CTSG"))]
        s12b[, genes_substrate := fifelse(has_substrate, "AGT", "")]
        s12b[, genes_receptor  := fifelse(has_receptor, "AGTR1", "")]
        s12b[, genes_protease  := fifelse(
            has_protease,
            paste(na.omit(c(fifelse(ace_step >= 0.05, "ACE", NA_character_),
                            fifelse(chymase_step >= 0.05, chymase_gene, NA_character_),
                            fifelse(renin_step >= 0.05, "REN", NA_character_))),
                  collapse = "+"), ""), by = ccc_group]
    }
    setnames(s12b, c("substrate", "receptor_AT1"),
             c("AGT_detection_fraction", "AGTR1_detection_fraction"),
             skip_absent = TRUE)
    setorder(s12b, -n_core_roles, -AGTR1_detection_fraction)
    write_part(s12b, "12B",
        "Circuit-role classification: completeness of the angiotensin II axis per cell type",
        supports = "Figure S10",
        sources = "agt_axis/_m/stats_data/ras_circuit_completeness.tsv",
        notes = paste("Roles are called at a donor-level detection threshold of",
                      "0.05. `n_core_roles` counts the three requirements",
                      "(substrate, AngII-generating protease, AT1 receptor) and is",
                      "computed here: the stored `n_steps_present` counts all five",
                      "step columns and therefore cannot support the",
                      "'no cell type carried more than one requirement' claim.",
                      "`chymase_gene` records which of CMA1/CTSG won the",
                      "max() the source script discarded."))
}

## ---- S12C: AGT rank stability, coexpression, target overlap -------------
c_bits <- list()
for (nm in c("agt_ligand_rank_bootstrap", "agt_ligand_coexpression",
             "agt_target_overlap")) {
    x <- read_src(AG("stats_data", paste0(nm, ".tsv")))
    if (!is.null(x)) c_bits[[nm]] <- x[, block := nm]
}
if (length(c_bits))
    write_part(rbindlist(c_bits, fill = TRUE), "12C",
        "AGT rank stability, ligand coexpression, and predicted-target overlap",
        supports = "Results (local RAS); cited in place of Figure S10",
        sources = c("agt_axis/_m/stats_data/agt_ligand_rank_bootstrap.tsv",
                    "agt_axis/_m/stats_data/agt_ligand_coexpression.tsv",
                    "agt_axis/_m/stats_data/agt_target_overlap.tsv"),
        notes = paste("CAVEAT on the coexpression block: rows with NA, and rows",
                      "with p_value == 0, arise where AGT is at or near zero",
                      "(DC2, monocyte and macrophage populations). AGT_SUMMARY.md",
                      "flags these as technical artifacts; they are retained for",
                      "completeness and must not be read as coexpression",
                      "evidence."))

## =========================================================================
## S13 -- disease associations
## =========================================================================
## Threshold labelling.
##
## 01.state_stats.R and 00.niche_index.py now stamp a `min_cells` column on every
## output and write the >=20 fit to `_mincells20` files. Outputs produced BEFORE
## that change carry no `min_cells` column and were computed at >=20 while
## occupying the unsuffixed (now "primary") filenames. Assuming the unsuffixed
## file is the >=10 fit would therefore relabel old >=20 results as the primary
## analysis -- exactly the silent-denominator problem this table exists to fix.
## So: trust the column, never the filename, and flag anything that lacks it.
PRE_HARMONIZATION <- FALSE
role_of <- function(x, sfx) {
    if ("min_cells" %in% names(x)) {
        mc <- unique(x$min_cells)[1]
        return(sprintf("%s (>=%s pericytes/donor)",
                       if (identical(as.integer(mc), 10L)) "PRIMARY" else "sensitivity", mc))
    }
    PRE_HARMONIZATION <<- TRUE
    "UNKNOWN THRESHOLD (output predates threshold harmonization; re-run upstream)"
}
stack_thresholds <- function(maker, label) {
    out <- list()
    for (sfx in c("", "_mincells20")) {
        x <- maker(sfx)
        if (is.null(x)) next
        x[, analysis_role := role_of(x, sfx)]
        out[[length(out) + 1L]] <- x
    }
    if (!length(out)) return(NULL)
    rbindlist(out, fill = TRUE)
}
st <- function() if (PRE_HARMONIZATION) "pending_upstream" else "complete"

## ---- S13A: cohort denominators ------------------------------------------
meta <- read_src(P("pericyte_states", "_m", "pericytes_states_metadata.tsv.gz"))
a_bits <- list()
if (!is.null(meta)) {
    map_dg <- function(lc) {
        lc <- as.character(lc)
        fcase(grepl("^Healthy", lc), "Healthy", lc == "COPD", "COPD",
              grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc,
                    ignore.case = TRUE), "Fibrotic_ILD", default = "Other")
    }
    meta[, age := suppressWarnings(as.numeric(age_or_mean_of_age_range))]
    meta[, disease_group := map_dg(lung_condition)]
    dn <- meta[, .(n_cells = .N, dataset = first(dataset), study = first(study),
                   disease_group = first(disease_group), sex = first(sex),
                   age = mean(age, na.rm = TRUE), smoking = first(smoking_status),
                   bmi = first(BMI)), by = donor_id]
    a_bits$den <- dn[, .(
        donors_in_atlas = .N,
        donors_ge10_pericytes = sum(n_cells >= 10),
        donors_ge20_pericytes = sum(n_cells >= 20),
        donors_ge10_complete_age_sex = sum(n_cells >= 10 & !is.na(age) & !is.na(sex)),
        donors_ge20_complete_age_sex = sum(n_cells >= 20 & !is.na(age) & !is.na(sex)),
        total_pericytes = sum(n_cells),
        pericytes_per_donor_median = as.numeric(median(n_cells)),
        pericytes_per_donor_min = min(n_cells),
        pericytes_per_donor_max = max(n_cells),
        n_studies = uniqueN(study), n_datasets = uniqueN(dataset),
        studies = paste(sort(unique(study)), collapse = "; "),
        n_missing_smoking = sum(is.na(smoking) | !nzchar(as.character(smoking)) |
                                    as.character(smoking) %in% c("nan", "NA", "unknown")),
        n_missing_bmi = sum(is.na(bmi)),
        n_missing_medication = .N   # HLCA carries no medication metadata at all
    ), by = disease_group][order(-donors_in_atlas)]
}
sm <- read_src(DA("smoking_availability_by_disease.tsv"))
if (!is.null(sm)) a_bits$smoking <- sm[, block := "smoking availability (disease_association)"]
if (length(a_bits))
    write_part(rbindlist(a_bits, fill = TRUE), "13A",
        "Analysis cohorts and denominators for the disease models",
        supports = "Figure 5; Figure S11",
        sources = c("pericyte_states/_m/pericytes_states_metadata.tsv.gz",
                    "disease_association/_m/mixed_model_forest/smoking_availability_by_disease.tsv"),
        notes = paste("Both donor filters are reported: >=10 pericytes is the",
                      "PRIMARY threshold and >=20 the sensitivity threshold.",
                      "`n_missing_medication` equals the group size because HLCA",
                      "carries no medication metadata at all (see",
                      "sensitivity/_m/stats_data/sensitivity_README.txt).",
                      "These analyses are underpowered and disease is partly",
                      "confounded with study: the null results below are NOT",
                      "evidence that composition is biologically invariant."))

## ---- S13B/C: composition models -----------------------------------------
label_adjusted <- function(dt) {
    if (is.null(dt) || !"p.value" %in% names(dt)) return(dt)
    setnames(dt, "p.value", "p_BH_within_level")
    dt[, adjustment := "BH within level (applied by emmeans::pairs)"][]
}

for (tag in c("state", "program")) {
    part <- if (tag == "state") "13B" else "13C"
    ttl <- if (tag == "state")
        "Six stable cluster fractions by disease group"
    else "Three dominant-program fractions by disease group"

    om <- stack_thresholds(function(sfx)
        read_src(PSS(paste0("composition_", tag, "_disease_anova_all", sfx, ".tsv"))), tag)
    lv <- if (!is.null(om)) unique(om$level) else character()

    emm_bits <- list(); ph_bits <- list()
    for (sfx in c("", "_mincells20")) {
        for (g in lv) {
            key <- gsub("[^A-Za-z0-9]+", "_", g)
            e <- read_src(PSS(paste0("composition_", tag, "_", key, "_emmeans", sfx, ".tsv")))
            p <- read_src(PSS(paste0("composition_", tag, "_", key, "_posthoc", sfx, ".tsv")))
            if (!is.null(e)) emm_bits[[paste(g, sfx)]] <-
                e[, `:=`(level = g, analysis_role = role_of(e, sfx),
                         block = "marginal means")]
            if (!is.null(p)) ph_bits[[paste(g, sfx)]] <-
                label_adjusted(p)[, `:=`(level = g, analysis_role = role_of(p, sfx),
                                         block = "pairwise contrasts")]
        }
    }
    all_bits <- c(list(if (!is.null(om)) om[, block := "omnibus disease test"]),
                  emm_bits, ph_bits)
    all_bits <- Filter(Negate(is.null), all_bits)
    if (length(all_bits))
        write_part(rbindlist(all_bits, fill = TRUE), part, ttl,
            supports = "Figure S11",
            sources = paste0("pericyte_states/_m/stats_data/composition_", tag, "_*"),
            status = st(),
            notes = paste("Donor-level ANCOVA frac ~ disease_group + age + sex.",
                          "`p_BH` in the omnibus block is BH across levels;",
                          "contrast P values are already BH-adjusted within level",
                          "by the source script and are NOT re-adjusted here."))
}

## ---- S13D: grouped injury-associated fraction ---------------------------
d_bits <- list()
for (sfx in c("", "_mincells20")) {
    e <- read_src(PSS(paste0("injury_fraction_emmeans", sfx, ".tsv")))
    p <- read_src(PSS(paste0("injury_fraction_posthoc", sfx, ".tsv")))
    if (!is.null(e)) d_bits[[paste0("e", sfx)]] <-
        e[, `:=`(analysis_role = role_of(e, sfx), block = "marginal means")]
    if (!is.null(p)) d_bits[[paste0("p", sfx)]] <-
        label_adjusted(p)[, `:=`(analysis_role = role_of(p, sfx),
                                 block = "pairwise contrasts")]
}
if (length(d_bits))
    write_part(rbindlist(d_bits, fill = TRUE), "13D",
        "Grouped injury-associated state fraction by disease group",
        supports = "Figure S11",
        sources = "pericyte_states/_m/stats_data/injury_fraction_*",
        status = st(),
        notes = paste("CORRECTED DEFINITION: the injury-associated group is the",
                      "activated/migratory program alone. Earlier summaries quoted",
                      "0.488 / 0.381 / 0.317 with P = 0.816, which were pre-relabel",
                      "values equal to basement-membrane + activated/migratory;",
                      "basement-membrane is a structural program and is",
                      "deliberately not counted as injury. The `injury_programs`",
                      "column records the definition actually used."))

## ---- S13E/F: primary and receptor-inclusive injury-stromal models -------
for (resp in c("injury_stromal_score", "injury_stromal_score_sens_agtr1")) {
    part <- if (grepl("sens", resp)) "13F" else "13E"
    ttl <- if (grepl("sens", resp))
        "Receptor-inclusive sensitivity composite (injury-stromal score + AGTR1+ fraction)"
    else "Primary injury-stromal disease model"
    bits <- list()
    for (sfx in c("", "_mincells20")) {
        for (kind in c("anova", "emmeans", "posthoc", "robust_coefs")) {
            x <- read_src(NI(paste0(resp, "_", kind, sfx, ".tsv")))
            if (is.null(x)) next
            ## car::Anova is written with row names, so fread lands the term
            ## labels in V1. Normalize to `term` rather than reconstructing them
            ## positionally, which would silently mislabel if the model changed.
            if ("V1" %in% names(x) && !"term" %in% names(x)) setnames(x, "V1", "term")
            bits[[paste(kind, sfx)]] <- x[, `:=`(block = kind,
                                                 analysis_role = role_of(x, sfx),
                                                 response = resp)]
        }
    }
    comp <- P("niche_index", "_m", "niche_index_components.txt")
    comp_txt <- if (file.exists(comp)) paste(readLines(comp), collapse = " | ") else NA
    if (length(bits))
        write_part(rbindlist(bits, fill = TRUE), part, ttl,
            supports = if (grepl("sens", resp)) "Figure S12" else "Figure 5A",
            sources = paste0("niche_index/_m/stats_data/", resp, "_*"),
            status = st(),
            notes = paste0("Model: lm(", resp, " ~ disease_group + age + sex). ",
                           "Components: ", comp_txt, ". ",
                           "The >=20 rows reproduce the previously published fit ",
                           "(F(2,31) = 5.97, P = 0.0064 for the primary score); ",
                           "the >=10 rows are the new primary analysis."))
}

## =========================================================================
## S14 -- robustness and limitations
## =========================================================================
loso <- read_src(SE("leave_one_study_out.tsv"))
if (!is.null(loso)) {
    ## "significant in 13 of 16" is quoted but the count was never stored.
    summ <- loso[, .(n_refits = .N, n_significant = sum(p < 0.05),
                     estimate_min = min(estimate), estimate_max = max(estimate),
                     n_donors_min = min(n), n_donors_max = max(n)), by = response]
    write_part(rbindlist(list(loso[, block := "per-refit"],
                              summ[, block := "summary"]), fill = TRUE), "14A",
        "Leave-one-study-out robustness of the donor-level disease effects",
        supports = "Figure S12",
        sources = "sensitivity/_m/stats_data/leave_one_study_out.tsv",
        notes = paste("Each study is dropped in turn and the model refit.",
                      "The summary rows give the significant-refit count, which",
                      "was quoted in the summaries but never stored."))
}

nib <- list()
for (resp in c("niche_index", "niche_index_sens_agtr1")) {
    for (kind in c("anova", "emmeans", "posthoc", "robust_coefs")) {
        for (sfx in c("", "_mincells20")) {
            x <- read_src(NI(paste0(resp, "_", kind, sfx, ".tsv")))
            if (is.null(x)) next
            if ("V1" %in% names(x) && !"term" %in% names(x)) setnames(x, "V1", "term")
            nib[[paste(resp, kind, sfx)]] <- x[, `:=`(
                response = resp, block = kind, analysis_role = role_of(x, sfx))]
        }
    }
}
if (length(nib))
    write_part(rbindlist(nib, fill = TRUE), "14B",
        "Net niche-stability index by disease group",
        supports = "Figure S12",
        sources = "niche_index/_m/stats_data/niche_index_*",
        status = st(),
        notes = paste("The net index is the stability composite minus the",
                      "injury-stromal composite, so it is not independent of",
                      "Table S13E and should not be treated as a second test."))

sb <- list()
for (nm in c("smoking_availability_by_disease", "smoking_main_effect_healthy",
             "covariate_robustness_emmeans")) {
    x <- read_src(SE(paste0(nm, ".tsv")))
    if (!is.null(x)) sb[[nm]] <- x[, block := nm]
}
strat <- SE("smoking_stratified_injury.tsv")
if (file.exists(strat) && file.size(strat) == 0)
    sb$note <- data.table(block = "smoking_stratified_injury",
                          note = paste("EMPTY BY CONSTRUCTION, not a pipeline",
                                       "failure: no fibrotic/ILD or other-disease",
                                       "donor has smoking status recorded, so a",
                                       "smoking-stratified disease contrast is",
                                       "inestimable. Documented in",
                                       "sensitivity/_m/stats_data/sensitivity_README.txt."))
if (length(sb))
    write_part(rbindlist(sb, fill = TRUE), "14C",
        "Smoking-data availability and covariate robustness",
        supports = "Figure S12",
        sources = "sensitivity/_m/stats_data/*.tsv",
        notes = paste("Smoking is recorded for a minority of healthy donors and",
                      "for no fibrotic/ILD donor, so it cannot be adjusted for in",
                      "the disease models. The comorbidity-adjusted model collapses",
                      "to 12 donors and is reported for completeness only."))

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()

if (PRE_HARMONIZATION)
    warning("Some S13/S14 sources predate the threshold harmonization and carry ",
            "no `min_cells` column. Those parts are marked pending_upstream and ",
            "their rows say UNKNOWN THRESHOLD. Re-run:\n",
            "  sbatch -D pericyte_states/_m pericyte_states/_h/step_1.sh\n",
            "  sbatch -D niche_index/_m niche_index/_h/step_0.sh\n",
            "  sbatch -D niche_index/_m niche_index/_h/step_1.sh\n",
            "  sbatch -D sensitivity/_m sensitivity/_h/step_0.sh", call. = FALSE)
