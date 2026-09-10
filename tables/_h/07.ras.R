## Supplementary Tables S12 (local renin-angiotensin machinery) and S13B/C/D
## (discrete pericyte composition by disease group).
##
## SPLIT 2026-09-10. This script used to be `07.ras_disease.R` and also built
## S13A, S13E, S13F and all of S14 -- the disease-association parts. Those moved
## to heart-gen/lung-pericyte-analysis along with `disease_association/`,
## `sensitivity/` and `niche_index/_h/01`, and now live there as
## `tables/_h/07.disease.R`. The parts that remain are the ones whose sources
## stayed: S12 comes from `agt_axis/`, and S13B/C/D from `pericyte_states/` --
## they are the tabular form of Figure S11, which is also still built here.
##
## The gaps at S13A and S13E/F are deliberate, and so is the absence of S14.
##
## Threshold convention throughout S13: >=10 pericytes per donor is PRIMARY and
## lives in the unsuffixed source files; >=20 is the sensitivity analysis in the
## `_mincells20` files. Both are stacked with a `min_cells` column so the two
## denominators can never be conflated the way they were when the modules
## disagreed silently.
##
## Composition posthoc P values are ALREADY BH-adjusted by
## emmeans::pairs(adjust = "BH") in the source scripts and are relabelled, not
## re-adjusted.
##
## Outputs: tsv/tableS12A.tsv .. S12C, tableS13B.tsv .. S13D

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

## P2-35. Both `write_part` calls below used to HARDCODE the model into their
## `notes` string, and both strings outlived the models they described: they
## still named `lm(... ~ disease_group + age + sex)` after P1-1/P1-2 moved `age`
## into a labelled `_ageadj` arm and added `(1 | study)`, and after P1-10 refit
## the >=20 arms. S13E even shipped a `model` column contradicting its own note,
## while S13B/S13C had no model column at all -- so the wrong note was their
## ONLY record of provenance.
##
## Derive it from the rows instead. Every source now carries `model` (added to
## the composition outputs in `pericyte_states/_h/01.state_stats.R` for this
## purpose), so the sentence cannot drift from the fit again.
model_sentence <- function(x) {
    if (!"model" %in% names(x)) return("Model: NOT RECORDED by the source script.")
    m <- unique(stats::na.omit(x$model)); m <- m[nzchar(m)]
    if (!length(m)) return("Model: NOT RECORDED by the source script.")
    if (!"arm" %in% names(x)) return(paste0("Model: ", paste(m, collapse = "; "), "."))
    ## Name the arm each formula belongs to; an unlabelled arm is the primary.
    ## A part can legitimately carry more than one model per arm -- S13E ships the
    ## guarded lmer alongside an HC3 comparison fit -- so name the blocks too,
    ## otherwise the sentence reads "primary = X; primary = Y" with no way to tell
    ## which rows are which.
    keys <- c("arm", "model", if ("block" %in% names(x)) "block")
    a <- unique(x[!is.na(model) & nzchar(model), ..keys])
    a[, arm_lab := fifelse(is.na(arm) | !nzchar(arm), "primary", sub("^_", "", arm))]
    a <- if ("block" %in% names(a))
        a[, .(blocks = paste(sort(unique(block)), collapse = "/")), by = .(arm_lab, model)]
        else a[, .(blocks = NA_character_), by = .(arm_lab, model)]
    paste0("Model (derived from the `model` column, not asserted): ",
           paste(sprintf("%s%s = %s", a$arm_lab,
                         fifelse(is.na(a$blocks), "", paste0(" [", a$blocks, "]")),
                         a$model), collapse = "; "), ".")
}

AG  <- function(...) P("agt_axis", "_m", ...)
PSS <- function(f) P("pericyte_states", "_m", "stats_data", f)

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
    ## Strata defined by the gene being scored carry no emmean and take no rank
    ## slot -- see CIRCULAR_UNITS in agt_axis/_h/01.ras_landscape_stats.R. Rank
    ## them with the rest and frank() would hand the artifact a position.
    s12a[, within_gene_rank := NA_integer_]
    s12a[!is.na(emmean),
         within_gene_rank := frank(-emmean, ties.method = "min"), by = gene]
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
                      "here for all populations the model could score; the stored",
                      "ras_top_celltypes.tsv keeps only the top 3.",
                      "Rows with circular_by_construction = TRUE are strata",
                      "DEFINED by the gene in that row -- the AT2 groups are split",
                      "on AGTR2 detectability, so their AGTR2 detection is 1.000",
                      "and 0.000 by construction. They are excluded before the",
                      "within-dataset standardization (the det stratum's raw AGTR2",
                      "is 117x the next cell type and would otherwise set the",
                      "scale for the whole gene), carry emmean = NA, and take no",
                      "rank. Their `detect` and `expr_raw` remain descriptive."))
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
        notes = paste("The coexpression block now carries an EXPRESSION FLOOR",
                      "(P2-18, 2026-09-08): a sender is tested only if it detects",
                      "AGT in >= 1% of its cells, computed cell-weighted over the",
                      "donors entering the fit. Read `tested`, `agt_detect_group`",
                      "and `excluded_reason`. 35 of 110 rows are tested; the 74",
                      "excluded by the floor and the 1 non-estimable row are",
                      "retained for completeness, are NOT in the BH family, and",
                      "must not be read as coexpression evidence. Before the",
                      "floor, 35 of the 42 BH-significant rows came from cell",
                      "types detecting AGT in under 1% of cells and sorted to the",
                      "top of the file; the 17 rows that printed `p_value == 0`",
                      "(Spearman underflow) were all among them and are now",
                      "flagged `p_underflow` and clamped to the double epsilon.",
                      "`agt_coexpression_floor_sweep.tsv` reports what floors of",
                      "0.005 / 0.01 / 0.02 / 0.05 would each have admitted.",
                      "CAVEAT on the rank block: `rank_median` is the centre of",
                      "the m-gene subsample distribution, NOT a bias-corrected",
                      "version of `rank_point`; read it against `rank_size_ref`,",
                      "the size-matched benchmark, and report the interval.",
                      "CAVEAT on the overlap block: `hyper_p_MISCALIBRATED`",
                      "assumes each ligand draws targets uniformly from the",
                      "24-gene shortlist, which is badly violated -- six targets",
                      "are used by 26-28 of the 30 ligands. Read `pair_pctile`",
                      "(rank among all ligand pairs in the same table) and",
                      "`degree_p` (degree-preserving permutation). The overlap",
                      "block is derived entirely from the curated NicheNet prior",
                      "and is not independent evidence from the coexpression",
                      "block."))

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
            notes = paste(model_sentence(rbindlist(all_bits, fill = TRUE)),
                          "Donor-level fractions. Only the PRIMARY arm is",
                          "assembled here; the age-complete `_ageadj` sensitivity",
                          "sits beside it in the source directory and is not",
                          "included, because `+ age` was a cohort filter on this",
                          "endpoint (P1-2), not a covariate.",
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


cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()

if (PRE_HARMONIZATION)
    warning("Some S13 sources predate the threshold harmonization and carry ",
            "no `min_cells` column. Those parts are marked pending_upstream and ",
            "their rows say UNKNOWN THRESHOLD. Re-run:\n",
            "  sbatch -D pericyte_states/_m pericyte_states/_h/step_1.sh",
            call. = FALSE)
