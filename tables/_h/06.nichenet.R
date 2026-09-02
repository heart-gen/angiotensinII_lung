## Supplementary Tables S10 (NicheNet ligand prioritization for the pericyte
## INJURY program) and S11 (the separate BM-target-set analysis).
##
## The two are kept apart deliberately: they use different target sets and
## different permutation budgets (1,000 vs 10,000 matched null sets), so their
## empirical P floors differ (9.99e-4 vs 9.999e-5) and quoting one floor for the
## other would misstate the test.
##
## S10 Part C reports BOTH donor-level specifications, labelled by figure:
##   lm(+dataset fixed)  -> Figure 4D   (beta = 0.50)
##   lmer(+1|dataset)    -> Figure S8C  (beta = 0.649)
## These are not competing versions of one model; see
## cell_communication/_h/05.donor_validation_stats.R.
##
## Outputs: tsv/tableS10A.tsv .. S10D, tableS11A.tsv .. S11C

suppressPackageStartupMessages({
    library(data.table); library(dplyr)
})
source("../_h/_tab_common.R")

CC <- function(...) P("cell_communication", "_m", ...)
NB <- function(f) P("basement_membrane", "_m", "nichenet_bm", f)
NF <- function(f) P("basement_membrane", "_m", "nichenet_fibrillar", f)

## Ligands prioritized in the manuscript text. `gene_program_detection.tsv` covers
## only 34 curated panel genes and is missing CCN2, TIMP2, COPA, VTN, MMP14,
## COL5A1 and AGT, so the class annotation for these is curated here by hand.
TOP11 <- c("TGFB2", "TGFB1", "CCN2", "TIMP2", "COPA", "VTN",
           "COL4A1", "COL1A1", "MMP14", "COL5A1", "AGT")
LIGAND_CLASS <- c(
    TGFB1 = "TGF-beta", TGFB2 = "TGF-beta", TGFB3 = "TGF-beta",
    CCN2 = "matricellular", THBS1 = "matricellular", VTN = "matricellular",
    POSTN = "matricellular", SPP1 = "matricellular",
    COL1A1 = "collagen/ECM", COL1A2 = "collagen/ECM", COL3A1 = "collagen/ECM",
    COL4A1 = "collagen/ECM", COL5A1 = "collagen/ECM", COL6A1 = "collagen/ECM",
    FN1 = "collagen/ECM", LAMB1 = "collagen/ECM",
    TIMP1 = "protease/remodeler", TIMP2 = "protease/remodeler",
    MMP2 = "protease/remodeler", MMP9 = "protease/remodeler",
    MMP14 = "protease/remodeler", SERPINE1 = "protease/remodeler",
    ADAMTS1 = "protease/remodeler",
    AGT = "renin-angiotensin substrate",
    COPA = "vesicular transport/secreted")

## =========================================================================
## S10A -- ligand ranking
## =========================================================================
la <- read_src(CC("nichenet", "ligand_activities_Pericytes.tsv"))
ef <- read_src(CC("expressed_fraction_per_group.tsv.gz"))
li <- read_src(CC("liana_into_receivers_main.tsv.gz"))

EXPR_THR <- 0.10   # matches cell_communication/_h/02.nichenet.R:32
if (!is.null(la)) {
    s10a <- copy(la)

    ## Expressed sender populations at the threshold the NicheNet run itself used.
    if (!is.null(ef)) {
        setnames(ef, 1, "gene")
        efm <- melt(ef, id.vars = "gene", variable.name = "population",
                    value.name = "frac")
        send <- efm[frac >= EXPR_THR, .(
            n_expressing_populations = .N,
            expressed_sender_populations = paste(population, collapse = "; "),
            max_expression_fraction = max(frac)), by = gene]
        s10a <- merge(s10a, send, by.x = "test_ligand", by.y = "gene", all.x = TRUE)
        pf <- efm[population == "Pericytes", .(gene, pericyte_expression_fraction = frac)]
        s10a <- merge(s10a, pf, by.x = "test_ligand", by.y = "gene", all.x = TRUE)
    }

    ## Cognate receptors expressed in pericytes, from the LIANA edges into pericytes.
    if (!is.null(li)) {
        rec <- li[target == "Pericytes",
                  .(cognate_receptors_in_pericytes = paste(sort(unique(receptor_complex)),
                                                           collapse = "; "),
                    n_cognate_receptors = uniqueN(receptor_complex)),
                  by = .(ligand_complex)]
        s10a <- merge(s10a, rec, by.x = "test_ligand", by.y = "ligand_complex",
                      all.x = TRUE)
    }

    s10a[, prioritized_top11 := test_ligand %in% TOP11]
    s10a[, ligand_class := LIGAND_CLASS[test_ligand]]
    s10a[is.na(ligand_class), ligand_class := "other"]
    setorder(s10a, rank)
    write_part(s10a, "10A",
        "NicheNet ligand prioritization for the pericyte injury program",
        supports = "Figure 4",
        sources = c("cell_communication/_m/nichenet/ligand_activities_Pericytes.tsv",
                    "cell_communication/_m/expressed_fraction_per_group.tsv.gz",
                    "cell_communication/_m/liana_into_receivers_main.tsv.gz"),
        notes = paste0("All ", nrow(s10a), " ranked candidate ligands. ",
                       "`aupr_corrected` is the ranking statistic. A sender ",
                       "population counts as expressing a ligand at fraction >= ",
                       EXPR_THR, ", the same threshold 02.nichenet.R used. ",
                       "`ligand_class` is curated by hand: the repo's ",
                       "gene_program_detection.tsv covers only 34 panel genes and ",
                       "is missing seven of the eleven prioritized ligands."))
}

## =========================================================================
## S10B -- matched-gene-set specificity
## =========================================================================
sp <- read_src(CC("nichenet", "nichenet_specificity_Pericytes.tsv"))
if (!is.null(sp)) {
    N_PERM <- 1000L   # 02b.nichenet_specificity.R:35
    ## p_emp = (1 + #{null >= obs}) / (N_PERM + 1); invert for the raw count, which
    ## the source never wrote out.
    sp[, n_null_ge_obs := round(p_emp * (N_PERM + 1) - 1)]
    sp[, n_permutations := N_PERM]
    sp[, prioritized_top11 := test_ligand %in% TOP11]
    setnames(sp, c("z", "p_emp", "p_emp_adj"),
             c("z_distance_from_null", "p_empirical", "p_empirical_BH"))
    setorder(sp, rank)
    write_part(sp, "10B",
        "Matched-gene-set specificity of the prioritized ligands (injury target set)",
        supports = "Figure S8A-B",
        sources = "cell_communication/_m/nichenet/nichenet_specificity_Pericytes.tsv",
        notes = paste0("Top 25 ligands against ", N_PERM, " matched random target ",
                       "sets. The empirical P floor is 1/(", N_PERM, "+1) = ",
                       signif(1 / (N_PERM + 1), 4), "; every prioritized ligand ",
                       "sits at the floor, i.e. no null replicate reached the ",
                       "observed AUPR. Do NOT compare this floor with Table S11, ",
                       "which used 10,000 permutations."))
}

## =========================================================================
## S10C -- donor-level validation
## =========================================================================
dv <- read_src(CC("donor_validation", "donor_validation_results.tsv"))
s10c_status <- "complete"
if (!is.null(dv)) {
    if (!"figure" %in% names(dv) || !any(grepl("lm\\(\\+dataset", dv$model))) {
        s10c_status <- "pending_upstream"
        cat("NOTE: donor_validation_results.tsv predates the Fig 4D lm and the\n",
            "  TGFB2-alone model. Re-run cell_communication steps 04 and 05.\n")
    }
    if (!"adj_ci_lo" %in% names(dv) && all(c("adj_estimate", "adj_se") %in% names(dv)))
        dv[, `:=`(adj_ci_lo = adj_estimate - 1.96 * adj_se,
                  adj_ci_hi = adj_estimate + 1.96 * adj_se)]
    dv[, sender_definition := paste("Fibroblast/mural senders (alveolar, adventitial,",
                                    "peribronchial, subpleural fibroblasts;",
                                    "myofibroblasts; vascular smooth muscle),",
                                    ">=10 sender cells per donor")]
    dv[, receiver_definition := paste("All pericytes (NOT an AGTR1+ split, which",
                                      "would be dropout-driven), >=5 pericytes per donor")]
    dv[, p_formatted := fmt_p(adj_p)]
    write_part(dv, "10C",
        "Donor-level validation of the predicted sender-to-pericyte association",
        supports = "Figure 4D; Figure S8C",
        sources = "cell_communication/_m/donor_validation/donor_validation_results.tsv",
        status = s10c_status,
        notes = paste("TWO SPECIFICATIONS, both published, on different figures:",
                      "lm(+dataset fixed) is Figure 4D, where the panel is a",
                      "partial-regression plot whose drawn slope IS this",
                      "coefficient; lmer(+1|dataset) is Figure S8C. They give",
                      "materially different estimates and must not be quoted",
                      "interchangeably. BH is applied within each model family.",
                      if (s10c_status != "complete")
                          "PENDING: re-run cell_communication steps 04-05." else ""))
}

## =========================================================================
## S10D -- receiver-definition robustness
## =========================================================================
rc <- list()
for (f in c("ligand_rank_matrix", "receiver_rank_correlation", "topk_overlap")) {
    x <- read_src(CC("receiver_concordance", paste0(f, ".tsv")))
    if (!is.null(x)) rc[[f]] <- x[, block := f]
}
if (length(rc))
    write_part(rbindlist(rc, fill = TRUE), "10D",
        "Robustness of ligand prioritization to the definition of the receiver population",
        supports = "Figure S9A-C",
        sources = "cell_communication/_m/receiver_concordance/*.tsv",
        notes = paste("Receiver schemes: whole pericytes, the three stable",
                      "dominant-program receivers, and the five CoGAPS-defined",
                      "receivers. `topk_overlap` uses k = 20."))

## =========================================================================
## S11 -- BM-restricted ligand-target analysis
## =========================================================================
bm  <- read_src(NB("ligand_activities_BM_Pericytes.tsv"))
cmp <- read_src(NB("bm_vs_frozen_ligand_ranking.tsv"))
NPERM_BM <- 10000L   # basement_membrane/_h/step_4.sh:29

if (!is.null(bm)) {
    s11a <- copy(bm)
    s11a[, fdr_significant := perm_p_BH < 0.05]
    s11a[, n_permutations := NPERM_BM]
    if (!is.null(cmp)) {
        s11a <- merge(s11a, cmp[, .(test_ligand, aupr_frozen, rank_frozen)],
                      by = "test_ligand", all.x = TRUE)
        s11a[, rank_change_bm_minus_injury := rank - rank_frozen]
    }
    s11a[, ligand_class := LIGAND_CLASS[test_ligand]]
    s11a[is.na(ligand_class), ligand_class := "other"]
    setorder(s11a, rank)
    write_part(s11a, "11A",
        "Ligand prioritization against the basement-membrane target set",
        supports = "Figure S9D",
        sources = c("basement_membrane/_m/nichenet_bm/ligand_activities_BM_Pericytes.tsv",
                    "basement_membrane/_m/nichenet_bm/bm_vs_frozen_ligand_ranking.tsv"),
        notes = paste0("All ", nrow(s11a), " candidate ligands against the ",
                       "basement-membrane genes expressed in pericytes (12 of the ",
                       "20-gene panel; read the count from the run log, do not ",
                       "assume it), with ", format(NPERM_BM, big.mark = ","),
                       " matched random target sets -- ten times the injury-set ",
                       "analysis, so the empirical P floor is 1/(",
                       format(NPERM_BM, big.mark = ","), "+1) = ",
                       signif(1 / (NPERM_BM + 1), 5), ". `rank_frozen` is the ",
                       "injury-program rank, so rank_change is positive when a ",
                       "ligand matters less for BM than for injury."))
}

## ---- S11D: BM vs fibrillar-collagen target set --------------------------
## Two runs of the SAME script (basement_membrane/_h/07.bm_nichenet_targets.R) on
## identical priors, receiver and expression threshold; only the target gene set
## differs. That is what makes the two rankings comparable at all. The comparable
## quantity is the permutation z, not the raw corrected AUPR: AUPR still scales
## with target-set size and prior connectivity, and the two sets differ in both.
fibact <- read_src(NF("ligand_activities_FIB_Pericytes.tsv"))
if (!is.null(bm) && !is.null(fibact)) {
    keep <- c("test_ligand", "rank", "aupr_corrected", "perm_z", "perm_p",
              "perm_p_BH")
    a <- bm[, intersect(keep, names(bm)), with = FALSE]
    b <- fibact[, intersect(keep, names(fibact)), with = FALSE]
    s11d <- merge(a, b, by = "test_ligand", suffixes = c("_bm", "_fibrillar"))
    s11d[, z_difference_bm_minus_fibrillar := perm_z_bm - perm_z_fibrillar]
    s11d[, rank_difference_bm_minus_fibrillar := rank_bm - rank_fibrillar]
    setorder(s11d, rank_bm)
    ok <- s11d[is.finite(rank_bm) & is.finite(rank_fibrillar)]
    rho <- if (nrow(ok) > 2)
        suppressWarnings(cor(ok$rank_bm, ok$rank_fibrillar, method = "spearman")) else
        NA_real_
    write_part(s11d, "11D",
        "Ligand prioritization toward basement-membrane versus fibrillar-collagen targets",
        supports = "Figure S17C",
        sources = c("basement_membrane/_m/nichenet_bm/ligand_activities_BM_Pericytes.tsv",
                    "basement_membrane/_m/nichenet_fibrillar/ligand_activities_FIB_Pericytes.tsv"),
        notes = paste0("Both columns come from the same script, priors, receiver ",
                       "and expression threshold; only geneset_oi differs, so the ",
                       "two rankings are comparable. Compare `perm_z`, not ",
                       "`aupr_corrected`: each z is standardised against its own ",
                       "target set's null of equally sized random gene sets drawn ",
                       "from the same receiver background, whereas AUPR still ",
                       "scales with target-set size. Spearman rank correlation ",
                       "between the two rankings: ", signif(rho, 3), " over ",
                       nrow(ok), " ligands. The fibrillar target set is small, so ",
                       "its individual estimates are noisy -- read it as a ",
                       "comparison of ligand ORDERING, not of absolute activity."))
}

## ---- S11B: BM vs injury ranking comparison ------------------------------
## The rank correlation and top-20 overlap exist only as text in
## bm_nichenet_README.txt; recompute so the table stands alone.
if (!is.null(cmp)) {
    ok <- cmp[is.finite(rank_bm) & is.finite(rank_frozen)]
    ct <- suppressWarnings(cor.test(ok$rank_bm, ok$rank_frozen, method = "spearman"))
    top_bm  <- ok[order(rank_bm)][1:20, test_ligand]
    top_inj <- ok[order(rank_frozen)][1:20, test_ligand]
    shared  <- intersect(top_bm, top_inj)
    s11b <- data.table(
        metric = c("Spearman rank correlation (BM vs injury)",
                   "n ligands compared", "top-20 overlap (count)",
                   "top-20 Jaccard", "shared top-20 ligands",
                   "BM-specific top-20 ligands", "injury-specific top-20 ligands"),
        value = c(sprintf("%.3f (P = %s)", ct$estimate, format(ct$p.value, digits = 3)),
                  as.character(nrow(ok)), as.character(length(shared)),
                  sprintf("%.3f", length(shared) / length(union(top_bm, top_inj))),
                  paste(sort(shared), collapse = "; "),
                  paste(sort(setdiff(top_bm, top_inj)), collapse = "; "),
                  paste(sort(setdiff(top_inj, top_bm)), collapse = "; ")))
    write_part(s11b, "11B",
        "Comparison of basement-membrane and injury-program ligand rankings",
        supports = "Figure S9D",
        sources = "basement_membrane/_m/nichenet_bm/bm_vs_frozen_ligand_ranking.tsv",
        notes = paste("RECOMPUTED for this table: these summary statistics existed",
                      "only as prose in bm_nichenet_README.txt. A modest rank",
                      "correlation with partial top-20 overlap is the evidence for",
                      "a partly distinct predicted regulatory architecture."))
}

## ---- S11C: ligand -> BM target links ------------------------------------
BM13 <- c("COL4A1", "COL4A2", "COL18A1", "LAMA3", "LAMA4", "LAMA5",
          "LAMB1", "LAMB2", "LAMC1", "NID1", "NID2", "HSPG2", "AGRN")
TARGET_CLASS <- c(COL4A1 = "collagen-associated", COL4A2 = "collagen-associated",
                  COL18A1 = "collagen-associated",
                  LAMA3 = "laminin", LAMA4 = "laminin", LAMA5 = "laminin",
                  LAMB1 = "laminin", LAMB2 = "laminin", LAMC1 = "laminin",
                  NID1 = "linker/proteoglycan", NID2 = "linker/proteoglycan",
                  HSPG2 = "linker/proteoglycan", AGRN = "linker/proteoglycan")

## Full ligand x target matrix, written by 07.bm_nichenet_targets.R --links-only.
## The older thresholded file is joined on as a flag rather than replaced, so the
## heatmap's 23 displayed links stay identifiable inside the complete matrix.
lt_full <- read_src(NB("ligand_target_matrix_BM_Pericytes.tsv"))
lt_thr  <- read_src(NB("ligand_target_links_BM_Pericytes.tsv"))
if (!is.null(lt_full)) {
    lt_full[, target_class := TARGET_CLASS[target]]
    if (!is.null(lt_thr)) {
        lt_thr[, key := paste(ligand, target)]
        lt_full[, shown_in_figure_heatmap := paste(ligand, target) %in% lt_thr$key]
    }
    ## Rank each target within a ligand and each ligand within a target: the
    ## interesting question in a dense matrix is relative preference, not the
    ## absolute weight, which varies by orders of magnitude between ligands.
    lt_full[, rank_target_within_ligand := frank(-regulatory_potential,
                                                 ties.method = "min"), by = ligand]
    lt_full[, rank_ligand_within_target := frank(-regulatory_potential,
                                                 ties.method = "min"), by = target]
    setorder(lt_full, ligand_rank, -regulatory_potential)
    n_lig <- uniqueN(lt_full$ligand); n_tgt <- uniqueN(lt_full$target)
    absent <- setdiff(BM13, unique(lt_full$target))
    write_part(lt_full, "11C",
        "Full ligand-to-basement-membrane-target regulatory-potential matrix",
        supports = "Figure S9D",
        sources = c("basement_membrane/_m/nichenet_bm/ligand_target_matrix_BM_Pericytes.tsv",
                    "basement_membrane/_m/nichenet_bm/ligand_target_links_BM_Pericytes.tsv"),
        notes = paste0(
            "Complete ", n_lig, " ligands x ", n_tgt, " targets = ", nrow(lt_full),
            " pairs, sliced directly from the NicheNet prior regulatory-potential ",
            "matrix. IMPORTANT: the target set is the ", n_tgt, " basement-membrane ",
            "genes expressed in pericytes above the 0.10 threshold, not all 13 -- ",
            if (length(absent))
                paste0(paste(absent, collapse = ", "),
                       " fall below the receiver expression threshold and were ",
                       "never in the tested gene set. ") else "",
            "`shown_in_figure_heatmap` marks the ",
            if (!is.null(lt_thr)) nrow(lt_thr) else 0,
            " links the Figure S9D heatmap displays, which are the pairs surviving ",
            "get_weighted_ligand_target_links(n = 200) plus the 0.33 visualization ",
            "cutoff. The matrix is dense, so absolute weights are only meaningful ",
            "relative to other pairs; the two rank columns give that comparison."))
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
