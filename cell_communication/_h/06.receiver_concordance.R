## Do the prioritized ligands depend on how the pericyte receiver is defined?
##
## NicheNet is run against three independent receiver definitions:
##   whole pericytes            _m/nichenet/ligand_activities_Pericytes.tsv
##   stable state programs      _m/nichenet_state/ligand_activities_Pericyte_<program>.tsv
##   CoGAPS dominant patterns   _m/nichenet_cogaps/ligand_activities_Pericyte_cg_<program>.tsv
##
## The stable-program and CoGAPS schemes are derived independently (supervised
## Leiden clustering vs unsupervised NMF), so agreement between them is a genuine
## robustness check rather than a re-description of the same partition.
##
## Outputs (../_m/receiver_concordance/):
##   ligand_rank_matrix.tsv        long: test_ligand x receiver -> rank, aupr_corrected
##   receiver_rank_correlation.tsv pairwise Spearman of ranks over shared ligands
##   topk_overlap.tsv              pairwise top-K overlap count + Jaccard
##
## Ranks are compared over the INTERSECTION of each pair's candidate ligands: the
## candidate set is receiver-specific (a ligand only enters if its receptor clears
## the expression threshold in that receiver), so a union-based comparison would
## score absent ligands as maximally discordant purely from receptor dropout.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

TOP_K <- 20

## The negative-control receiver: a different cell type, deliberately.
CONTROL_RECEIVER <- "AT2_AGTR2det"

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
MDIR <- file.path(ROOT, "cell_communication", "_m")
OUT  <- file.path(MDIR, "receiver_concordance")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- collect every ligand-activity ranking on disk ----------------------
sources <- c(
    file.path(MDIR, "nichenet"),
    file.path(MDIR, "nichenet_state"),
    file.path(MDIR, "nichenet_cogaps")
)
files <- unlist(lapply(sources, function(d)
    if (dir.exists(d)) list.files(d, "^ligand_activities_.*\\.tsv$", full.names = TRUE)
    else character(0)))
if (!length(files)) stop("no ligand_activities_*.tsv found under ", MDIR)

## Pericyte receivers are the comparison of interest. AT2_AGTR2det -- a different
## CELL TYPE, not a different receiver definition for the same cells -- was
## previously excluded for that reason. It is now retained and labelled as the
## NEGATIVE CONTROL (defect P2-3, 2026-09-07).
##
## WHY. Every pericyte-vs-pericyte pair scores rho >= 0.96, including two
## biologically opposed programs (contractile vs vascular-stabilizing, 0.962) and
## one CoGAPS pattern against ALL pericytes pooled (0.998). That is presented as
## robustness, but it is equally consistent with the ranking being driven by the
## shared expressed-gene background and a tissue-agnostic prior rather than by
## receiver-specific biology -- and the design as it stood could not separate the
## two readings.
##
## AT2 is the discriminating comparison precisely BECAUSE it is a different cell
## type. It is not offered as a receiver definition; it is the answer to "what
## does this statistic look like when the receivers really are different?"
##   - pericyte-vs-AT2 markedly lower than pericyte-vs-pericyte -> the invariance
##     is receiver-specific biology, and the robustness reading holds.
##   - pericyte-vs-AT2 also ~0.96 -> the rankings are largely receiver-INDEPENDENT
##     and no state-stratified signaling panel can be read as state-specific.
act <- bind_rows(lapply(files, function(f) {
    d <- data.table::fread(f, data.table = FALSE)
    d$scheme <- basename(dirname(f))
    d
})) |>
    filter(grepl("^Pericyte", receiver) | receiver == CONTROL_RECEIVER) |>
    mutate(scheme = recode(scheme,
                           nichenet        = "whole",
                           nichenet_state  = "stable_program",
                           nichenet_cogaps = "cogaps"),
           scheme = ifelse(receiver == CONTROL_RECEIVER,
                           "NEGATIVE CONTROL -- different cell type", scheme),
           lineage = ifelse(receiver == CONTROL_RECEIVER, "AT2", "pericyte")) |>
    select(test_ligand, receiver, scheme, lineage, aupr_corrected, rank)

receivers <- sort(unique(act$receiver))
has_control <- CONTROL_RECEIVER %in% receivers
message("receiver definitions: ", paste(receivers, collapse = ", "))
if (!has_control)
    warning("negative-control receiver '", CONTROL_RECEIVER, "' not found -- ",
            "the P2-3 control cannot be computed and the invariance result ",
            "remains uninterpretable. Run 02.nichenet.R for that receiver.")

data.table::fwrite(act |> arrange(receiver, rank), sep = "\t",
                   file.path(OUT, "ligand_rank_matrix.tsv"))

## ---- pairwise rank agreement -------------------------------------------
pairs_df <- as.data.frame(t(combn(receivers, 2)),
                          stringsAsFactors = FALSE) |>
    setNames(c("receiver_a", "receiver_b"))

by_rec <- split(act, act$receiver)

concordance <- bind_rows(lapply(seq_len(nrow(pairs_df)), function(i) {
    a <- by_rec[[pairs_df$receiver_a[i]]]
    b <- by_rec[[pairs_df$receiver_b[i]]]
    j <- inner_join(a |> select(test_ligand, rank_a = rank, aupr_a = aupr_corrected),
                    b |> select(test_ligand, rank_b = rank, aupr_b = aupr_corrected),
                    by = "test_ligand")
    if (nrow(j) < 10) return(NULL)
    ct <- suppressWarnings(cor.test(j$rank_a, j$rank_b, method = "spearman"))
    data.frame(receiver_a = pairs_df$receiver_a[i],
               receiver_b = pairs_df$receiver_b[i],
               n_shared_ligands = nrow(j),
               spearman_rho = unname(ct$estimate),
               p_value = ct$p.value,
               pearson_aupr = cor(j$aupr_a, j$aupr_b))
}))
concordance$p_BH <- p.adjust(concordance$p_value, method = "BH")

## Label every pair by what it tests, so a reader cannot take the control for a
## receiver definition or read the within-pericyte rows as if they were controlled.
lin <- setNames(act$lineage[!duplicated(act$receiver)],
                act$receiver[!duplicated(act$receiver)])
concordance$comparison <- ifelse(
    lin[concordance$receiver_a] == lin[concordance$receiver_b],
    "within-pericyte -- receiver DEFINITION differs, cells do not",
    "NEGATIVE CONTROL -- different cell type")
data.table::fwrite(concordance |> arrange(comparison, desc(spearman_rho)), sep = "\t",
                   file.path(OUT, "receiver_rank_correlation.tsv"))

## ---- pairwise top-K overlap --------------------------------------------
topk <- lapply(by_rec, function(d) d$test_ligand[order(d$rank)][seq_len(min(TOP_K, nrow(d)))])
overlap <- bind_rows(lapply(seq_len(nrow(pairs_df)), function(i) {
    a <- topk[[pairs_df$receiver_a[i]]]; b <- topk[[pairs_df$receiver_b[i]]]
    inter <- length(intersect(a, b))
    data.frame(receiver_a = pairs_df$receiver_a[i],
               receiver_b = pairs_df$receiver_b[i],
               k = TOP_K, n_shared = inter,
               jaccard = inter / length(union(a, b)),
               shared_ligands = paste(sort(intersect(a, b)), collapse = ","))
}))
data.table::fwrite(overlap |> arrange(desc(n_shared)), sep = "\t",
                   file.path(OUT, "topk_overlap.tsv"))

cat("\n== Pairwise rank concordance across receiver definitions ==\n")
print(concordance |> arrange(desc(spearman_rho)), row.names = FALSE)
cat("\n== Top-", TOP_K, " overlap ==\n", sep = "")
print(overlap |> select(-shared_ligands) |> arrange(desc(n_shared)), row.names = FALSE)

## ---- P2-3: the negative-receiver control, decided in code ---------------
## The verdict is derived at runtime and written to disk rather than left to a
## reader's judgement, so it cannot drift from the numbers it describes.
if (has_control) {
    w <- concordance$spearman_rho[concordance$comparison ==
        "within-pericyte -- receiver DEFINITION differs, cells do not"]
    b <- concordance$spearman_rho[concordance$comparison ==
        "NEGATIVE CONTROL -- different cell type"]
    ## TWO comparisons, because they answer different questions and disagree here.
    ##   tail gap  = min(within) - max(control). "Can any within-pericyte pair be
    ##               mistaken for a cross-cell-type pair?" Driven by the extremes.
    ##   central   = median(within) - median(control), with a Wilcoxon on the two
    ##               distributions. "Do the two kinds of comparison differ at all?"
    ## Reporting only the tail gap understates a real separation; reporting only
    ## the central gap hides that the distributions nearly touch.
    gap <- min(w) - max(b)
    central <- median(w) - median(b)
    wt <- suppressWarnings(wilcox.test(w, b))

    sep_central <- central > 0.10 && wt$p.value < 0.05
    sep_tail    <- gap > 0.10

    verdict <- if (sep_central && sep_tail)
        paste0("SEPARATES ON BOTH: within-pericyte rho (median ", round(median(w), 3),
               ") exceeds pericyte-vs-AT2 (median ", round(median(b), 3),
               ") centrally and at the tails (gap ", round(gap, 3),
               "). The ranking responds to receiver identity; the within-pericyte ",
               "invariance is robustness, not a method floor.")
    else if (sep_central)
        paste0("SEPARATES CENTRALLY, OVERLAPS AT THE TAILS: within-pericyte median ",
               round(median(w), 3), " vs pericyte-vs-AT2 median ", round(median(b), 3),
               " (difference ", round(central, 3), ", Wilcoxon p = ",
               signif(wt$p.value, 3), ") -- so the ranking DOES respond to receiver ",
               "identity and this is not a pure method floor. BUT the lowest ",
               "within-pericyte pair (", round(min(w), 3),
               ") is only ", round(gap, 3), " above the highest cross-cell-type pair (",
               round(max(b), 3), "), so the MOST discordant pericyte pairs are barely ",
               "more concordant than a different cell type. Read the aggregate ",
               "prioritization as receiver-responsive; do NOT read any INDIVIDUAL ",
               "state-stratified pair as resolved.")
    else if (gap > 0 || central > 0)
        paste0("MARGINAL: within-pericyte median ", round(median(w), 3),
               " vs control median ", round(median(b), 3), " (Wilcoxon p = ",
               signif(wt$p.value, 3), "), tail gap ", round(gap, 3),
               ". Too weak to license state-stratified claims.")
    else
        paste0("NO SEPARATION: pericyte-vs-AT2 rho reaches ", round(max(b), 3),
               ", at or above the within-pericyte range (min ", round(min(w), 3),
               "). The rankings are largely RECEIVER-INDEPENDENT -- ",
               "state-stratified signaling panels show near-identical rankings ",
               "relabelled, and must not be read as state-specific.")

    summary_df <- data.frame(
        n_within_pericyte_pairs = length(w),
        n_control_pairs = length(b),
        within_pericyte_rho_min = min(w),
        within_pericyte_rho_median = median(w),
        within_pericyte_rho_max = max(w),
        control_rho_min = min(b),
        control_rho_median = median(b),
        control_rho_max = max(b),
        separation_gap_tail = gap,
        separation_gap_central = central,
        wilcoxon_p = wt$p.value,
        n_within_below_0.96 = sum(w < 0.96),
        verdict = verdict)
    data.table::fwrite(summary_df, sep = "\t",
                       file.path(OUT, "negative_receiver_control.tsv"))

    cat("\n== P2-3 NEGATIVE-RECEIVER CONTROL ==\n")
    cat(sprintf("  within-pericyte  rho: min %.3f  median %.3f  max %.3f  (%d pairs)\n",
                min(w), median(w), max(w), length(w)))
    cat(sprintf("  pericyte vs AT2  rho: min %.3f  median %.3f  max %.3f  (%d pairs)\n",
                min(b), median(b), max(b), length(b)))
    cat(sprintf("  tail gap    (min within - max control): %.3f\n", gap))
    cat(sprintf("  central gap (median - median):          %.3f  (Wilcoxon p = %.3g)\n",
                central, wt$p.value))
    cat(sprintf("  within-pericyte pairs below rho 0.96:   %d of %d\n\n  %s\n",
                sum(w < 0.96), length(w), verdict))
}

cat("\nReproducibility information:\n"); sessioninfo::session_info()
