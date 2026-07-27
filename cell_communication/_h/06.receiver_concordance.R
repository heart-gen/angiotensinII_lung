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

## Only pericyte receivers belong in this comparison; the main run also emits
## AT2_AGTR2det, which is a different cell type, not a different receiver
## definition for the same cells.
act <- bind_rows(lapply(files, function(f) {
    d <- data.table::fread(f, data.table = FALSE)
    d$scheme <- basename(dirname(f))
    d
})) |>
    filter(grepl("^Pericyte", receiver)) |>
    mutate(scheme = recode(scheme,
                           nichenet        = "whole",
                           nichenet_state  = "stable_program",
                           nichenet_cogaps = "cogaps")) |>
    select(test_ligand, receiver, scheme, aupr_corrected, rank)

receivers <- sort(unique(act$receiver))
message("receiver definitions: ", paste(receivers, collapse = ", "))

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
data.table::fwrite(concordance |> arrange(desc(spearman_rho)), sep = "\t",
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

cat("\nReproducibility information:\n"); sessioninfo::session_info()
