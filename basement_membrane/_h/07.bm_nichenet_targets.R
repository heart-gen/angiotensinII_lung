## Which niche ligands carry regulatory potential toward BASEMENT-MEMBRANE targets?
##
## ADDITIVE AND NON-DESTRUCTIVE BY DESIGN. cell_communication/_h/02.nichenet.R:78
## builds its gene set as `intersect(unique(unlist(TARGET_PROGRAM)), background)`
## -- it UNIONS every category into one target set. Adding a `basement_membrane`
## category there would therefore change the published ligand ranking, i.e. the
## "pericytes receive fibroblast-derived TGF-beta" result, its AUPR values and
## figure_ccc_nichenet; and 02b.nichenet_specificity.R holds a second hardcoded
## copy of the program that would silently desynchronize.
##
## So this script runs the same priors, the same receiver and the same expression
## threshold, but with geneset_oi = the BM panel alone. The published result is
## untouched and the two rankings can be compared directly.
##
## Question answered: do the same TGF-beta ligands that drive the pericyte ECM
## program also rank top for basement-membrane deposition, or is BM driven by a
## different set of signals?
##
## A permutation null is included because a ligand can rank highly simply by being
## well connected in the prior network; without the null, a BM ranking is not
## interpretable.

suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(nichenetr)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
PRIORS   <- parse_arg("--priors", "../../cell_communication/_m/nichenet_priors")
LIANADIR <- parse_arg("--liana-dir", "../../cell_communication/_m")
FRAC_FILE <- parse_arg("--frac-file", "expressed_fraction_main.tsv.gz")
PANELS   <- parse_arg("--panels", "./bm_panel_genes.tsv")
FROZEN   <- parse_arg("--frozen-activities",
                      "../../cell_communication/_m/nichenet/ligand_activities_Pericytes.tsv")
OUTDIR   <- parse_arg("--outdir", "./nichenet_bm")
EXPR_THR <- as.numeric(parse_arg("--expr-thr", "0.10"))
RECEIVER <- parse_arg("--receiver", "Pericytes")
NPERM    <- as.integer(parse_arg("--nperm", "1000"))
SEED     <- as.integer(parse_arg("--seed", "13"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

ligand_target_matrix <- readRDS(file.path(PRIORS, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(PRIORS, "lr_network_human_21122021.rds")) |>
    distinct(from, to)

frac <- data.table::fread(file.path(LIANADIR, FRAC_FILE), data.table = FALSE)
rownames(frac) <- frac[[1]]; frac[[1]] <- NULL

panels <- data.table::fread(PANELS)
bm_genes <- unique(panels[panels$panel == "basement_membrane", ]$gene)

all_ligands <- unique(lr_network$from)
all_receptors <- unique(lr_network$to)
genes_expressed_in <- function(group, thr = EXPR_THR) {
    if (!group %in% colnames(frac)) return(character(0))
    rownames(frac)[frac[[group]] >= thr]
}

expressed_receiver <- genes_expressed_in(RECEIVER)
background <- intersect(expressed_receiver, rownames(ligand_target_matrix))
geneset_oi <- intersect(bm_genes, background)
message("BM targets expressed in ", RECEIVER, ": ", length(geneset_oi), "/",
        length(bm_genes), " -> ", paste(geneset_oi, collapse = ", "))
if (length(geneset_oi) < 5)
    stop("fewer than 5 BM genes in the receiver background; ranking would be unstable")

sender_groups <- setdiff(colnames(frac), RECEIVER)
expressed_senders <- unique(unlist(lapply(sender_groups, genes_expressed_in)))
expressed_ligands <- intersect(all_ligands, expressed_senders)
expressed_receptors <- intersect(all_receptors, expressed_receiver)
potential_ligands <- lr_network |>
    filter(from %in% expressed_ligands, to %in% expressed_receptors) |>
    pull(from) |> unique()
potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))
message("potential ligands: ", length(potential_ligands))

activities <- predict_ligand_activities(
    geneset = geneset_oi, background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands) |>
    arrange(desc(aupr_corrected)) |>
    mutate(rank = row_number(), receiver = RECEIVER, target_set = "basement_membrane")

## ---- permutation null: is a high BM rank more than prior-network connectivity? ----
## Random gene sets of the same size drawn from the receiver background.
## Parallel over the allocated cores: at NPERM = 10000 the serial loop does not fit
## in any reasonable walltime. L'Ecuyer-CMRG + mc.set.seed gives reproducible
## parallel streams, so the run is still deterministic given SEED and NCORES.
NCORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
message("Permutation null: ", NPERM, " random gene sets of size ",
        length(geneset_oi), " on ", NCORES, " cores")
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)
perm_list <- parallel::mclapply(seq_len(NPERM), function(i) {
    rnd <- sample(background, length(geneset_oi))
    a <- suppressMessages(predict_ligand_activities(
        geneset = rnd, background_expressed_genes = background,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands))
    setNames(a$aupr_corrected, a$test_ligand)[activities$test_ligand]
}, mc.cores = NCORES, mc.set.seed = TRUE)
failed <- !vapply(perm_list, is.numeric, logical(1))
if (any(failed))
    stop(sum(failed), " of ", NPERM, " permutations failed; first error: ",
         conditionMessage(attr(perm_list[[which(failed)[1]]], "condition")))
perm <- matrix(unlist(perm_list), nrow = nrow(activities))

obs <- activities$aupr_corrected
activities$perm_mean <- rowMeans(perm, na.rm = TRUE)
activities$perm_sd <- apply(perm, 1, sd, na.rm = TRUE)
activities$perm_z <- (obs - activities$perm_mean) / activities$perm_sd
activities$perm_p <- (rowSums(perm >= obs, na.rm = TRUE) + 1) / (NPERM + 1)
activities$perm_p_BH <- p.adjust(activities$perm_p, method = "BH")

write.table(activities,
            file.path(OUTDIR, paste0("ligand_activities_BM_", RECEIVER, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- comparison against the frozen (unmodified) ECM-program ranking ----------
if (file.exists(FROZEN)) {
    frozen <- data.table::fread(FROZEN, data.table = FALSE)
    cmp <- frozen |>
        select(test_ligand, aupr_frozen = aupr_corrected, rank_frozen = rank) |>
        inner_join(activities |>
                       select(test_ligand, aupr_bm = aupr_corrected,
                              rank_bm = rank, perm_p_BH),
                   by = "test_ligand")
    rho <- suppressWarnings(cor(cmp$rank_frozen, cmp$rank_bm, method = "spearman"))
    top20 <- length(intersect(head(frozen$test_ligand[order(frozen$rank)], 20),
                              head(activities$test_ligand[order(activities$rank)], 20)))
    tgf <- activities |> filter(test_ligand %in% c("TGFB1", "TGFB2", "TGFB3")) |>
        select(test_ligand, rank, aupr_corrected, perm_z, perm_p_BH)
    write.table(cmp, file.path(OUTDIR, "bm_vs_frozen_ligand_ranking.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    msg <- c(
        "BM-target ligand ranking vs the frozen ECM-program ranking",
        sprintf("Spearman rank correlation: %.3f", rho),
        sprintf("Top-20 overlap: %d/20", top20),
        "TGF-beta ligands in the BM ranking:",
        paste(utils::capture.output(print(tgf)), collapse = "\n"),
        "",
        "Top 15 ligands toward BM targets:",
        paste(utils::capture.output(print(
            head(activities[, c("test_ligand", "rank", "aupr_corrected",
                                "perm_z", "perm_p_BH")], 15))), collapse = "\n"))
    writeLines(msg, file.path(OUTDIR, "bm_nichenet_README.txt"))
    message(paste(msg, collapse = "\n"))
}

## ---- ligand -> BM target regulatory potential -------------------------------
best <- activities |> slice_max(aupr_corrected, n = 25) |> pull(test_ligand)
links <- best |>
    lapply(get_weighted_ligand_target_links, geneset = geneset_oi,
           ligand_target_matrix = ligand_target_matrix, n = 200) |>
    bind_rows() |> drop_na()
write.table(links, file.path(OUTDIR, paste0("ligand_target_links_BM_", RECEIVER, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(links)) {
    vis <- prepare_ligand_target_visualization(
        ligand_target_df = links, ligand_target_matrix = ligand_target_matrix,
        cutoff = 0.33)
    ord_l <- intersect(rev(best), colnames(vis))
    vis <- vis[rownames(vis), ord_l, drop = FALSE]
    p <- make_heatmap_ggplot(t(vis), "Prioritized ligands",
                             "Basement-membrane targets", color = "purple",
                             legend_position = "top", x_axis_position = "top",
                             legend_title = "Regulatory potential") +
        theme(axis.text.x = element_text(angle = 90, hjust = 0))
    ggsave(file.path(OUTDIR, "ligand_target_heatmap_BM.pdf"), p, width = 10, height = 7)
    ggsave(file.path(OUTDIR, "ligand_target_heatmap_BM.png"), p, width = 10, height = 7,
           dpi = 300)
}

sessioninfo::session_info()
