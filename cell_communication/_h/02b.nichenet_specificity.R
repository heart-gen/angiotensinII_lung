## NicheNet specificity control: gene-set permutation null.
##
## The main result (02.nichenet.R) prioritizes ligands (TGFB2/TGFB1, ...) by how
## well their regulatory-potential signature predicts a CURATED pericyte injury /
## activation target program. A skeptic's question: would *any* gene set of the
## same size score these ligands as highly? If so, the "ligand activity" is just
## a property of the ligand's connectivity, not of the injury program.
##
## Test: hold the receiver, background, and potential-ligand set fixed (exactly as
## in the main run for the primary receiver = all Pericytes). Recompute ligand
## activity for N_PERM random gene sets drawn from the receiver's expressed
## background, each the same size as the real geneset_oi. For each prioritized
## ligand compare its observed aupr_corrected against this null:
##   empirical p = (1 + #{null >= obs}) / (N_PERM + 1)
##   z          = (obs - mean_null) / sd_null
## A specific ligand->program link has obs far in the right tail (small p, large z).
##
## Output (to --outdir): nichenet_specificity_<receiver>.tsv  + a forest-style png.
suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(nichenetr); library(dplyr); library(tidyr); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
PRIORS   <- parse_arg("--priors", "../_m/nichenet_priors")
LIANADIR <- parse_arg("--liana-dir", "../_m")
OUTDIR   <- parse_arg("--outdir", "../_m/nichenet")
EXPR_THR <- as.numeric(parse_arg("--expr-thr", "0.10"))
FRAC_FILE <- parse_arg("--frac-file", "expressed_fraction_main.tsv.gz")
RECEIVER <- parse_arg("--receiver", "Pericytes")
N_PERM   <- as.integer(parse_arg("--n-perm", "1000"))
TOPN     <- as.integer(parse_arg("--top-ligands", "15"))
SEED     <- as.integer(parse_arg("--seed", "13"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

## curated target program -- IDENTICAL to 02.nichenet.R
TARGET_PROGRAM <- c(
    "COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","SPARC","LOX","POSTN",
    "MMP2","MMP3","MMP9","MMP14","TIMP1","TIMP2","LUM","DCN","BGN","FBN1",
    "ACTA2","TAGLN","MYH11","CNN1","MYL9","TPM1","TPM2",
    "ADAMTS1","THBS1","SERPINE1","CCN2","CCN1","PDGFRB","ITGA5","ITGB1",
    "IL6","CXCL1","CXCL2","CXCL8","CXCL10","CCL2","ICAM1","VCAM1","NFKBIA","SOD2","CCL20",
    "MKI67","PCNA","TOP2A","CCND1","CCNB1","BIRC5")

ligand_target_matrix <- readRDS(file.path(PRIORS, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(PRIORS, "lr_network_human_21122021.rds")) |> distinct(from, to)

frac <- data.table::fread(file.path(LIANADIR, FRAC_FILE), data.table = FALSE)
rownames(frac) <- frac[[1]]; frac[[1]] <- NULL
genes_expressed_in <- function(group, thr = EXPR_THR) {
    if (!group %in% colnames(frac)) return(character(0))
    rownames(frac)[frac[[group]] >= thr]
}

## --- reproduce the main run's setup for the receiver ------------------------
sender_groups <- setdiff(colnames(frac), RECEIVER)
expressed_receiver <- genes_expressed_in(RECEIVER)
expressed_senders  <- unique(unlist(lapply(sender_groups, genes_expressed_in)))
background <- intersect(expressed_receiver, rownames(ligand_target_matrix))
geneset_oi <- intersect(TARGET_PROGRAM, background)
expressed_ligands   <- intersect(unique(lr_network$from), expressed_senders)
expressed_receptors <- intersect(unique(lr_network$to), expressed_receiver)
potential_ligands <- lr_network |>
    filter(from %in% expressed_ligands, to %in% expressed_receptors) |>
    pull(from) |> unique() |> intersect(colnames(ligand_target_matrix))
message(sprintf("receiver=%s | background=%d | geneset_oi=%d | potential_ligands=%d | n_perm=%d",
                RECEIVER, length(background), length(geneset_oi),
                length(potential_ligands), N_PERM))
stopifnot(length(geneset_oi) >= 5, length(potential_ligands) >= 3)

activity_vec <- function(geneset, ligands) {
    predict_ligand_activities(
        geneset = geneset, background_expressed_genes = background,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = ligands) |>
        select(test_ligand, aupr_corrected)
}

## observed activity over the full prioritized set
obs_full <- activity_vec(geneset_oi, potential_ligands) |> rename(obs_aupr = aupr_corrected)

## A ligand's corrected-AUPR depends only on its own target signature vs the gene
## set (independent of which other ligands are scored), so the permutation null is
## run only over the top observed ligands -- identical result, ~10x faster.
null_ligands <- obs_full |> slice_max(obs_aupr, n = max(TOPN, 25)) |> pull(test_ligand)
obs <- obs_full |> filter(test_ligand %in% null_ligands)
message(sprintf("null over %d top ligands (of %d) x %d perms",
                length(null_ligands), length(potential_ligands), N_PERM))

## null: random background gene sets, same size as geneset_oi
k <- length(geneset_oi)
null_mat <- matrix(NA_real_, nrow = nrow(obs), ncol = N_PERM,
                   dimnames = list(obs$test_ligand, NULL))
for (b in seq_len(N_PERM)) {
    gs <- sample(background, k)
    av <- activity_vec(gs, null_ligands)
    null_mat[av$test_ligand, b] <- av$aupr_corrected
    if (b %% 100 == 0) message("  perm ", b, "/", N_PERM)
}

res <- obs |>
    rowwise() |>
    mutate(
        null_mean = mean(null_mat[test_ligand, ], na.rm = TRUE),
        null_sd   = sd(null_mat[test_ligand, ], na.rm = TRUE),
        z         = (obs_aupr - null_mean) / null_sd,
        p_emp     = (1 + sum(null_mat[test_ligand, ] >= obs_aupr, na.rm = TRUE)) / (N_PERM + 1)
    ) |>
    ungroup() |>
    arrange(desc(obs_aupr)) |>
    mutate(rank = row_number(), receiver = RECEIVER,
           p_emp_adj = p.adjust(p_emp, "BH"))

data.table::fwrite(res, file.path(OUTDIR, sprintf("nichenet_specificity_%s.tsv", RECEIVER)),
                   sep = "\t")
cat("\n== Top prioritized ligands vs random-geneset null ==\n")
print(head(res, TOPN))

## figure: top ligands, observed AUPR vs null mean +/- 2 sd
plt <- head(res, TOPN) |>
    mutate(test_ligand = reorder(test_ligand, obs_aupr))
p <- ggplot(plt, aes(y = test_ligand)) +
    geom_errorbarh(aes(xmin = null_mean - 2 * null_sd, xmax = null_mean + 2 * null_sd),
                   height = 0.3, colour = "grey60") +
    geom_point(aes(x = null_mean), colour = "grey50", size = 2, shape = 1) +
    geom_point(aes(x = obs_aupr, colour = p_emp_adj < 0.05), size = 3) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey40"),
                        name = "BH p < 0.05") +
    labs(x = "AUPR (corrected): observed vs random-geneset null (mean +/- 2 sd)",
         y = "", title = sprintf("NicheNet ligand-program specificity (%s)", RECEIVER)) +
    theme_bw(base_size = 12)
ggsave(file.path(OUTDIR, sprintf("nichenet_specificity_%s.pdf", RECEIVER)), p, width = 7, height = 5)
ggsave(file.path(OUTDIR, sprintf("nichenet_specificity_%s.png", RECEIVER)), p, width = 7, height = 5, dpi = 300)

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
