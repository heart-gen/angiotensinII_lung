## NicheNet ligand->target modeling (canonical nichenetr).
##
## Connects predicted incoming ligands (expressed in niche senders, receptors
## expressed in the receiver) to a curated pericyte injury/activation target
## program. This is the mechanistic link that converts "pericytes exist" into
## "pericytes receive identifiable signals that drive ECM / activation /
## migration / inflammatory programs".
##
## Receiver focus: ALL pericytes ("Pericytes") in the main analysis, and the
## functional-state receivers ("Pericyte_<program>") when run on the state-scheme
## fraction table; AGTR2-detectable AT2 epithelium ("AT2_AGTR2det"). Receivers and
## the per-group expressed-fraction table are passed in so the same script serves
## the main and state-stratified schemes.
## Inputs: NicheNet v2 priors + expressed_fraction_<scheme>.tsv.gz (from 01).

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
PRIORS  <- parse_arg("--priors", "../_m/nichenet_priors")
LIANADIR <- parse_arg("--liana-dir", "../_m")
OUTDIR  <- parse_arg("--outdir", "../_m/nichenet")
EXPR_THR <- as.numeric(parse_arg("--expr-thr", "0.10"))
FRAC_FILE <- parse_arg("--frac-file", "expressed_fraction_per_group.tsv.gz")
RECEIVERS <- strsplit(parse_arg("--receivers", "Pericytes,AT2_AGTR2det"), ",")[[1]]
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## Curated pericyte injury/activation target program (human symbols).
TARGET_PROGRAM <- list(
    ECM_remodeling = c("COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","SPARC",
                       "LOX","POSTN","MMP2","MMP3","MMP9","MMP14","TIMP1","TIMP2",
                       "LUM","DCN","BGN","FBN1"),
    activation_contractile = c("ACTA2","TAGLN","MYH11","CNN1","MYL9","TPM1","TPM2"),
    migration = c("ADAMTS1","THBS1","SERPINE1","CCN2","CCN1","PDGFRB","ITGA5","ITGB1"),
    inflammatory = c("IL6","CXCL1","CXCL2","CXCL8","CXCL10","CCL2","ICAM1","VCAM1",
                     "NFKBIA","SOD2","CCL20"),
    proliferation = c("MKI67","PCNA","TOP2A","CCND1","CCNB1","BIRC5")
)

## ---- Load priors --------------------------------------------------------
ligand_target_matrix <- readRDS(file.path(PRIORS, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(PRIORS, "lr_network_human_21122021.rds")) |>
    distinct(from, to)
weighted_networks <- readRDS(file.path(PRIORS, "weighted_networks_nsga2r_final.rds"))

frac <- data.table::fread(file.path(LIANADIR, FRAC_FILE), data.table = FALSE)
rownames(frac) <- frac[[1]]; frac[[1]] <- NULL

all_ligands <- unique(lr_network$from)
all_receptors <- unique(lr_network$to)

genes_expressed_in <- function(group, thr = EXPR_THR) {
    if (!group %in% colnames(frac)) return(character(0))
    rownames(frac)[frac[[group]] >= thr]
}

run_receiver <- function(receiver) {
    message("=== NicheNet receiver: ", receiver, " ===")
    sender_groups <- setdiff(colnames(frac), receiver)

    expressed_receiver <- genes_expressed_in(receiver)
    expressed_senders <- unique(unlist(lapply(sender_groups, genes_expressed_in)))
    if (length(expressed_receiver) < 50) {
        message("  receiver has too few expressed genes; skipping"); return(NULL)
    }

    background <- intersect(expressed_receiver, rownames(ligand_target_matrix))

    geneset_oi <- intersect(unique(unlist(TARGET_PROGRAM)), background)
    message("  geneset_oi (expressed target-program genes): ", length(geneset_oi))

    expressed_ligands <- intersect(all_ligands, expressed_senders)
    expressed_receptors <- intersect(all_receptors, expressed_receiver)
    potential_ligands <- lr_network |>
        filter(from %in% expressed_ligands, to %in% expressed_receptors) |>
        pull(from) |> unique()
    potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))
    message("  potential ligands: ", length(potential_ligands))
    if (length(potential_ligands) < 3 || length(geneset_oi) < 5) {
        message("  insufficient ligands/targets; skipping"); return(NULL)
    }

    ligand_activities <- predict_ligand_activities(
        geneset = geneset_oi, background_expressed_genes = background,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands) |>
        arrange(desc(aupr_corrected)) |>
        mutate(rank = row_number(), receiver = receiver)
    write.table(ligand_activities, file.path(OUTDIR, paste0("ligand_activities_", receiver, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    best_ligands <- ligand_activities |> slice_max(aupr_corrected, n = 30) |> pull(test_ligand)

    active_links <- best_ligands |>
        lapply(get_weighted_ligand_target_links, geneset = geneset_oi,
               ligand_target_matrix = ligand_target_matrix, n = 200) |>
        bind_rows() |> drop_na()
    write.table(active_links, file.path(OUTDIR, paste0("ligand_target_links_", receiver, ".tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    ## Ligand -> target heatmap
    vis <- prepare_ligand_target_visualization(
        ligand_target_df = active_links, ligand_target_matrix = ligand_target_matrix,
        cutoff = 0.33)
    order_ligands <- intersect(rev(best_ligands), colnames(vis))
    order_targets <- rownames(vis)
    vis <- vis[order_targets, order_ligands, drop = FALSE]
    p <- make_heatmap_ggplot(t(vis), "Prioritized ligands", "Target genes",
                             color = "purple",
                             legend_position = "top", x_axis_position = "top",
                             legend_title = "Regulatory potential") +
        theme(axis.text.x = element_text(angle = 90, hjust = 0))
    ggsave(file.path(OUTDIR, paste0("ligand_target_heatmap_", receiver, ".pdf")),
           p, width = 12, height = 8)
    ggsave(file.path(OUTDIR, paste0("ligand_target_heatmap_", receiver, ".png")),
           p, width = 12, height = 8, dpi = 300)

    ## Top-ligand activity barplot
    pa <- ligand_activities |> slice_max(aupr_corrected, n = 25) |>
        mutate(test_ligand = reorder(test_ligand, aupr_corrected)) |>
        ggplot(aes(test_ligand, aupr_corrected)) +
        geom_col(fill = "#3B6FB6") + coord_flip() +
        labs(x = "Ligand", y = "AUPR (corrected)",
             title = paste("Ligand activity ->", receiver)) +
        theme_bw(base_size = 12)
    ggsave(file.path(OUTDIR, paste0("ligand_activity_", receiver, ".pdf")), pa, width = 5, height = 7)
    ggsave(file.path(OUTDIR, paste0("ligand_activity_", receiver, ".png")), pa, width = 5, height = 7, dpi = 300)

    ligand_activities
}

for (recv in RECEIVERS) {
    tryCatch(run_receiver(recv),
             error = function(e) message("Receiver ", recv, " failed: ", conditionMessage(e)))
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
