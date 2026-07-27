## Supplementary figure S9: does the predicted signaling depend on how the pericyte
## receiver is defined, and does the basement-membrane program have its own senders?
##
## The primary NicheNet result treats all pericytes as one receiver. Two objections:
## (i) the ranking might be an artifact of that pooling, and (ii) if pericyte
## programs are distinct, the ligands predicted to drive them should differ. Panels
## A-C address (i) by re-running the prioritization against independently derived
## receiver definitions -- supervised stable-program clusters and unsupervised CoGAPS
## patterns -- and panel D addresses (ii) with the basement-membrane-restricted run.
##
##   A  ligand rank across the nine receiver definitions
##   B  fibroblast/myeloid -> pericyte TGF-beta and ECM edges, both receiver schemes
##   C  pericyte-derived CoGAPS patterns (nP=8) projected across the niche
##   D  basement-membrane-restricted ligand ranking and its permutation significance
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

CCC <- P("cell_communication", "_m")
BM  <- P("basement_membrane", "_m", "nichenet_bm")
CG  <- P("pericyte_cogaps", "_m")

TOP_PER_RECEIVER <- 8   # union of each receiver's top-N defines the panel-A rows
RANK_CAP <- 40          # ranks beyond this are visually equivalent ("not prioritized")

PROG_NICE <- c(vascular_stabilizing = "Vascular-stabilizing",
               synthetic_contractile = "Synthetic/contractile",
               inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
               activated_migratory = "Activated/migratory",
               basement_membrane = "Basement-membrane")

## Receiver display names, ordered whole -> stable program -> CoGAPS pattern.
rec_label <- function(x) {
    x <- as.character(x)
    out <- ifelse(x == "Pericytes", "All pericytes",
           ifelse(grepl("^Pericyte_cg_", x),
                  paste0(PROG_NICE[sub("^Pericyte_cg_", "", x)], " (CoGAPS)"),
                  paste0(PROG_NICE[sub("^Pericyte_", "", x)], " (state)")))
    unname(out)
}
rec_scheme <- function(x) fifelse(x == "Pericytes", "All",
                          fifelse(grepl("^Pericyte_cg_", x), "CoGAPS", "State"))

## ===== Panel A -- ligand rank across receiver definitions =================
rank_mat <- fread(file.path(CCC, "receiver_concordance", "ligand_rank_matrix.tsv"))
conc     <- fread(file.path(CCC, "receiver_concordance", "receiver_rank_correlation.tsv"))

keep_ligands <- unique(rank_mat[rank <= TOP_PER_RECEIVER, test_ligand])
ra <- rank_mat[test_ligand %in% keep_ligands]
## A ligand absent from a receiver's candidate set (its receptor is not expressed
## there) is shown as missing rather than as a bad rank -- those are different claims.
ra <- dcast(ra, test_ligand ~ receiver, value.var = "rank") |>
    melt(id.vars = "test_ligand", variable.name = "receiver", value.name = "rank")
ra[, rank_capped := pmin(rank, RANK_CAP)]
ra[, scheme := rec_scheme(receiver)]
ra[, receiver_lab := rec_label(receiver)]

rec_order <- unique(ra[order(match(scheme, c("All", "State", "CoGAPS")), receiver_lab),
                       receiver_lab])
ra[, receiver_lab := factor(receiver_lab, levels = rec_order)]
## Ligands ordered by their rank under the primary (all-pericyte) definition.
lig_order <- ra[receiver == "Pericytes"][order(rank), test_ligand]
lig_order <- c(lig_order, setdiff(unique(ra$test_ligand), lig_order))
ra[, test_ligand := factor(test_ligand, levels = rev(lig_order))]

rho_min <- min(conc$spearman_rho)
pA <- ggplot(ra, aes(receiver_lab, test_ligand, fill = rank_capped)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(rank <= TOP_PER_RECEIVER, rank, "")),
              size = 1.7, colour = "white") +
    ## direction = 1 so rank 1 is DARK: the prioritized ligands must be the salient
    ## ones, and the white rank labels need a dark tile to sit on.
    scale_fill_viridis_c(option = "mako", direction = 1, na.value = "grey88",
                         name = "Rank", breaks = c(1, 10, 20, 30, RANK_CAP),
                         labels = c("1", "10", "20", "30", paste0(">=", RANK_CAP))) +
    labs(x = NULL, y = NULL,
         caption = sprintf(paste("all %d receiver pairs: Spearman rho >= %.2f;",
                                 "grey = ligand not in that receiver's candidate set"),
                           nrow(conc), rho_min)) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          panel.grid = element_blank(), panel.border = element_blank(),
          legend.position = "right", legend.key.width = unit(3, "mm"),
          legend.key.height = unit(6, "mm"), legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6),
          plot.caption = element_text(size = 5.8, colour = "grey30", hjust = 0))

## ===== Panel B -- fibroblast/myeloid -> pericyte TGF-beta & ECM edges =====
FIBRO <- c("Adventitial fibroblasts", "Alveolar fibroblasts", "Myofibroblasts",
           "Peribronchial fibroblasts", "Subpleural fibroblasts")
MYELOID <- c("Alveolar macrophages", "Interstitial macrophages", "Classical monocytes",
             "Non-classical monocytes", "DC2")
## MMP14 is deliberately absent: it is panel D's top basement-membrane ligand, but
## LIANA's ligand-receptor resource carries no MMP14 ligand edge at all (0 rows), so
## including it would show an empty column rather than a negative result. The two
## panels answer different questions -- D is NicheNet regulatory potential over a
## target set, B is receptor-mediated ligand-receptor scoring.
LIGS <- c("TGFB1", "TGFB2", "TGFB3", "COL1A1", "COL4A1", "FN1", "TIMP2")

liana <- rbindlist(list(
    fread(file.path(CCC, "liana_into_receivers_state.tsv.gz"))[, scheme := "State"],
    fread(file.path(CCC, "liana_into_receivers_cogaps.tsv.gz"))[, scheme := "CoGAPS"]
), fill = TRUE)
liana <- liana[source %chin% c(FIBRO, MYELOID) &
               grepl("^Pericyte", target) &
               ligand_complex %chin% LIGS]
## Collapse the disease strata: this panel is about which senders reach pericytes at
## all, not about a disease contrast (that is the disease figure's job).
## Faceting is by receiver SCHEME rather than by individual receiver: the nine
## receiver names are far too long to serve as an axis, and the claim this panel
## makes is that the same senders reach pericytes under either definition -- which
## receiver-scheme columns show directly. Per-receiver resolution is panel A's job.
## Within a scheme, edges are summarised across its pericyte receivers by the
## strongest (max) score, i.e. "does this sender reach pericytes at all".
edges <- liana[, .(lrscore = max(lrscore, na.rm = TRUE),
                   spec = min(specificity_rank, na.rm = TRUE)),
               by = .(source, ligand_complex, scheme)]
edges[, sender_class := fifelse(source %chin% FIBRO, "Fibroblast", "Myeloid")]
edges[, ligand_complex := factor(ligand_complex, levels = LIGS)]
edges[, scheme := factor(scheme, levels = c("State", "CoGAPS"),
                         labels = c("Stable-program receivers", "CoGAPS receivers"))]

pB <- ggplot(edges, aes(ligand_complex, source, colour = lrscore,
                        size = -log10(spec + 1e-6))) +
    geom_point() +
    facet_grid(sender_class ~ scheme, scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(option = "rocket", direction = -1, name = "LR score") +
    scale_size_continuous(range = c(0.4, 2.6), name = expression(-log[10]~specificity)) +
    labs(x = NULL, y = NULL) +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 5.5),
          strip.text.x = element_text(size = 6, face = "bold"),
          strip.text.y = element_text(size = 6, face = "bold"),
          panel.grid.major = element_line(linewidth = 0.15),
          legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5))

## ===== Panel C -- CoGAPS nP=8 patterns projected across the niche =========
## Ported unchanged from manuscript_mechanism_figure.R (formerly figureS_cogaps_transfer
## panel B, and before that CCC main panel D). One row per PATTERN, labelled
## "Program (Pn)"; columns ordered by the fibrillar/disease anchor pattern.
pp <- fread(file.path(CG, "projected_pattern_by_celltype_np8.tsv"))
ann <- fread(file.path(CCC, "cogaps_receiver_annotation_np8.tsv"))
pats <- intersect(ann$pattern, names(pp)); ann <- ann[match(pats, ann$pattern)]
mat <- as.matrix(pp[, ..pats]); rownames(mat) <- pp$cell_type
z <- scale(mat)
anchor <- ann$pattern[which.max(ann$rho_fibroblast_like)]
ct_order <- rownames(z)[order(z[, anchor], decreasing = TRUE)]
prog_pri <- c("basement_membrane", "synthetic_contractile", "inflammatory",
              "vascular_stabilizing", "activated_migratory", "fibroblast_like")
selfrho <- vapply(seq_len(nrow(ann)), function(i)
    ann[[paste0("rho_", ann$assigned_program[i])]][i], numeric(1))
ann <- ann[order(match(ann$assigned_program, prog_pri), ann$pattern != anchor, -selfrho)]
rlab <- sprintf("%s (%s)", PROG_NICE[ann$assigned_program], sub("Pattern_", "P", ann$pattern))
names(rlab) <- ann$pattern
long <- as.data.table(as.table(z)); setnames(long, c("cell_type", "pattern", "z"))
long[, row := factor(rlab[as.character(pattern)], levels = rev(rlab))]
long[, cell_type := factor(cell_type, levels = ct_order)]
long <- long[!is.na(row)]

pC <- ggplot(long, aes(cell_type, row, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "z (within\npattern)") +
    labs(x = NULL, y = NULL) + theme_ms() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
          axis.text.y = element_text(size = 6), panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "right", legend.key.size = unit(3, "mm"),
          legend.title = element_text(size = 6), legend.text = element_text(size = 5.5))

## ===== Panel D -- basement-membrane-restricted ranking ====================
bm <- fread(file.path(BM, "ligand_activities_BM_Pericytes.tsv"))[order(rank)]
n_cand <- nrow(bm)
bmt <- head(bm, 15)
bmt[, test_ligand := factor(test_ligand, levels = rev(test_ligand))]
bmt[, sig := fifelse(perm_p_BH < 0.05, "FDR < 0.05", "n.s.")]

cmpr <- fread(file.path(BM, "bm_vs_frozen_ligand_ranking.tsv"))
rho_bm <- cor(cmpr$rank_frozen, cmpr$rank_bm, method = "spearman")
shared20 <- length(intersect(cmpr[order(rank_frozen)][seq_len(20), test_ligand],
                             cmpr[order(rank_bm)][seq_len(20), test_ligand]))

pD <- ggplot(bmt, aes(aupr_corrected, test_ligand)) +
    geom_segment(aes(x = 0, xend = aupr_corrected, y = test_ligand, yend = test_ligand),
                 linewidth = 0.4, colour = "grey75") +
    geom_point(aes(colour = sig), size = 2) +
    scale_colour_manual(values = c(`FDR < 0.05` = "#D55E00", `n.s.` = "grey55"),
                        name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
    annotate("text", x = Inf, y = 1.5, hjust = 1.02, vjust = 0, size = 2.2,
             colour = "grey20",
             label = sprintf("vs injury-program ranking:\nSpearman rho = %.3f, %d/20 top ligands shared\n%d candidate ligands",
                             rho_bm, shared20, n_cand)) +
    labs(x = "Corrected AUPR (basement-membrane target set)", y = NULL) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6),
          legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6))

## ---- assemble ------------------------------------------------------------
fig <- (pA | pB) / (pC | pD) +
    plot_layout(heights = c(1.25, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_receiver_robustness", fig, 11.0, 9.0)

cat("Wrote figureS_receiver_robustness to", OUT, "\n")
cat(sprintf("  panel A: %d ligands x %d receivers; min pairwise rho = %.3f\n",
            length(keep_ligands), length(rec_order), rho_min))
cat(sprintf("  panel D: %d BM candidates, rho vs injury ranking = %.3f, %d/20 shared\n",
            n_cand, rho_bm, shared20))
cat("\nReproducibility information:\n"); sessioninfo::session_info()
