## projectR transfer: project the pericyte CoGAPS patterns onto the niche.
##
## Question: are the latent programs CoGAPS learned in pericytes -- the
## vascular-stabilizing pattern and the injury (inflammatory / fibroblast-like)
## patterns -- MIRRORED in the neighbouring niche cell types? A "yes" turns the
## correlative LIANA/NicheNet edge into a COORDINATED multicellular program:
## e.g. fibroblasts/myeloid loading the same injury pattern whose target genes
## the pericytes express.
##
## Method: projectR (Bioconductor) performs a linear projection of new expression
## onto a learned loadings (genes x patterns) matrix. We project cell-type x donor
## pseudobulk (04.niche_pseudobulk.py) so projected pattern usage is per-donor and disease-aware.
## projectR is preferred; if it is not installed (PSC compute nodes lack internet,
## so install it on a login node into ./.Rlib first) the script falls back to the
## identical ordinary-least-squares projection so the job still completes.
##
## Outputs (to --outdir):
##   - projected_patterns_per_sample_npN.tsv   sample-level projected weights + meta
##   - projected_pattern_by_celltype_npN.tsv   mean projected weight per cell_type
##   - projected_pattern_celltype_emmeans_npN.tsv  which cell types load each pattern
##   - injury_pattern_disease_npN.tsv          disease effect on injury-pattern usage
##   - figures/projected_pattern_heatmap_npN.{pdf,png}
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr); library(ggplot2)
    library(lme4); library(lmerTest); library(emmeans)
})
.libPaths(c("./.Rlib", .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
PB        <- parse_arg("--pseudobulk", "../_m/niche_pseudobulk_logmean.tsv.gz")
SAMPLES   <- parse_arg("--samples", "../_m/niche_pseudobulk_samples.tsv")
NP        <- as.integer(parse_arg("--npatterns", "5"))
LOADINGS  <- parse_arg("--loadings", sprintf("../_m/feature_loadings_np%d.tsv.gz", NP))
ANNOT     <- parse_arg("--annotation",
                       sprintf("../../cell_communication/_m/cogaps_receiver_annotation_np%d.tsv", NP))
OUTDIR    <- parse_arg("--outdir", "../_m")
dir.create(file.path(OUTDIR, "figures"), showWarnings = FALSE, recursive = TRUE)

save_gg <- function(fn, p, w, h)
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)

## ---- load loadings (genes x patterns) and pseudobulk (genes x samples) --------
load_df <- as.data.frame(fread(LOADINGS)); rownames(load_df) <- load_df[[1]]; load_df[[1]] <- NULL
L <- as.matrix(load_df)                                  # genes x patterns
pb_df <- as.data.frame(fread(PB)); rownames(pb_df) <- pb_df[[1]]; pb_df[[1]] <- NULL
Y <- as.matrix(pb_df)                                    # genes x samples

shared <- intersect(rownames(L), rownames(Y))
message(sprintf("Shared genes for projection: %d (loadings %d, pseudobulk %d)",
                length(shared), nrow(L), nrow(Y)))
L <- L[shared, , drop = FALSE]; Y <- Y[shared, , drop = FALSE]

## ---- project ------------------------------------------------------------------
proj <- NULL
if (requireNamespace("projectR", quietly = TRUE)) {
    message("Using projectR::projectR")
    proj <- projectR::projectR(data = Y, loadings = L, full = FALSE)  # patterns x samples
} else {
    message("projectR not installed -> OLS projection (install projectR on a login node for the canonical call)")
    ## (L'L)^-1 L'Y : ordinary least squares, patterns x samples (= projectR linear core)
    proj <- solve(crossprod(L), crossprod(L, Y))
}
proj <- as.data.frame(t(proj))                            # samples x patterns
proj$sample <- rownames(proj)

## ---- pattern -> curated program labels ---------------------------------------
ann <- tryCatch(fread(ANNOT), error = function(e) NULL)
prog_of <- setNames(colnames(L), colnames(L))
if (!is.null(ann) && all(c("pattern", "assigned_program") %in% names(ann))) {
    m <- ann$assigned_program[match(colnames(L), ann$pattern)]
    prog_of <- setNames(ifelse(is.na(m), colnames(L), paste0(colnames(L), " (", m, ")")),
                        colnames(L))
}

meta <- fread(SAMPLES)
df <- merge(proj, meta, by = "sample")
pat_cols <- colnames(L)
fwrite(df, file.path(OUTDIR, sprintf("projected_patterns_per_sample_np%d.tsv", NP)), sep = "\t")

## ---- (A) which cell types load each pattern (mean + donor-aware emmeans) ------
by_ct <- df |>
    group_by(cell_type) |>
    summarise(across(all_of(pat_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
fwrite(by_ct, file.path(OUTDIR, sprintf("projected_pattern_by_celltype_np%d.tsv", NP)), sep = "\t")

emm_rows <- list()
for (p in pat_cols) {
    d <- df; d$w <- d[[p]]
    d <- d |> mutate(cell_type = factor(cell_type), donor_id = factor(donor_id))
    fit <- tryCatch(
        suppressMessages(lmer(w ~ cell_type + (1 | donor_id), data = d)),
        error = function(e) lm(w ~ cell_type, data = d))
    e <- as.data.frame(emmeans(fit, ~ cell_type))
    e$pattern <- p; e$program <- prog_of[[p]]
    emm_rows[[p]] <- e
}
fwrite(rbindlist(emm_rows, fill = TRUE),
       file.path(OUTDIR, sprintf("projected_pattern_celltype_emmeans_np%d.tsv", NP)), sep = "\t")

## ---- (B) does injury-pattern usage rise with disease, per cell type? ----------
map_disease_group <- function(lc) {
    lc <- as.character(lc)
    case_when(grepl("^Healthy", lc) ~ "Healthy",
              lc %in% c("COPD") ~ "COPD",
              grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
              TRUE ~ "Other")
}
df <- df |> mutate(disease_group = factor(map_disease_group(disease_group)))
## injury patterns = those whose annotation maps to an injury program
injury_progs <- c("inflammatory", "fibroblast_like", "activated_migratory")
injury_pat <- if (!is.null(ann)) ann$pattern[ann$assigned_program %in% injury_progs] else character(0)
dx_rows <- list()
for (p in intersect(injury_pat, pat_cols)) {
    for (ct in unique(df$cell_type)) {
        d <- df |> filter(cell_type == ct) |> mutate(w = .data[[p]])
        if (n_distinct(d$disease_group) < 2 || nrow(d) < 8) next
        fit <- tryCatch(lm(w ~ disease_group, data = d), error = function(e) NULL)
        if (is.null(fit)) next
        co <- summary(fit)$coefficients
        for (term in rownames(co)[-1]) {
            dx_rows[[paste(p, ct, term)]] <- data.frame(
                pattern = p, program = prog_of[[p]], cell_type = ct, term = term,
                estimate = co[term, 1], se = co[term, 2], p = co[term, 4], n = nrow(d))
        }
    }
}
if (length(dx_rows))
    fwrite(rbindlist(dx_rows, fill = TRUE),
           file.path(OUTDIR, sprintf("injury_pattern_disease_np%d.tsv", NP)), sep = "\t")

## ---- (C) heatmap: cell_type x pattern (mean projected weight, column z-scored) -
hm <- by_ct |> tibble::column_to_rownames("cell_type") |> as.matrix()
hmz <- scale(hm)                                          # z within pattern (column)
colnames(hmz) <- prog_of[colnames(hmz)]
long <- as.data.frame(as.table(hmz)); names(long) <- c("cell_type", "pattern", "z")
p_hm <- ggplot(long, aes(pattern, cell_type, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(x = "CoGAPS pattern (pericyte-learned)", y = "Niche cell type",
         fill = "z (within\npattern)",
         title = sprintf("projectR transfer of pericyte patterns onto the niche (nP=%d)", NP)) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_gg(file.path(OUTDIR, "figures", sprintf("projected_pattern_heatmap_np%d", NP)), p_hm, 9, 7)

cat("\nReproducibility information:\n"); Sys.time(); sessioninfo::session_info()
