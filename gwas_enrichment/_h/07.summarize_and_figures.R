#!/usr/bin/env Rscript
## ===========================================================================
## 07.summarize_and_figures.R — collate MAGMA gene-property + S-LDSC partitioned
## heritability into one cell-state x trait enrichment summary, and plot:
##   (A) heatmap of MAGMA signed -log10(p) over annotations x traits
##   (B) forest of MAGMA enrichment (BETA +/- SE) for the 5 pericyte states across
##       lung-disease traits (the headline: injury-niche state enriched, with the
##       smoking and negative-control traits as specificity guards)
##   (C) MAGMA vs S-LDSC concordance for the focal annotations
## Reads MAGMA magma/geneprop/*.gsa.out and S-LDSC results/<key>/<trait>.results.
## ===========================================================================
.libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib", .libPaths()))
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})

M <- "../_m"
OUT <- file.path(M, "summary"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
STATES <- c("vascular_stabilizing", "inflammatory", "synthetic_contractile",
            "activated_migratory", "fibroblast_like")
LUNG_TRAITS <- c("COPD_UKB", "COPD_COPDGene", "IPF_UKB", "ASTHMA",
                 "FEV1", "FVC", "FEV1FVC")
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 300)
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h)
}

## ---- MAGMA gene-property -------------------------------------------------
magma_files <- list.files(file.path(M, "magma/geneprop"), "\\.gsa\\.out$", full.names = TRUE)
read_gsa <- function(f) {
    dt <- tryCatch(fread(f, skip = "VARIABLE"), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    base <- sub("\\.gsa\\.out$", "", basename(f))
    # trait = leading token before first "."; annotation may include the joint tag
    trait <- sub("\\..*$", "", base)
    dt[, trait := trait]
    dt[, model := ifelse(grepl("JOINT_states", base), "joint", "marginal")]
    dt[, .(trait, model, annotation = VARIABLE, NGENES, BETA, SE, P)]
}
magma <- rbindlist(lapply(magma_files, read_gsa), fill = TRUE)
if (nrow(magma)) {
    magma[, signed_logp := -log10(P) * sign(BETA)]
    magma[, fdr := p.adjust(P, "BH"), by = model]
    fwrite(magma, file.path(OUT, "magma_geneprop_summary.tsv"), sep = "\t")
}

## ---- S-LDSC partitioned heritability (focal annotation = last category) --
ldsc_files <- list.files(file.path(M, "results"), "\\.results$",
                         full.names = TRUE, recursive = TRUE)
read_ldsc <- function(f) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    row <- dt[.N]   # appended custom annotation is the last category
    z <- row[["Coefficient_z-score"]]
    data.table(annotation = basename(dirname(f)),
               trait = sub("\\.results$", "", basename(f)),
               enrichment = row[["Enrichment"]],
               enrichment_p = row[["Enrichment_p"]],
               coef_z = z, coef_p = pnorm(z, lower.tail = FALSE))
}
ldsc <- rbindlist(lapply(ldsc_files, read_ldsc), fill = TRUE)
if (nrow(ldsc)) {
    ldsc[, coef_fdr := p.adjust(coef_p, "BH")]
    fwrite(ldsc, file.path(OUT, "sldsc_summary.tsv"), sep = "\t")
}

## ---- (A) MAGMA heatmap ---------------------------------------------------
if (nrow(magma)) {
    h <- magma[model == "marginal"]
    h[, annotation := factor(annotation, levels = unique(annotation))]
    pHeat <- ggplot(h, aes(trait, annotation, fill = signed_logp)) +
        geom_tile(colour = "white") +
        geom_text(data = h[P < 0.05],
                  aes(label = ifelse(fdr < 0.05, "**", "*")), size = 3) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "signed\n-log10 P") +
        labs(x = NULL, y = NULL, title = "MAGMA gene-property enrichment") +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_fig("magma_enrichment_heatmap", pHeat, 8, 7)
}

## ---- (B) Pericyte-state forest across lung-disease traits ----------------
if (nrow(magma)) {
    fr <- magma[model == "marginal" & annotation %in% STATES &
                trait %in% LUNG_TRAITS]
    if (nrow(fr)) {
        fr[, annotation := factor(annotation, levels = rev(STATES))]
        pFor <- ggplot(fr, aes(BETA, annotation, colour = trait)) +
            geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
            geom_pointrange(aes(xmin = BETA - 1.96 * SE, xmax = BETA + 1.96 * SE),
                            position = position_dodge(width = 0.6), size = 0.3) +
            labs(x = "MAGMA enrichment beta (+/- 95% CI)", y = NULL,
                 title = "Pericyte-state heritability enrichment (lung traits)") +
            theme_bw(base_size = 9)
        save_fig("pericyte_state_forest", pFor, 7.5, 3.5)
    }
}

## ---- (C) MAGMA vs S-LDSC concordance -------------------------------------
if (nrow(magma) && nrow(ldsc)) {
    m <- magma[model == "marginal", .(annotation, trait, magma_logp = -log10(P))]
    l <- ldsc[, .(annotation, trait, ldsc_logp = -log10(coef_p))]
    cc <- merge(m, l, by = c("annotation", "trait"))
    if (nrow(cc) > 2) {
        rho <- suppressWarnings(cor(cc$magma_logp, cc$ldsc_logp,
                                    method = "spearman", use = "complete.obs"))
        pCon <- ggplot(cc, aes(magma_logp, ldsc_logp)) +
            geom_point(size = 0.8, alpha = 0.6) +
            geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, colour = "#B2182B") +
            labs(x = "MAGMA -log10 P", y = "S-LDSC -log10 P",
                 title = sprintf("Method concordance (Spearman rho=%.2f)", rho)) +
            theme_bw(base_size = 9)
        save_fig("magma_vs_sldsc_concordance", pCon, 4, 4)
        fwrite(cc, file.path(OUT, "method_concordance.tsv"), sep = "\t")
    }
}

cat("[done] summary + figures in", OUT, "\n")
sessionInfo()
