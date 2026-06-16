#!/usr/bin/env Rscript
## ===========================================================================
## 06.coloc_agtr1.R — colocalization at the AGTR1 locus (chr3, GRCh38) between
## COPD GWAS and lung cis-eQTL.  Tests whether the COPD risk signal and the
## genetic control of AGTR1 expression share a causal variant (high PP4 ->
## AGTR1 expression mediates genetic COPD risk in lung).
##
## Primary:   COPDGene COPD GWAS  x  COPDGene lung cis-eQTL (Saferali)  via coloc.abf
##            (coloc.abf uses beta^2/varbeta, so COPDGene's unsigned |beta| is valid)
## Support:   GTEx v11 Lung — is AGTR1 a significant eGene? (arrow, if available)
##
## AGTR1 = ENSG00000144891 (chr3).  AGTR2 = ENSG00000180772 (chrX) is NOT testable
## here — autosomal/EUR genetics only — and is reported as such.
## ===========================================================================
.libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib", .libPaths()))
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})
if (!requireNamespace("coloc", quietly = TRUE))
    stop("Install coloc into .Rlib on a login node: ",
         "Rscript -e 'install.packages(\"coloc\", lib=\".Rlib\", repos=\"https://cloud.r-project.org\")'")

AGTR1_ENSG <- "ENSG00000144891"
GWAS  <- "../../inputs/gwas/_m/COPD_COPDGene.phs000179.pha004496.txt"
EQTL  <- "../../inputs/gwas/_m/COPDGene_lung_eQTL.Saferali_nhw376.txt"
GTEX_LUNG <- "/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11/GTEx_Analysis_v11_eQTL/Lung.v11.eQTLs.signif_pairs.parquet"
GTEX_EGENES <- "/ocean/projects/bio250020p/shared/resources/public-data/gtex_v11/GTEx_Analysis_v11_eQTL/Lung.v11.eGenes.txt.gz"
N_GWAS   <- 20066    # COPDGene COPD-Primary meta-analysis
N_EQTL   <- 376      # Saferali NHW lung eQTL
OUT <- "../_m/coloc"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 300)
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h)
}

## ---- GWAS at AGTR1 locus -------------------------------------------------
g <- fread(GWAS, skip = "ID\tAnalysis")
setnames(g, old = grep("beta", names(g), ignore.case = TRUE, value = TRUE)[1], new = "absbeta")
setnames(g, c("SNP ID", "Chr ID", "Chr Position", "P-value", "SE"),
         c("rsid", "chr", "pos", "P", "SE"), skip_absent = TRUE)
g <- g[chr == 3 & !is.na(absbeta) & !is.na(SE) & SE > 0]
## AGTR1 cis window (GRCh38 ~148.40-148.45 Mb) ± 500 kb
g <- g[pos > 147.9e6 & pos < 149.0e6]
g[, varbeta := SE^2]
g <- g[!duplicated(rsid) & rsid %like% "^rs"]
cat(sprintf("GWAS SNPs in AGTR1 window: %d\n", nrow(g)))

## ---- eQTL for AGTR1 ------------------------------------------------------
e <- fread(EQTL)
e <- e[gene == AGTR1_ENSG]
e[, se := beta / tstat]
e <- e[is.finite(se) & se != 0]
e[, varbeta := se^2]
setnames(e, "rsID", "rsid")
e <- e[!duplicated(rsid) & rsid %like% "^rs"]
cat(sprintf("eQTL (Saferali) cis-SNPs for AGTR1: %d\n", nrow(e)))

shared <- intersect(g$rsid, e$rsid)
cat(sprintf("Shared rsIDs for coloc: %d\n", length(shared)))

run <- function() {
    if (length(shared) < 20) {
        cat("Too few shared SNPs for reliable coloc; reporting NA.\n")
        return(data.table(test = "COPDGene_GWAS_x_COPDGene_lung_eQTL",
                          nsnps = length(shared), PP0 = NA, PP1 = NA, PP2 = NA,
                          PP3 = NA, PP4 = NA))
    }
    g2 <- g[rsid %in% shared][order(rsid)]
    e2 <- e[rsid %in% shared][order(rsid)]
    d1 <- list(beta = g2$absbeta, varbeta = g2$varbeta, snp = g2$rsid,
               type = "cc", N = N_GWAS, s = 0.5)
    d2 <- list(beta = e2$beta, varbeta = e2$varbeta, snp = e2$rsid,
               type = "quant", N = N_EQTL)
    res <- coloc::coloc.abf(d1, d2)
    s <- as.list(res$summary)
    data.table(test = "COPDGene_GWAS_x_COPDGene_lung_eQTL",
               nsnps = s$nsnps, PP0 = s$PP.H0.abf, PP1 = s$PP.H1.abf,
               PP2 = s$PP.H2.abf, PP3 = s$PP.H3.abf, PP4 = s$PP.H4.abf)
}
coloc_res <- run()
fwrite(coloc_res, file.path(OUT, "coloc_AGTR1_summary.tsv"), sep = "\t")
print(coloc_res)

## ---- GTEx v11 Lung eGene support (orthogonal) ----------------------------
gtex_note <- "GTEx v11 Lung eGene status: not evaluated"
if (requireNamespace("arrow", quietly = TRUE) && file.exists(GTEX_EGENES)) {
    eg <- tryCatch(fread(GTEX_EGENES), error = function(e) NULL)
    if (!is.null(eg)) {
        gcol <- grep("gene_id", names(eg), value = TRUE)[1]
        qcol <- grep("qval", names(eg), value = TRUE)[1]
        hit <- eg[get(gcol) %like% AGTR1_ENSG]
        if (nrow(hit))
            gtex_note <- sprintf("GTEx v11 Lung: AGTR1 eGene qval=%.3g (%s)",
                                 hit[[qcol]][1], ifelse(hit[[qcol]][1] < 0.05,
                                 "significant", "n.s."))
    }
}
cat(gtex_note, "\n")
writeLines(c(gtex_note,
             "AGTR2 (ENSG00000180772, chrX): NOT testable — autosomal/EUR genetics only."),
           file.path(OUT, "coloc_AGTR1_notes.txt"))

## ---- Figures: regional + PP4 --------------------------------------------
reg <- rbind(
    data.table(rsid = g$rsid, pos = g$pos, logp = -log10(g$P), track = "COPD GWAS"),
    data.table(rsid = e$rsid, pos = NA_real_,
               logp = -log10(e$pvalue), track = "Lung AGTR1 eQTL"),
    fill = TRUE)
pReg <- ggplot(reg[track == "COPD GWAS"], aes(pos / 1e6, logp)) +
    geom_point(size = 0.6, alpha = 0.5, colour = "#2166AC") +
    labs(x = "chr3 position (Mb, GRCh38)", y = expression(-log[10](P)),
         title = "COPD GWAS — AGTR1 locus") + theme_bw(base_size = 9)

pp <- data.table(h = factor(c("H0", "H1", "H2", "H3", "H4"),
                            levels = c("H0", "H1", "H2", "H3", "H4")),
                 pp = as.numeric(coloc_res[, .(PP0, PP1, PP2, PP3, PP4)]))
pPP <- ggplot(pp, aes(h, pp, fill = h == "H4")) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("grey70", "#B2182B"), guide = "none") +
    labs(x = "coloc hypothesis", y = "posterior prob.",
         title = sprintf("PP4 = %.2f (shared causal variant)", coloc_res$PP4)) +
    ylim(0, 1) + theme_bw(base_size = 9)

save_fig("coloc_AGTR1", (pReg | pPP) + plot_annotation(tag_levels = "A"), 7.5, 3.2)
cat("[done] coloc outputs in", OUT, "\n")
sessionInfo()
