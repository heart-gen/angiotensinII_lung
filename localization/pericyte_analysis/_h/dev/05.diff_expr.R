#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(qvalue)
    library(ggplot2)
    library(data.table)
    library(SingleCellExperiment)
})

## Functions
source("../_h/registration_pseudobulk.R")

save_qc_plots <- function(fit, coef_name, res_out, outdir) {
    ## MA plot
    pdf(file.path(outdir, "MAplot_agtr1_detect.pdf"), width = 6, height = 5)
    plotMA(fit, coef = coef_name, main = "AGTR1+ vs AGTR1- (pericyte pseudobulk)")
    abline(h = 0, col = "grey50")
    dev.off()

    ## Volcano plot
    vol <- res_out; vol$neg_log10_p <- -log10(vol$P.Value + 1e-300)

    pdf(file.path(outdir, "volcanoPlot_agtr1_detect.pdf"), width = 6, height = 5)
    ggplot(vol, aes(x = logFC, y = neg_log10_p)) +
        geom_point(alpha = 0.5, size = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "AGTR1+ vs AGTR1- (pericyte pseudobulk)",
             x = "log2 fold-change (AGTR1+ vs AGTR1-)", y = "-log10(p-value)") +
        theme_bw()
    dev.off()
}

add_agtr1_expression <- function(sce, gene = "AGTR1") {
    rowData(sce)$gene_id <- rownames(sce)
    rownames(sce) <- rowData(sce)$feature_name
    x <- assay(sce, "X")[gene, , drop = TRUE]
    x <- as.numeric(x)
    sce[[paste0(gene, "_expr")]] <- x
    sce[[paste0(gene, "_detect")]] <- as.integer(x > 0)
    return(sce)
}

build_pseudobulk <- function(sce, opt, outdir) {
    pheno_data <- as.data.frame(colData(sce))
    keep_cells <- pheno_data[[opt$celltype_col]] == opt$pericyte_label &
    !is.na(pheno_data[[opt$agtr1_col]]) & !is.na(pheno_data[[opt$donor_col]])
    
    sce_sub   <- sce[, keep_cells]
    pheno_sub <- as.data.frame(colData(sce_sub))    
    pheno_sub$agtr1_class <- ifelse(pheno_sub[[opt$agtr1_col]] == 1,
                                    "AGTR1_pos", "AGTR1_neg")

                                        # Count cells per donor × class
    dt_cells  <- as.data.table(pheno_sub)
    dt_counts <- dt_cells[, .N, by = .(donor = get(opt$donor_col), agtr1_class)]

                                        # Keep donors with >= min_cells_per_class
    wide_counts <- dcast(dt_counts, donor ~ agtr1_class, value.var = "N", fill = 0)
    wide_counts[, keep_donor := (AGTR1_pos >= opt$min_cells_per_class &
                                 AGTR1_neg >= opt$min_cells_per_class)]
    kept_donors <- wide_counts[keep_donor == TRUE, donor]
    message("Donors kept (>= ", opt$min_cells_per_class, " cells per class): ",
            length(kept_donors))
    
                                        # Summary report of included/excluded donors
    summary_dt <- copy(wide_counts)
    summary_dt[, status := ifelse(keep_donor, "included", "excluded")]
    fwrite(summary_dt, file = file.path(outdir, "pseudobulk_donor_summary.tsv"),
           sep = "\t")

                                        # Restrict to kept donors
    keep_cells_pb <- pheno_sub[[opt$donor_col]] %in% kept_donors
    sce_pb        <- sce_sub[, keep_cells_pb]
    pheno_pb      <- as.data.frame(colData(sce_pb))
    pheno_pb$agtr1_class <- ifelse(pheno_pb[[opt$agtr1_col]] == 1,
                                   "AGTR1_pos", "AGTR1_neg")

                                        # Create sample_id = donor_agtr1class
    pheno_pb$sample_id <- paste(pheno_pb[[opt$donor_col]], pheno_pb$agtr1_class, sep = "_")
    dt_pb_meta <- as.data.table(pheno_pb)
    pb_samples <- dt_pb_meta[, .(donor = get(opt$donor_col), agtr1_class, sex = sex,
                                 age = age_or_mean_of_age_range, disease = disease, 
                                 n_cells = .N), by = sample_id] |> dplyr::distinct()

                                        # Aggregate counts: genes x sample_id
    counts_mat <- get_counts_matrix(sce_pb)
    pb_counts  <- t(rowsum(t(counts_mat), group = pheno_pb$sample_id))
    pb_counts  <- pb_counts[, pb_samples$sample_id, drop = FALSE]
    pb_samples$lib_size <- colSums(pb_counts) # Library sizes

                                        # Filter low expression
    bad <- which(colSums(pb_counts) == 0)
    pb_samples[bad, ]

    
                                        # Save pseudobulk counts + sample table
    fwrite(as.data.table(pb_counts, keep.rownames = "gene_id"),
           file = file.path(outdir, "pseudobulk_counts.tsv.gz"), sep = "\t")
    fwrite(pb_samples, file = file.path(outdir, "pseudobulk_sample_table.tsv"),
           sep = "\t")

                                        # Donor counts per class
    donors_agtr1_pos <- unique(pb_samples[agtr1_class == "AGTR1_pos"]$donor)
    donors_agtr1_neg <- unique(pb_samples[agtr1_class == "AGTR1_neg"]$donor)

    list(pb_counts = pb_counts, pb_samples = pb_samples,
         kept_donors = kept_donors, donors_agtr1_pos = donors_agtr1_pos,
         donors_agtr1_neg = donors_agtr1_neg)
}

run_de_analysis <- function(pb_counts, pb_samples, opt, outdir) {
    dge <- DGEList(counts = pb_counts)
    sample_meta <- pb_samples[
        match(colnames(pb_counts), pb_samples$sample_id),
        ]
    
    dge$samples <- cbind(dge$samples, sample_meta)

                                        # Filter low-expressed genes
    cpm_mat <- edgeR::cpm(dge)
    keep_genes <- rowSums(cpm_mat > 1) >= opt$min_cpm_samples
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

    message("Genes kept after CPM filtering: ", nrow(dge))

                                        # Recompute normalization
    dge <- calcNormFactors(dge)

                                        # Design matrix
    pb_samples$agtr1_class <- factor(
        pb_samples$agtr1_class, levels = c("AGTR1_neg", "AGTR1_pos")
    )

    design_data <- data.frame(
        agtr1_class = pb_samples$agtr1_class, sex = pb_samples[["sex"]],
        age = pb_samples[["age_or_mean_of_age_range"]],
        disease = pb_samples[["disease"]]
    )

    design <- model.matrix(
        ~ agtr1_class + sex + age + disease, data = design_data
    )

    colnames(design) <- make.names(colnames(design))

                                        # Initial voom for dup corr
    v0 <- voom(dge, design = design, plot = FALSE)
    corfit <- duplicateCorrelation(v0, design, block = pb_samples$donor)
    message("Estimated intra-donor correlation: ",
            round(corfit$consensus, 3))

                                        # Voom with correlation
    v <- voom(
        dge, design = design, plot = FALSE, block = pb_samples$donor,
        correlation = corfit$consensus
    )

    fit <- lmFit(
        v, design, block = pb_samples$donor, correlation = corfit$consensus
    )
    fit <- eBayes(fit, robust = TRUE)

                                        # Extract DE table
    coef_name <- "agtr1_classAGTR1_pos"
    if (!coef_name %in% colnames(coef(fit))) {
        stop("Could not find coefficient ", coef_name, " in fit.")
    }

    res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
    ## TODO: calculate SE
    res$gene_id <- rownames(res)

                                        # Storey q-values + pi0
    qobj <- qvalue(res$P.Value)
    res$qvalue <- qobj$qvalues
    pi0  <- qobj$pi0

    write.table(
        data.frame(pi0 = pi0), file = file.path(outdir, "de_qvalue_pi0.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    res_out <- res[, c("gene_id", "logFC", "AveExpr", "t", "P.Value",
                       "adj.P.Val", "qvalue", "B")] ## Add SE
    colnames(res_out)[colnames(res_out) == "adj.P.Val"] <- "FDR"

    fwrite(
        as.data.table(res_out),
        file = file.path(outdir, "de_agtr1_pos_vs_neg.tsv.gz"),
        sep = "\t"
    )

    list(fit = fit, coef_name = coef_name, res_out = res_out,
         pi0 = pi0, n_genes = nrow(dge))
}

write_summary_report <- function(opt, outdir, kept_donors, donors_agtr1_neg,
                                 donors_agtr1_pos, n_genes, pi0) {
    report_lines <- c(
        "# Pericyte AGTR1 pseudobulk DE summary", "",
        paste0("- Input SCE/H5AD: `", opt$sce, "`"),
        paste0("- Output dir: `", outdir, "`"), "",
        "## Pseudobulk construction",
        paste0("- Pericyte label: `", opt$pericyte_label, "` in `",
               opt$celltype_col, "`"),
        paste0("- AGTR1 column: `", opt$agtr1_col,
               "` → classes: AGTR1_neg / AGTR1_pos"),
        paste0("- Minimum cells per donor × class: ", opt$min_cells_per_class),
        "", paste0("- Donors included: ", length(kept_donors)),
        paste0("- Donors with AGTR1_neg pseudobulk: ", length(donors_agtr1_neg)),
        paste0("- Donors with AGTR1_pos pseudobulk: ", length(donors_agtr1_pos)),
        "", "## DE analysis", paste0("- Genes tested (post-CPM filtering): ", n_genes),
        paste0("- Design: ~ agtr1_class + sex + age + disease"),
        "- Random effect approximation: limma duplicateCorrelation (block = donor)",
        paste0("- π0 estimate (Storey): ", round(pi0, 3)), "")

    writeLines(report_lines, con = file.path(outdir, "pseudobulk_DE_summary.md"))
}

#### Main
                                        # Create output directory
outdir <- "diff_expr"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

                                        # Load data
fn    <- "../../celltype_subset/_m/pericyte.hlca_core.dataset.h5ad"
adata <- zellkonverter::readH5AD(fn)
sce   <- as(adata, "SingleCellExperiment")

                                        # Prepare data
sce <- add_agtr1_expression(sce, gene = "AGTR1")
pheno_data <- as.data.frame(colData(sce)) |>
    dplyr::mutate(agtr1_class = ifelse(AGTR1_detect == 1, "AGTR1_pos", "AGTR1_neg"))
sce$agtr1_class <- pheno_data$agtr1_class
sce$age    <- sce$age_or_mean_of_age_range

                                        # Count cells in each class
class_counts <- table(pheno_data$agtr1_class)
print(class_counts)

                                        # Generate pseudobulk
sce_pseudo <- registration_pseudobulk(
    sce, "subclusters", "donor_id", c("sex", "age"), min_ncells = 10
)

                                        # Build pseudobulk
pb <- build_pseudobulk(sce = sce, opt = opt, outdir = outdir)

                                        # Run DE analysis
de <- run_de_analysis(pb_counts = pb$pb_counts, pb_samples = pb$pb_samples,
                      opt = opt, outdir = outdir)

                                        # Generate QC plot
save_qc_plots(fit = de$fit, coef_name = de$coef_name, res_out = de$res_out,
              outdir = outdir)

                                        # Write report
write_summary_report(opt = opt, outdir = outdir, kept_donors = pb$kept_donors,
                     donors_agtr1_neg = pb$donors_agtr1_neg,
                     donors_agtr1_pos = pb$donors_agtr1_pos,
                     n_genes = de$n_genes, pi0 = de$pi0)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
