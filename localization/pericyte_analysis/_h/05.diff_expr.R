#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(edgeR)
  library(limma)
  library(qvalue)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
    make_option(c("--sce"), type = "character", help = "Input AnnData .h5ad"),
    make_option(c("--outdir"), type = "character", default = "results", help = "Output directory"),
    make_option(c("--celltype_col"), type = "character", default = "subclusters",
                help = "colData column with cell type labels [default %default]"),
    make_option(c("--pericyte_label"), type = "character", default = "Pericytes",
                help = "Label used for pericytes in celltype_col [default %default]"),
    make_option(c("--agtr1_col"), type = "character", default = "AGTR1_detect",
                help = "colData column indicating AGTR1 detection (0/1) [default %default]"),
    make_option(c("--donor_col"), type = "character", default = "donor_id",
                help = "colData column for donor ID [default %default]"),
    make_option(c("--min_cells_per_class"), type = "integer", default = 30,
                help = "Minimum cells per donor × class to keep [default %default]"),
    make_option(c("--min_donors_per_class"), type = "integer", default = 5,
                help = "Minimum donors per class required [default %default]"),
    make_option(c("--min_cpm_samples"), type = "integer", default = 5,
                help = "Genes must have CPM > 1 in at least this many samples [default %default]")
)

## Functions
load_sce_object <- function(path) {
    ext <- tolower(tools::file_ext(path))
    message("Loading input from: ", path)
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        stop(
            "To read .h5ad files you must install the 'zellkonverter' package ",
            "(and its dependencies such as 'basilisk')."
        )
    }
    return(zellkonverter::readH5AD(path))
}

## TODO: auto detect raw counts (see cerebral project)
get_counts_matrix <- function(sce) {
    available <- SummarizedExperiment::assayNames(sce)
    if ("counts" %in% available) {
        return(as.matrix(SummarizedExperiment::assay(sce, "counts")))
    }
    if ("X" %in% available) {
        return(as.matrix(SummarizedExperiment::assay(sce, "X")))
    }
}

build_pseudobulk <- function(sce, opt, outdir) {
    cd <- as.data.frame(SummarizedExperiment::colData(sce))
    required_cols <- c(opt$celltype_col, opt$agtr1_col, opt$donor_col,
                       "sex", "age_or_mean_of_age_range", "disease")
    missing_cols <- setdiff(required_cols, colnames(cd))
    if (length(missing_cols) > 0) {
        stop("Missing required metadata columns: ", 
             paste(missing_cols, collapse = ", "))
    }
    keep_cells <- cd[[opt$celltype_col]] == opt$pericyte_label &
    !is.na(cd[[opt$agtr1_col]]) & !is.na(cd[[opt$donor_col]])
    
    sce_sub <- sce[, keep_cells]
    cd_sub <- as.data.frame(SummarizedExperiment::colData(sce_sub))
    
    if (ncol(sce_sub) == 0) {
        stop("No pericyte cells found after filtering.")
    }
    cd_sub$agtr1_class <- ifelse(cd_sub[[opt$agtr1_col]] == 1,
                                 "AGTR1_pos", "AGTR1_neg")

    # Count cells per donor × class
    dt_cells <- as.data.table(cd_sub)
    dt_counts <- dt_cells[, .N, by = .(donor = get(opt$donor_col), agtr1_class)]

    # Keep donors with >= min_cells_per_class in BOTH classes
    wide_counts <- dcast(dt_counts, donor ~ agtr1_class, value.var = "N", fill = 0)

    if (!all(c("AGTR1_pos", "AGTR1_neg") %in% colnames(wide_counts))) {
        wide_counts$AGTR1_pos <- wide_counts[["AGTR1_pos"]] %||% 0L
        wide_counts$AGTR1_neg <- wide_counts[["AGTR1_neg"]] %||% 0L
    }

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

    if (length(kept_donors) < opt$min_donors_per_class) {
        warning(
            "Fewer than ", opt$min_donors_per_class,
            " donors retained with >= ", opt$min_cells_per_class, " cells per class.")
    }

    # Restrict to kept donors
    keep_cells_pb <- cd_sub[[opt$donor_col]] %in% kept_donors
    sce_pb <- sce_sub[, keep_cells_pb]
    cd_pb <- as.data.frame(SummarizedExperiment::colData(sce_pb))
    cd_pb$agtr1_class <- ifelse(cd_pb[[opt$agtr1_col]] == 1,
                                "AGTR1_pos", "AGTR1_neg")

    # Create sample_id = donor_agtr1class
    cd_pb$sample_id <- paste(cd_pb[[opt$donor_col]], cd_pb$agtr1_class, sep = "_")

    # Sample table (one row per pseudobulk sample)
    dt_pb_meta <- as.data.table(cd_pb)
    pb_samples <- dt_pb_meta[, .(donor = get(opt$donor_col), agtr1_class, sex = "sex",
                                 age = "age_or_mean_of_age_range", disease = "disease", 
                                 n_cells = .N), by = sample_id]

    # Aggregate counts: genes x sample_id
    counts_mat <- get_counts_matrix(sce_pb, opt$assay)
    if (nrow(counts_mat) == 0L || ncol(counts_mat) == 0L) {
        stop("Counts matrix is empty after filtering.")
    }

    pb_counts <- t(rowsum(t(counts_mat), group = cd_pb$sample_id))
    pb_counts <- pb_counts[, pb_samples$sample_id, drop = FALSE] # ensure order
    pb_samples$lib_size <- colSums(pb_counts) # Library sizes
    
    # Save pseudobulk counts + sample table
    fwrite(as.data.table(pb_counts, keep.rownames = "gene_id"),
           file = file.path(outdir, "pseudobulk_counts.tsv.gz"), sep = "\t")

    fwrite(pb_samples, file = file.path(outdir, "pseudobulk_sample_table.tsv"),
           sep = "\t")

  ## Donor counts per class
  donors_agtr1_pos <- unique(pb_samples[agtr1_class == "AGTR1_pos"]$donor)
  donors_agtr1_neg <- unique(pb_samples[agtr1_class == "AGTR1_neg"]$donor)

  message("Donors with AGTR1_pos pseudobulk: ", length(donors_agtr1_pos))
  message("Donors with AGTR1_neg pseudobulk: ", length(donors_agtr1_neg))

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
    cpm_mat <- cpm(dge)
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

make_qc_plots <- function(fit, coef_name, res_out, outdir) {
    ## MA plot
    pdf(file.path(outdir, "MA_plot_agtr1_pos_vs_neg.pdf"), width = 6, height = 5)
    plotMA(fit, coef = coef_name, main = "AGTR1+ vs AGTR1- (pericyte pseudobulk)")
    abline(h = 0, col = "grey50")
    dev.off()

    png(file.path(outdir, "MA_plot_agtr1_pos_vs_neg.png"),
        width = 6, height = 5, units = "in", res = 300)
    plotMA(fit, coef = coef_name, main = "AGTR1+ vs AGTR1- (pericyte pseudobulk)")
    abline(h = 0, col = "grey50")
    dev.off()

    ## Volcano plot
    vol <- res_out
    vol$neg_log10_p <- -log10(vol$P.Value + 1e-300)

    pdf(file.path(outdir, "volcano_agtr1_pos_vs_neg.pdf"), width = 6, height = 5)
    ggplot(vol, aes(x = logFC, y = neg_log10_p)) +
        geom_point(alpha = 0.5, size = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "AGTR1+ vs AGTR1- (pericyte pseudobulk)",
             x = "log2 fold-change (AGTR1+ vs AGTR1-)", y = "-log10(p-value)") +
        theme_bw()
    dev.off()

    png(file.path(outdir, "volcano_agtr1_pos_vs_neg.png"), width = 6,
        height = 5, units = "in", res = 300)
    ggplot(vol, aes(x = logFC, y = neg_log10_p)) +
        geom_point(alpha = 0.5, size = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "AGTR1+ vs AGTR1- (pericyte pseudobulk)",
             x = "log2 fold-change (AGTR1+ vs AGTR1-)", y = "-log10(p-value)") +
        theme_bw()
    dev.off()
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
                                        # Parse inputs
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$sce)) { stop("Both --sce is required.") }

                                        # Create output directory
outdir <- opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

                                        # Load data
sce <- load_sce_object(opt$sce)

                                        # Build pseudobulk
pb <- build_pseudobulk(sce = sce, opt = opt, outdir = outdir)

                                        # Run DE analysis
de <- run_de_analysis(pb_counts = pb$pb_counts, pb_samples = pb$pb_samples,
                      opt = opt, outdir = outdir)

                                        # Generate QC plot
make_qc_plots(fit = de$fit, coef_name = de$coef_name, res_out = de$res_out,
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
