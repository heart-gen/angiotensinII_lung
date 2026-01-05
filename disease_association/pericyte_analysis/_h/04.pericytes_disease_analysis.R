suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggplot2)
    library(rstatix)
    library(data.table)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h) {
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(fn = "results/clustered_data.h5ad", pred_threshold = 0.75) {
    sce <- zellkonverter::readH5AD(fn)
    sce <- sce[, colData(sce)$prediction_confidence > pred_threshold]    
    return(sce)
}

first_non_na <- function(x) {
    x <- x[!is.na(x)]
    return(if(length(x) == 0) NA else x[[1]])
}

prepare_agtr1_donor_table <- function(
    sce, cell_type_key = "predicted_labels", donor_key = "patient",
    disease_key = "disease", gene = "AGTR1",
    min_cells_per_donor_celltype = 2, min_donors_per_celltype = 2
) {
    stopifnot(gene %in% rownames(sce))
    cd <- as.data.frame(colData(sce))

                                        # Extract AGTR1 logcounts for all cells
    assay_name <- ifelse("logcounts" %in% assayNames(sce),
                         "logcounts", "counts")
    ag <- as.numeric(assay(sce, assay_name)[gene, ])
    cd$AGTR1_logcounts <- ag

    cols_needed <- c(donor_key, cell_type_key, disease_key, "AGTR1_logcounts")
    cols_present <- cols_needed[cols_needed %in% colnames(cd)]
    cd <- cd[, cols_present, drop = FALSE]

                                        # Rename to standard names for downstream code
    cd <- cd |>
      dplyr::rename(
          donor_id = all_of(donor_key),
          subcluster = all_of(cell_type_key),
          disease = all_of(disease_key)
      )

    if (disease_key %in% cols_present) {
        cd <- cd |> dplyr::rename(disease = all_of(disease_key))
    }

                                        # Per donor x subcluster summaries
    donor_celltype <- cd |>
        group_by(donor_id, subcluster) |>
        summarise(n_cells = n(), AGTR1_mean = mean(AGTR1_logcounts, na.rm = TRUE),
                  frac_AGTR1_pos = mean(AGTR1_logcounts > 0, na.rm = TRUE),
                  disease = first_non_na(disease), .groups = "drop") |>
        tidyr::drop_na(AGTR1_mean) |>
        filter(n_cells >= min_cells_per_donor_celltype) |>
        group_by(subcluster) |>
        filter(n() >= min_donors_per_celltype) |>
        ungroup()

    return(donor_celltype)
}

get_agtr1_enriched_celltypes <- function(donor_celltype, metric = "mean_expr") {
    summ <- donor_celltype |>
      group_by(subcluster) |>
      summarise(
          n_donors = n(), mean_expr = mean(AGTR1_mean, na.rm = TRUE),
          frac_expr = mean(frac_AGTR1_pos, na.rm = TRUE),
          .groups = "drop"
      ) |>
      arrange(desc(.data[[metric]]))
    return(summ)
}

is_agtr1_positive <- function(x, threshold = 0) {
    x > threshold
}

quantify_pericytes_by_disease <- function(donor_celltype, agtr1_threshold = 0) {
    donor_celltype |>
        filter(!is.na(disease)) |>
        mutate(agtr1_positive = is_agtr1_positive(
                   AGTR1_mean, threshold = agtr1_threshold)) |>
        group_by(disease) |>
        summarise(
            n_donors = n_distinct(donor_id), n_pericyte_donors = n(),
            n_agtr1_pos = sum(agtr1_positive),
            frac_agtr1_pos = n_agtr1_pos / n_pericyte_donors,
            mean_agtr1 = mean(AGTR1_mean, na.rm = TRUE),
            .groups = "drop") |>
        arrange(desc(frac_agtr1_pos))
}

run_disease_agtr1_by_celltype <- function(donor_celltype){
    donor_celltype |>
        group_by(subcluster) |>
        group_modify(~{
            dat <- .x |> filter(!is.na(disease))
            n   <- nrow(dat)
            k   <- n_distinct(dat$disease)
            if (k < 2) return(tibble())
            kw  <- kruskal.test(AGTR1_mean ~ disease, data = dat)
            eta_sq_h <- (unname(kw$statistic) - k + 1) / (n - k)
            kw_tbl <- tibble(
                n_donors = n, n_diseases = k, kw_stat = unname(kw$statistic),
                kw_pval = kw$p.value, kw_eta_sq = eta_sq_h
            )
            dunn_tbl <- dat |>
                dunn_test(AGTR1_mean ~ disease, p.adjust.method = "fdr") |>
                mutate(effsize_r = statistic / sqrt(n),
                       effsize_mag = case_when(
                           abs(effsize_r) < 0.1 ~ "negligible",
                           abs(effsize_r) < 0.3 ~ "small",
                           abs(effsize_r) < 0.5 ~ "medium",
                           TRUE ~ "large"
                       ), test = "dunn")
            bind_cols(kw_tbl, dunn_tbl)
        }) |>
        ungroup() |>
        mutate(fdr = p.adjust(kw_pval, method = "fdr")) |>
        arrange(kw_pval, p.adj)
}

plot_disease_agtr1 <- function(
    donor_celltype, outdir, filename = "disease_vs_agtr1_by_celltype") {
    dfp <- donor_celltype |> 
        mutate(disease = factor(disease, levels = c("Control", "IPF", "COPD")))

    bxp <- ggboxplot(dfp, x = "disease", y = "AGTR1_mean", fill = "disease",
                     palette = "npg", add = "jitter", facet.by = "subcluster",
                     scales = "free_x", add.params = list(alpha=0.8), xlab = "",
                     ylab = "Normalized Expression (AGTR1)", legend="none",
                     outlier.shape=NA, ncol=5, panel.labs.font=list(face='bold'),
                     ggtheme = theme_pubr(base_size = 15, border=TRUE)) +
        rotate_x_text(angle = 45, hjust = 1) +
        geom_pwc(method = "dunn_test", p.adjust.method = "fdr",
                 label = "p.format")

    save_ggplots(file.path(outdir, filename), bxp, w = 12, h = 5)
}

disease_agtr1_analysis <- function(
    sce, outdir = "all_cells", cell_type_key = "subcluster",
    donor_key = "patient", gene = "AGTR1", disease_key ="disease",
    enrichment_metric = "mean_expr", min_cells_per_donor_celltype = 5,
    min_donors_per_celltype = 3
    ) {
    if (!gene %in% rownames(sce)) {
        rownames(sce) <- rowData(sce)[, "HGNC_EnsemblAlt_GeneID"]
    }

    if (is.null(colnames(sce))) {
        colnames(sce) <- colData(sce)[, "CellBarcode_Identity"]
    }

                                        # Build donor x subcluster table for AGTR1
    donor_celltype <- prepare_agtr1_donor_table(
        sce, cell_type_key = cell_type_key, donor_key = donor_key,
        gene = gene, disease_key = disease_key,
        min_cells_per_donor_celltype = min_cells_per_donor_celltype,
        min_donors_per_celltype = min_donors_per_celltype
    )

                                        # Enrichment of subclusters
    enriched <- get_agtr1_enriched_celltypes(
        donor_celltype, metric = enrichment_metric
    )
    peri_freq <- quantify_pericytes_by_disease(donor_celltype)

                                        # Correlation results + FDR
    stats <- run_disease_agtr1_by_celltype(donor_celltype)

                                        # Save outputs
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    fwrite(donor_celltype, file.path(outdir, "donor_metadata.tsv"), sep = "\t")
    fwrite(enriched, file.path(outdir, "agtr1_enriched_celltypes.tsv"), sep = "\t")
    fwrite(peri_freq, file.path(outdir, "agtr1_freq_perictyes.tsv"), sep = "\t")
    fwrite(stats, file.path(outdir, "disease_agtr1_kruskal_by_celltype.tsv"),
           sep = "\t")

                                        # Plot
    plot_disease_agtr1(donor_celltype, outdir)
    return(list(donor_celltype = donor_celltype, enriched = enriched,
                stats = stats, peri_freq = prei_freq))
}

#### Main
sce <- load_data()

res <- disease_agtr1_analysis(
    sce, outdir = "pericyte_subclusters", cell_type_key = "predicted_labels",
    enrichment_metric = "mean_expr", min_cells_per_donor_celltype = 2
)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
