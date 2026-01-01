suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggplot2)
    library(rstatix)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    sce <- zellkonverter::readH5AD("stroma.hlca_full.dataset.h5ad",
                                   reader = "R", use_hdf5 = TRUE)
    return(sce)
}

drop_zero_count_cells <- function(sce) {
    cnt  <- assay(sce, "counts")
    keep <- Matrix::colSums(cnt) > 0
    return(sce[, keep])
}

ensure_logcounts <- function(sce) {
    if ("logcounts" %in% assayNames(sce)) { return(sce) }
    if (!"counts" %in% assayNames(sce)) { stop("No counts assay") }

    sce <- drop_zero_count_cells(sce)
    sce <- scuttle::logNormCounts(sce)
    return(sce)
}

first_non_na <- function(x) {
    x <- x[!is.na(x)]
    return(if(length(x) == 0) NA else x[[1]])
}

prepare_agtr1_donor_table <- function(
    sce, cell_type_key = "cell_type", donor_key = "patient",
    age_key = "age_or_mean_of_age_range", sex_key = "sex",
    disease_key = "disease", ethnicity_key = "self_reported_ethnicity",
    gene = "AGTR1", min_cells_per_donor_celltype = 20,
    min_donors_per_celltype = 3
) {
    stopifnot(gene %in% rownames(sce))
    cd <- as.data.frame(colData(sce))

                                        # Extract AGTR1 logcounts for all cells
    ##sce <- ensure_logcounts(sce)
    ag  <- as.numeric(assay(sce, "X")[gene, ])
    cd$AGTR1_logcounts <- ag

    cols_needed <- c(donor_key, cell_type_key, age_key, sex_key, 
                     disease_key, ethnicity_key, "AGTR1_logcounts")
    cols_present <- cols_needed[cols_needed %in% colnames(cd)]
    cd <- cd[, cols_present, drop = FALSE]

                                        # Rename to standard names for downstream code
    cd <- cd |>
      dplyr::rename(
          donor_id = all_of(donor_key),
          cell_type = all_of(cell_type_key),
          age = all_of(age_key),
          disease = all_of(disease_key)
      ) |>
      filter(age > 20)

    if (sex_key %in% cols_present)      cd <- cd |> dplyr::rename(sex = all_of(sex_key))
    if (disease_key %in% cols_present)  cd <- cd |> dplyr::rename(disease = all_of(disease_key))
    if (ethnicity_key %in% cols_present)cd <- cd |> dplyr::rename(ethnicity = all_of(ethnicity_key))

                                        # Per donor x cell_type summaries
    donor_celltype <- cd |>
        group_by(donor_id, cell_type) |>
        summarise(
            n_cells = n(), AGTR1_mean = mean(AGTR1_logcounts, na.rm = TRUE),
            frac_AGTR1_pos = mean(AGTR1_logcounts > 0, na.rm = TRUE),
            age = mean(age, na.rm = TRUE), disease = first_non_na(disease),
            sex = if ("sex" %in% colnames(cd)) first_non_na(sex) else NA,
            ethnicity=if ("ethnicity" %in% colnames(cd)) first_non_na(ethnicity) else NA,
            .groups = "drop"
        ) |>
        tidyr::drop_na(age, AGTR1_mean) |>
        filter(n_cells >= min_cells_per_donor_celltype) |>
        group_by(cell_type) |>
        filter(n() >= min_donors_per_celltype) |>
        ungroup()

    return(donor_celltype)
}

filter_non_cancer_diseases <- function(df) {
    cancer_patterns <- c("carcinoma", "adenocarcinoma", "cancer")
    return(filter(df, !grepl(paste(cancer_patterns, collapse="|"),
                             disease, ignore.case = TRUE)))
}

get_agtr1_enriched_celltypes <- function(
    donor_celltype, metric = c("mean_expr", "frac_expr"), top_n = 5) {
    metric <- match.arg(metric)
    summ <- donor_celltype |>
      group_by(cell_type) |>
      summarise(
          n_donors = n(), mean_expr = mean(AGTR1_mean, na.rm = TRUE),
          frac_expr = mean(frac_AGTR1_pos, na.rm = TRUE),
          .groups = "drop"
      ) |>
      arrange(desc(.data[[metric]])) |>
      filter(cell_type != "Unknown") |>
      slice_head(n = top_n)
    return(summ)
}

run_disease_agtr1_by_celltype <- function(donor_celltype, keep_celltypes){
    donor_celltype |>
        filter(cell_type %in% keep_celltypes) |>
        group_by(cell_type) |>
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
    donor_celltype, keep_celltypes, outdir, filename = "disease_vs_agtr1_by_celltype") {
    dfp <- donor_celltype |>
        filter(cell_type %in% keep_celltypes) |>
        mutate(cell_type = forcats::fct_reorder(cell_type, AGTR1_mean, .fun = median))

    bxp <- ggboxplot(dfp, x = "disease", y = "AGTR1_mean", fill = "disease",
                     palette = "jco", add = "jitter", facet.by = "cell_type",
                     scales = "free", add.params = list(alpha=0.5, size=1),
                     panel.labs.font=list(face='bold'), xlab = "",
                     ylab = "Normalized Expression (AGTR1)", legend = "none",
                     ncol=5, ggtheme = theme_pubr(base_size = 15, border=TRUE)) +
        rotate_x_text(angle = 45, hjust = 1) +
        geom_pwc(p.adjust.method = "fdr", label = "p.format",
                 method = "dunn_test", tip.length = 0)

    save_ggplots(file.path(outdir, filename), bxp, w = 15, h = 5)
}

disease_agtr1_analysis <- function(
    sce, outdir = "all_cells", cell_type_key = "cell_type",
    donor_key = "patient", age_key = "age_or_mean_of_age_range",
    gene = "AGTR1", top_n_celltypes = 5, disease_key = "disease",
    enrichment_metric = "mean_expr",
    min_cells_per_donor_celltype = 20, min_donors_per_celltype = 3
) {
    if (!gene %in% rownames(sce)) {
        print("Set rownames")
        rownames(sce) <- rowData(sce)[, "feature_name"]
    }

                                        # Build donor x cell_type table for AGTR1
    print("Building donor by cell type table")
    donor_celltype <- prepare_agtr1_donor_table(
        sce, cell_type_key = cell_type_key, donor_key = donor_key,
        age_key = age_key, gene = gene, disease_key = disease_key,
        min_cells_per_donor_celltype = min_cells_per_donor_celltype,
        min_donors_per_celltype = min_donors_per_celltype
    ) |> filter_non_cancer_diseases()

                                        # Select enriched cell types descriptively
    print("Calculating enriched cell types")
    enriched <- get_agtr1_enriched_celltypes(
        donor_celltype, metric = enrichment_metric, top_n = top_n_celltypes
    )
    keep_celltypes <- enriched$cell_type

                                        # Correlation results + FDR
    print("Perform correlation")
    stats <- run_disease_agtr1_by_celltype(donor_celltype, keep_celltypes)

                                        # Save outputs
    print("Write data")
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(donor_celltype,
                       file.path(outdir, "donor_metadata.tsv"), sep = "\t")
    data.table::fwrite(enriched, file.path(outdir, "agtr1_enriched_celltypes.tsv"),
                       sep = "\t")
    data.table::fwrite(stats, file.path(outdir, "disease_agtr1_kruskal_by_celltype.tsv"),
                       sep = "\t")

                                        # Plot
    plot_disease_agtr1(donor_celltype, keep_celltypes, outdir)
    return(list(donor_celltype = donor_celltype, enriched = enriched, stats = stats))
}

#### Main
print("Load data")
sce <- load_data()

res <- disease_agtr1_analysis(
    sce, outdir = "mean_expr", cell_type_key = "cell_type",
    enrichment_metric = "mean_expr", top_n_celltypes = 5,
    min_cells_per_donor_celltype = 10
)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
