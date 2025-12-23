suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggplot2)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    fn  <- here::here("inputs/hlca/_m/hlca_core.rds")
    sce <- as.SingleCellExperiment(readRDS(fn))
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_coarse_for_GWAS_and_modeling
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type
    return(sce)
}

get_rowData <- function(){
    fn  <- here::here("inputs/hlca/_m/hlca_core.h5ad")
    sce <- zellkonverter::readH5AD(fn)
    dfx <- rowData(sce); rm(sce)
    return(dfx)
}

add_rowdata <- function(sce){
    dfx          <- get_rowData()
    sgenes       <- intersect(rownames(dfx), rownames(sce))
    sce          <- sce[sgenes, ]
    rowData(sce) <- dfx[sgenes,]
    return(sce)
}

add_qc <- function(sce){
    is_mito <- grep("MT-", rowData(sce)$feature_name)
    sce     <- scuttle::addPerCellQCMetrics(sce,
                                            subsets=list(mito=is_mito))
    sce     <- scuttle::addPerFeatureQCMetrics(sce)
    return(sce)
}

filter_qc <- function(sce){
                                        # Mitochondria percentage
    print(summary(sce$subsets_mito_percent))
    qc_stats <- scuttle::perCellQCFilters(sce)
    print(colSums(as.matrix(qc_stats)))
                                        # Discard outliers
    sce      <- sce[, !qc_stats$discard]
    return(sce)
}

first_non_na <- function(x) {
    x <- x[!is.na(x)]
    return(if(length(x) == 0) NA else x[[1]])
}

prepare_age_agtr1_donor_table <- function(
    sce, cell_type_key = "cell_type", donor_key = "patient",
    age_key = "age_or_mean_of_age_range", sex_key = "sex",
    disease_key = "disease", ethnicity_key = "self_reported_ethnicity",
    gene = "AGTR1", min_cells_per_donor_celltype = 20, 
    min_donors_per_celltype = 3
) {
    stopifnot(gene %in% rownames(sce))
    cd <- as.data.frame(colData(sce))

                                        # Extract AGTR1 logcounts for all cells
    ag <- as.numeric(logcounts(sce)[gene, ])
    cd$AGTR1_logcounts <- ag

    cols_needed <- c(donor_key, cell_type_key, age_key, sex_key, 
                     disease_key, ethnicity_key, "AGTR1_logcounts")
    cols_present <- cols_needed[cols_needed %in% colnames(cd)]
    cd <- cd[, cols_present, drop = FALSE]

                                        # Rename to standard names for downstream code
    cd <- cd |>
      rename(
          donor_id = all_of(donor_key),
          cell_type = all_of(cell_type_key),
          age = all_of(age_key)
      )

    if (sex_key %in% cols_present)      cd <- cd |> rename(sex = all_of(sex_key))
    if (disease_key %in% cols_present)  cd <- cd |> rename(disease = all_of(disease_key))
    if (ethnicity_key %in% cols_present)cd <- cd |> rename(ethnicity = all_of(ethnicity_key))

                                        # Per donor x cell_type summaries
    donor_celltype <- cd |>
      group_by(donor_id, cell_type) |>
      summarise(
          n_cells = n(), AGTR1_mean = mean(AGTR1_logcounts, na.rm = TRUE),
          frac_AGTR1_pos = mean(AGTR1_logcounts > 0, na.rm = TRUE),
          age = mean(age, na.rm = TRUE),
          sex = if ("sex" %in% colnames(cd)) first_non_na(sex) else NA,
          disease = if ("disease" %in% colnames(cd)) first_non_na(disease) else NA,
          ethnicity = if ("ethnicity" %in% colnames(cd)) first_non_na(ethnicity) else NA,
        .groups = "drop"
      ) |>
      tidyr::drop_na(age, AGTR1_mean)

                                        # Cell to donor QC
    donor_celltype <- donor_celltype |>
      filter(n_cells >= min_cells_per_donor_celltype)

                                        # Require enough donors per cell type
    donor_celltype <- donor_celltype |>
      group_by(cell_type) |>
      filter(n() >= min_donors_per_celltype) |>
      ungroup()

    return(donor_celltype)
}

get_agtr1_enriched_celltypes <- function(
    donor_celltype, metric = c("mean_expr", "frac_expr"), top_n = 10) {
    metric <- match.arg(metric)
    summ <- donor_celltype |>
      group_by(cell_type) |>
      summarise(
          n_donors = n(), mean_expr = mean(AGTR1_mean, na.rm = TRUE),
          frac_expr = mean(frac_AGTR1_pos, na.rm = TRUE),
          .groups = "drop"
      ) |>
      arrange(desc(.data[[metric]])) |>
      slice_head(n = top_n)
    return(summ)
}

run_age_agtr1_by_celltype <- function(donor_celltype, keep_celltypes) {
    res <- donor_celltype |>
      filter(cell_type %in% keep_celltypes) |>
      group_by(cell_type) |>
      summarise(
          n_donors = n(),
          rho = suppressWarnings(cor(age, AGTR1_mean, method = "spearman")),
          pval = suppressWarnings(cor.test(age, AGTR1_mean,
                                           method = "spearman")$p.value),
          .groups = "drop"
      ) |>
      mutate(fdr = p.adjust(pval, method = "fdr")) |>
      arrange(fdr, pval)
    return(res)
}

plot_age_agtr1_facets <- function(
    donor_celltype, keep_celltypes, outdir, filename = "age_vs_agtr1_by_celltype") {
    dfp <- donor_celltype |>
      filter(cell_type %in% keep_celltypes) |>
      mutate(cell_type = forcats::fct_reorder(cell_type, AGTR1_mean, .fun = median))

    sca <- ggscatter(dfp, x = "age", y = "AGTR1_mean", add = "loess",
                     facet.by = "cell_type", scales = "free_y",
                     conf.int = TRUE, cor.coef = TRUE, alpha = 0.6,
                     size = 1.2, xlab = "Donor age (years)",
                     ylab = "Donor-mean AGTR1 (logcounts)", legend = "none",
                     ggtheme = theme_pubr(base_size = 15))

    save_ggplots(file.path(outdir, filename), sca, w = 12, h = 8)
}

age_agtr1_analysis <- function(
      sce, outdir = "all_cells", cell_type_key = "cell_type",
      donor_key = "patient", age_key = "age_or_mean_of_age_range",
      gene = "AGTR1", top_n_celltypes = 10, write_donor = FALSE,
      enrichment_metric = c("mean_expr", "frac_expr"),
      min_cells_per_donor_celltype = 20, min_donors_per_celltype = 3
) {
    enrichment_metric <- match.arg(enrichment_metric)

    sce <- add_rowdata(sce)
    sce <- add_qc(sce)
    sce <- filter_qc(sce)
    rownames(sce) <- rowData(sce)[, "feature_name"]

                                        # Build donor x cell_type table for AGTR1
    donor_celltype <- prepare_age_agtr1_donor_table(
        sce, cell_type_key = cell_type_key, donor_key = donor_key,
        age_key = age_key, gene = gene,
        min_cells_per_donor_celltype = min_cells_per_donor_celltype,
        min_donors_per_celltype = min_donors_per_celltype
    )

                                        # Select enriched cell types descriptively
    enriched <- get_agtr1_enriched_celltypes(
        donor_celltype, metric = enrichment_metric, top_n = top_n_celltypes
    )
    keep_celltypes <- enriched$cell_type

                                        # Correlation results + FDR
    stats <- run_age_agtr1_by_celltype(donor_celltype, keep_celltypes)

                                        # Save outputs
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    if (write_donor) {
        data.table::fwrite(donor_celltype,
                           file.path(outdir, "donor_metadata.tsv"), sep = "\t")
    }
    data.table::fwrite(enriched, file.path(outdir, "agtr1_enriched_celltypes.tsv"),
                       sep = "\t")
    data.table::fwrite(stats, file.path(outdir, "age_agtr1_spearman_by_celltype.tsv"),
                       sep = "\t")

                                        # Plot
    plot_age_agtr1_facets(donor_celltype, keep_celltypes, outdir)
    return(list(donor_celltype = donor_celltype, enriched = enriched, stats = stats))
}

#### Main
                                        # Top 10 cell types by mean AGTR1 expression
res <- age_agtr1_analysis(
    outdir = "mean_expr", cell_type_key = "cell_type",
    enrichment_metric = "mean_expr", top_n_celltypes = 10,
    min_cells_per_donor_celltype = 20, min_donors_per_celltype = 3,
    write_donor = TRUE
)

                                        # Top 10 cell types by fraction AGTR1-positive
res <- age_agtr1_analysis(
    outdir = "frac_expr", cell_type_key = "cell_type",
    enrichment_metric = "frac_expr", top_n_celltypes = 10
)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
