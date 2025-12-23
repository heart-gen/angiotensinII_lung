## --- Add: donor-level age ~ AGTR1 analysis (within AGTR1-enriched cell types) ---

suppressPackageStartupMessages({
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
})

first_non_na <- function(x) {
    x <- x[!is.na(x)]
    if(length(x) == 0) NA else x[[1]]
}

prepare_age_agtr1_donor_table <- function(
    sce, cell_type_key = "cell_type", donor_key = "patient",
    age_key = "age_or_mean_of_age_range", sex_key = "sex",
    disease_key = "disease", ethnicity_key = "self_reported_ethnicity",
    gene = "AGTR1", min_cells_per_donor_celltype = 20, 
    min_donors_per_celltype = 3
) {
    stopifnot(gene %in% rownames(sce))
    cd <- as.data.frame(SummarizedExperiment::colData(sce))

    # Extract AGTR1 logcounts for all cells
    ag <- as.numeric(SingleCellExperiment::logcounts(sce)[gene, ])
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

    # Require enough cells per donor within a cell type
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
          pval = suppressWarnings(cor.test(age, AGTR1_mean, method = "spearman")$p.value),
          .groups = "drop"
      ) |>
      mutate(fdr = p.adjust(pval, method = "fdr")) |>
      arrange(fdr, pval)
    return(res)
}

plot_age_agtr1_facets <- function(
    donor_celltype, keep_celltypes, outdir, filename = "age_vs_agtr1_by_celltype"){
    dfp <- donor_celltype |>
      filter(cell_type %in% keep_celltypes) |>
      mutate(cell_type = forcats::fct_reorder(cell_type, AGTR1_mean, .fun = median))

    p <- ggplot(dfp, aes(x = age, y = AGTR1_mean)) +
      geom_point(alpha = 0.6, size = 1.2) +
      geom_smooth(method = "loess", se = TRUE, color = "black") +
      facet_wrap(~ cell_type, scales = "free_y") +
      theme_classic(base_size = 14) +
      labs(x = "Donor age (years)", y = "Donor-mean AGTR1 (logcounts)",
           title = "Age-associated AGTR1 expression in AGTR1-enriched cell types")

    save_ggplots(file.path(outdir, filename), p, w = 12, h = 8)
}

age_agtr1_analysis <- function(
  outdir = "all_cells",
  cell_type_key = "cell_type",
  donor_key = "patient",
  age_key = "age_or_mean_of_age_range",
  gene = "AGTR1",
  top_n_celltypes = 10,
  enrichment_metric = c("mean_expr", "frac_expr"),
  min_cells_per_donor_celltype = 20,
  min_donors_per_celltype = 3
){
  enrichment_metric <- match.arg(enrichment_metric)

  # Use the full HLCA object (QC-filtered) rather than memAGTR() (which only has AGTR1/AGTR2 rows)
  sce <- memSC()
  sce <- add_rowdata(sce)
  sce <- add_qc(sce)
  sce <- filter_qc(sce)
  rownames(sce) <- rowData(sce)[, "feature_name"]

  # Build donor x cell_type table for AGTR1
  donor_celltype <- prepare_age_agtr1_donor_table(
    sce,
    cell_type_key = cell_type_key,
    donor_key = donor_key,
    age_key = age_key,
    gene = gene,
    min_cells_per_donor_celltype = min_cells_per_donor_celltype,
    min_donors_per_celltype = min_donors_per_celltype
  )

  # Select enriched cell types descriptively
  enriched <- get_agtr1_enriched_celltypes(
    donor_celltype,
    metric = enrichment_metric,
    top_n = top_n_celltypes
  )
  keep_celltypes <- enriched$cell_type

  # Correlation results + FDR
  stats <- run_age_agtr1_by_celltype(donor_celltype, keep_celltypes)

  # Save outputs
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(enriched, file.path(outdir, "agtr1_enriched_celltypes.tsv"), sep = "\t")
  data.table::fwrite(stats, file.path(outdir, "age_agtr1_spearman_by_celltype.tsv"), sep = "\t")

  # Plot
  plot_age_agtr1_facets(donor_celltype, keep_celltypes, outdir)
  invisible(list(donor_celltype = donor_celltype, enriched = enriched, stats = stats))
}

#### Main
# Example 1: top 10 cell types by mean AGTR1 expression
age_agtr1_analysis(
  outdir = "all_cells",
  cell_type_key = "cell_type",
  enrichment_metric = "mean_expr",
  top_n_celltypes = 10,
  min_cells_per_donor_celltype = 20,
  min_donors_per_celltype = 3
)

# Example 2: top 10 cell types by fraction AGTR1-positive
# age_agtr1_analysis(
#   outdir = "all_cells",
#   cell_type_key = "cell_type",
#   enrichment_metric = "frac_expr",
#   top_n_celltypes = 10
# )

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
