## This script extracts the cell type proportions from R variable
suppressPackageStartupMessages({
    library(dplyr)
    library(SingleCellExperiment)
})

load_data <- function(){
    fn  <- here::here("inputs/hlca/_m/hlca_core.rds")
    sce <- Seurat::as.SingleCellExperiment(readRDS(fn))
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_coarse_for_GWAS_and_modeling
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type
    return(sce)
}

#### Main
colData(load_data()) |> as.data.frame() |>
    select(c("compartment", "cell_type")) |> distinct() |>
    mutate(clean_cell_type = janitor::make_clean_names(cell_type)) |>
    data.table::fwrite("hlca_lung.celltype_mapping.tsv", sep='\t')

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
