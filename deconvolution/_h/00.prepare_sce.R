## Prepare data for deconvolution using HLCA version 2.

suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(BayesPrism)
    library(SingleCellExperiment)
})

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
memSC <- memoise::memoise(load_data)

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
    ## Find outliers
                                        # Mitochondria percentage
    print(summary(sce$subsets_mito_percent))
    ## The percentage is very small, so will not add additional filtering
    qc_stats <- scuttle::perCellQCFilters(sce)
    print(colSums(as.matrix(qc_stats)))
                                        # Discard outliers
    sce      <- sce[, !qc_stats$discard]
    return(sce)
}

#### Main
                                        # Prepare SCE
sce <- memSC()
sce <- add_rowdata(sce)
                                        # Add QC
sce <- add_qc(sce)
                                        # Filter QC
sce <- filter_qc(sce)
                                        # Save results
save(sce, file="scRNA_HLCA_version2.RData")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
