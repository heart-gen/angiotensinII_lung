## Subset lung single-cell data

suppressPackageStartupMessages({
    library(here)
    library(SingleCellExperiment)
})

subset_N_save <- function(){
					# Load data
    fn <- here("inputs/hlca/_m/hlca_full.h5ad")
    sce <- zellkonverter::readH5AD(fn)
    names(assays(sce)) <- c("counts", "soupX")
                                        # Remove missing annotations (level 4)
    sce <- sce[, !is.na(colData(sce)$ann_level_4)]
                                        # Remove unknown cells
    sce <- sce[, !colData(sce)$ann_level_4 %in% c("None", "Unknown")]
                                        # Update annotation
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_level_4
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type
    ##sce <- scuttle::logNormCounts(sce)
    ##sce <- scran::computeSumFactors(sce)
                                        # Filter celltypes
    peri <- sce[,colData(sce)$subclusters %in% c("Pericytes")]
    stro <- sce[,colData(sce)$compartment %in% c("Stroma")]
                                        # Write as H5AD
    zellkonverter::writeH5AD(peri,
                             file="peri.hlca_full.dataset.h5ad")
    zellkonverter::writeH5AD(stro,
                             file="stroma.hlca_full.dataset.h5ad")
}

subset_at2 <- function(){
					# Load data
    fn <- here("inputs/hlca/_m/hlca_core.h5ad")
    sce <- zellkonverter::readH5AD(fn)
    names(assays(sce)) <- c("counts")
                                        # Update annotation
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_level_4
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type
                                        # Filter celltypes
    at2  <- sce[,colData(sce)$subclusters %in% c("AT2")]
    peri <- sce[,colData(sce)$subclusters %in% c("Pericytes")]
    stro <- sce[,colData(sce)$compartment %in% c("Stroma")]
                                        # Write as H5AD
    zellkonverter::writeH5AD(at2,
                             file="at2.hlca_core.dataset.h5ad")
    zellkonverter::writeH5AD(peri,
                             file="peri.hlca_core.dataset.h5ad")
    zellkonverter::writeH5AD(stro,
                             file="stroma.hlca_core.dataset.h5ad")
}

#### Main
subset_N_save()
subset_at2()
    
#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

