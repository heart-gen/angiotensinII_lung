## Generate IPF dataset

suppressPackageStartupMessages({
    library(here)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Main

                                        # Load counts
fn1 = here("inputs/lungmap/_m/GSE161382_counts_matrix.txt.gz")
counts <- data.table::fread(fn1) |> tibble::column_to_rownames("UID")
                                        # Load meta data
fn2 = here("inputs/lungmap/_m/GSE161382_metadata.txt.gz")
meta <- data.table::fread(fn2) |> tibble::column_to_rownames("V1")
                                        # Generate RSE object
sce <- SingleCellExperiment(list(counts=counts),
                            colData=meta,
                            metadata=list(study="GSE161382"))

                                        # Update labels
colData(sce)$cell_type <- colData(sce)$celltype
colData(sce)$location  <- colData(sce)$lineage
colData(sce)$patient   <- colData(sce)$donor
colLabels(sce)         <- colData(sce)$celltype

                                        # Write as H5AD
zellkonverter::writeH5AD(sce, file="lungmap_dataset.h5ad")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
