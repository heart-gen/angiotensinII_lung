## Generate IPF dataset

suppressPackageStartupMessages({
    library(here)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Main

fn1 = here::here("inputs/lungmap/_m/GSE161382_counts_matrix.txt.gz")
counts <- data.table::fread(fn1) |> tibble::column_to_rownames("UID")
## Load meta data
fn2 = here::here("inputs/lungmap/_m/GSE161382_metadata.txt.gz")
meta <- data.table::fread(fn2) |> tibble::column_to_rownames("V1")
## Generate RSE object
sce <- SingleCellExperiment(list(counts=counts),
                            colData=meta,
                            metadata=list(study="GSE161382"))
sce <- scuttle::logNormCounts(sce)
sce <- scran::computeSumFactors(sce)
colLabels(sce) <- colData(sce)$celltype

fn = here::here("inputs/lungmap/_m/GSE161382_UMAP_coord.tsv.gz")
umap <- data.table::fread(fn) |> tibble::column_to_rownames("V1")
reducedDims(sce) <- list(UMAP=umap)
                                        # Write as H5AD
zellkonverter::writeH5AD(sce, file="lungmap_dataset.h5ad")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
