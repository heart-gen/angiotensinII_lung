## Generate IPF dataset

suppressPackageStartupMessages({
    library(here)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Main

fn0    <- here("inputs/ipf/_m/GSE136831_AllCells.GeneIDs.txt.gz")
gnames <- data.table::fread(fn0)
                                        # Counts
fn1    <- here("inputs/ipf/_m/GSE136831_RawCounts_Sparse.mtx.gz")
counts <- Matrix::readMM(gzfile(fn1))
                                        # Meta data
fn2  <- here("inputs/ipf/_m",
             "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
meta <- data.table::fread(fn2)
                                        # Generate RSE object
sce  <- SingleCellExperiment(list(counts=counts),
                             colData=meta, rowData=gnames,
                             metadata=list(study="GSE136831"))
rownames(sce) <- gnames$HGNC_EnsemblAlt_GeneID
colnames(sce) <- meta$CellBarcode_Identity
counts(sce)   <- as(counts(sce), "CsparseMatrix")
                                        # Update labels
colData(sce)$cell_type <- colData(sce)$Manuscript_Identity
colData(sce)$location  <- colData(sce)$CellType_Category
colData(sce)$disease   <- colData(sce)$Disease_Identity
colData(sce)$patient   <- colData(sce)$Subject_Identity
colLabels(sce)         <- colData(sce)$Subclass_Cell_Identity
                                        # Write as H5AD
zellkonverter::writeH5AD(sce, file="ipf_dataset.h5ad")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
