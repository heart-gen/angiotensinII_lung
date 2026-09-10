## Build `inputs/ipf/_m/ipf_dataset.h5ad` from the GEO GSE136831 (Adams/Kaminski)
## download: 312,928 cells x 45,947 genes of RAW counts.
##
## RELOCATED 2026-09-10. This script used to be
## `disease_association/ipf_analysis/_h/02.generate_ipf_data.R`, exempted by name
## from that module's retirement because four current modules read its output
## (upstream defect P1-18). When the disease analyses moved to
## heart-gen/lung-pericyte-analysis, three of those four readers went with them,
## leaving `basement_membrane/_h/step_3.sh` -- the COPD basement-membrane arm --
## reading a path with no builder behind it. Rather than reach into the other
## repository, the builder came here, where the GEO download it consumes already
## lives. The disease repository carries its own copy.
##
## Run from inputs/ipf/_m via `../_h/step_1.sh`. The download must exist first:
## `bash ../_h/submit.sh` (login node -- compute nodes have no internet).


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
