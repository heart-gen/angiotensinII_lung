## Analyze Lung 10X scRNA-seq data for angiotensin II.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggvenn)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Functions
                                        # Basic plotting helper function
save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

                                        # Prepare single cell data for angiotensin II
load_data <- function(){
    ## Load counts
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
    return(sce)
}
memSC <- memoise::memoise(load_data)

add_qc <- function(sce){
    is.mito <- grep("^MT", rownames(sce))
    sce <- scuttle::addPerCellQC(sce, subsets=list(Mito=is.mito))
    sce <- scuttle::addPerFeatureQC(sce)
    return(sce)
}

add_umap <- function(sce){
    fn = here::here("inputs/lungmap/_m/GSE161382_UMAP_coord.tsv.gz")
    umap <- data.table::fread(fn) |> tibble::column_to_rownames("V1")
    reducedDims(sce) <- list(UMAP=umap)
    return(sce)
}

filter_qc <- function(sce){
    ## Remove cells with > 25% mitochondria and library size < 1000
    qc.lib <- sce$sum > 1000
    lib = paste("Removing", sum(!qc.lib),
                "cells due to small (< 1000) lib size!")
    print(lib)
    qc.mito <- sce$subsets_Mito_percent < 25
    mito = paste("Removing", sum(!qc.mito),
                 "cells due large percentage (> 25%) of mitochondria genes!")
    print(mito)
    discard <- qc.mito | qc.lib
    sce <- sce[, discard]
    return(sce)
}

extract_angiotensinII <- function(){
                                        # Load data
    sce <- memSC()
    sce <- add_umap(sce)
                                        # Add QC
    sce <- add_qc(sce)
                                        # Filter QC
    sce <- filter_qc(sce)
                                        # Select angiotensin II genes
    genes = c("AGTR1", "AGTR2")
    angiotensin2 <- sce[genes,]
    angiotensin2$age <- factor(angiotensin2$age, levels=c("31wk","3yr","31yr"))
    return(angiotensin2)
}
memAGTR <- memoise::memoise(extract_angiotensinII)

                                        # Prepare data functions
prepare_counts <- function(){
    angiotensin2 <- memAGTR()
    df = counts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Donor=colData(angiotensin2)$donor,
               Celltype=colData(angiotensin2)$celltype,
               Lineage=colData(angiotensin2)$lineage,
               Age=colData(angiotensin2)$age) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Donor, Celltype, Lineage, Age),
                            names_to="Gene_Name", values_to="Counts") |>
        mutate_if(is.character, as.factor) |> as.data.frame()
    return(df)
}

subset_stromal <- function(){
    outdir = "all_cells"
    sce <- memAGTR()[,colData(memAGTR())$lineage == "Mesenchymal"]
    df.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
    df.seurat$celltype = forcats::fct_rev(factor(df.seurat$celltype))
    df.seurat$active.ident <- df.seurat$celltype
    Idents(df.seurat) <- "celltype"
                                        # feature plots
    pp = FeaturePlot(df.seurat, features=c("AGTR1"), split.by="age",
                     max.cutoff=3, cols=c("grey", "red"),
                     reduction="UMAP")
    save_ggplots(paste0(outdir,"/featureplot.mesenchymal.celltype"),
                 pp, 12, 5)
                                        # dot plots
    ncols <- ggpubr::get_palette("aaas", 3)
    qq = DotPlot(df.seurat, features = c("AGTR1"), cols=ncols,
                 split.by="age") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot.mesenchymal.celltype"),
                 qq, 6, 10)
                                        # dim plot
    rr = DimPlot(df.seurat, label=FALSE, repel=FALSE)
    save_ggplots(paste0(outdir, "/umap.mesenchymal.celltype"),
                 rr, 7, 6)
}

subset_epithelial <- function(){
    outdir = "all_cells"
    sce <- memAGTR()[,colData(memAGTR())$lineage == "Epithelial"]
    df.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")
    df.seurat$celltype = forcats::fct_rev(factor(df.seurat$celltype))
    df.seurat$active.ident <- df.seurat$celltype
    Idents(df.seurat) <- "celltype"
                                        # feature plots
    pp = FeaturePlot(df.seurat, features=c("AGTR2"), split.by="age",
                     max.cutoff=3, cols=c("grey", "red"),
                     reduction="UMAP")
    save_ggplots(paste0(outdir,"/featureplot.epithelial.celltype"),
                 pp, 12, 5)
                                        # dot plots
    ncols <- ggpubr::get_palette("aaas", 3)
    qq = DotPlot(df.seurat, features = c("AGTR2"), cols=ncols,
                 split.by="age") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot.epithelial.celltype"),
                 qq, 6, 10)
                                        # dim plot
    rr = DimPlot(df.seurat, label=FALSE, repel=FALSE)
    save_ggplots(paste0(outdir, "/umap.epithelial.celltype"),
                 rr, 7, 6)
}

#### Main
subset_stromal()
subset_epithelial()

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
