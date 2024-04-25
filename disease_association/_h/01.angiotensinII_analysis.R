## Analyze Lung 10X scRNA-seq data for angiotensin II.
suppressPackageStartupMessages({
    library(here)
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

                                        # Prepare single cell data
                                        # for angiotensin II
load_data <- function(){# Gene annotation
    fn0 <- here("inputs/ipf/_m/GSE136831_AllCells.GeneIDs.txt.gz")
    gnames <- data.table::fread(fn0)
                                        # Counts
    fn1 <- here("inputs/ipf/_m/GSE136831_RawCounts_Sparse.mtx.gz")
    counts <- Matrix::readMM(gzfile(fn1))
                                        # Meta data
    fn2 <- here("inputs/ipf/_m",
                "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
    meta <- data.table::fread(fn2)
                                        # Generate RSE object
    sce <- SingleCellExperiment(list(counts=counts),
                                colData=meta, rowData=gnames,
                                metadata=list(study="GSE136831"))
    rownames(sce)  <- gnames$HGNC_EnsemblAlt_GeneID
    colnames(sce)  <- meta$CellBarcode_Identity
    counts(sce)    <- as(counts(sce), "CsparseMatrix")
                                        # Update labels
    colData(sce)$cell_type <- colData(sce)$Manuscript_Identity
    colData(sce)$location  <- colData(sce)$CellType_Category
    colData(sce)$disease   <- colData(sce)$Disease_Identity
    colData(sce)$patient   <- colData(sce)$Subject_Identity
    colLabels(sce)         <- colData(sce)$Subclass_Cell_Identity
                                        # Normalize
    sce <- scran::computeSumFactors(sce, clusters=sce$cell_type)
    sce <- scuttle::logNormCounts(sce)
    logcounts(sce) <- as(logcounts(sce), "CsparseMatrix")
    return(sce)
}
memSC <- memoise::memoise(load_data)

add_qc <- function(sce){
    is.mito <- grep("^MT", rownames(sce))
    sce     <- scuttle::addPerCellQC(sce, subsets=list(Mito=is.mito))
    sce     <- scuttle::addPerFeatureQC(sce)
    return(sce)
}

filter_qc <- function(sce){
                                        # Add QC
    sce <- add_qc(sce)
                                        # Select cells library size < 1000
    qc.lib <- sce$nUMI > 1000
    lib    <- paste("Removing", sum(!qc.lib),
                    "cells due to small (< 1000) lib size!")
    print(lib)
                                        # Select cells with > 25% mitochondria
    qc.mito <- sce$subsets_Mito_percent < 25
    mito    <- paste("Removing", sum(!qc.mito),
                     "cells due large percentage (> 25%) of mitochondria genes!")
    print(mito)
                                        # Remove cells based on lib size
                                        # and percent mitochondria
    discard <- qc.mito | qc.lib
    sce <- sce[, discard]
    return(sce)
}

extract_angiotensinII <- function(){
                                        # Load data
    sce <- memSC()
                                        # Filter QC
    sce <- filter_qc(sce)
                                        # Select angiotensin II genes
    genes = c("AGTR1", "AGTR2")
    angiotensin2 <- sce[genes,]
    return(angiotensin2)
}
memAGTR <- memoise::memoise(extract_angiotensinII)

venn_diagrams <- function(){
    outfile = "angiotensinII_venn_diagram"
    x = list(
        AGTR1 = counts(memAGTR()) |> t() |> as.data.frame() |>
            tibble::rownames_to_column("Cell_ID") |> filter(AGTR1 > 0) |>
            select(Cell_ID) |> unlist(),
        AGTR2 = counts(memAGTR()) |> t() |> as.data.frame() |>
            tibble::rownames_to_column("Cell_ID") |> filter(AGTR2 > 0) |>
            select(Cell_ID) |> unlist()
    )
    vv <- ggvenn(x, fill_color=get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_ggplots(tolower(outfile), vv, 5, 5)
}

plot_dotplot <- function(outdir, FILTER=FALSE){
    df.seurat <- as.Seurat(memAGTR(), counts = "counts", data = "logcounts")
    df.seurat$cell_type    <- forcats::fct_rev(factor(df.seurat$cell_type))
    df.seurat$disease      <- forcats::fct_rev(factor(df.seurat$disease))
    df.seurat$active.ident <- df.seurat$cell_type
    Idents(df.seurat)      <- "cell_type"
    if(FILTER){
        df.seurat <- subset(x=df.seurat, subset=AGTR1 >0 | AGTR2 > 0)
        h = 6; w = 6
    } else {
        h = 10; w = 6
    }
    pp = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2")) +
        RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_cell_annotation"),
                 pp, w, h)
    qq = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="disease") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_disease"), qq, 5, 5)
}

prepare_data <- function(){
    angiotensin2 <- memAGTR()
    df <- logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Disease=colData(angiotensin2)$disease) |>
        mutate_if(is.character, as.factor) |>
        tidyr::pivot_longer(-c(Cell_Annotation, Patient, Disease),
                            names_to="Gene_Name",
                            values_to="Normalized Expression")
    levels(df$Disease) <- c("Control", "COPD", "IPF")
    return(df)
}

prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df <- logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Disease=colData(angiotensin2)$disease) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Cell_Annotation, Patient, Disease),
                            names_to="Gene_Name",
                            values_to="Normalized Expression")
    levels(df$Disease) <- c("Control", "COPD", "IPF")
    return(df)
}

plot_box_stats <- function(dt, xlab, outdir){
    mycomps <- list(c("Control", "COPD"), c("Control", "IPF"),
                    c("COPD", "IPF"))
    outfile = paste0(outdir,"/normalized_expression_boxplot_", xlab)
    bxp = dt |> group_by_at(c("Patient", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, facet.by="Gene_Name",
                  add="jitter", legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=2) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp, 6, 6)
}

plot_box_stats_celltype <- function(dt, outdir){
    mycomps   <- list(c("Control", "COPD"), c("Control", "IPF"),
                      c("COPD", "IPF"));
    xlab      <- "Disease"
    celltypes <- c("ATI", "ATII", "Fibroblast", "Pericyte", "Myofibroblast")
    outfile = paste0(outdir,"/normalized_expression_boxplot_celltype_by_", xlab)
    bxp = dt |> group_by_at(c("Patient",xlab,"Cell_Annotation","Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        filter(Cell_Annotation %in% celltypes) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, add="jitter",
                  facet.by=c("Gene_Name", "Cell_Annotation"), ylim=c(0, 4.25),
                  legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=4) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp, 12, 8)
}

plot_box_stats_celltype_all <- function(dt, outdir){
    mycomps <- list(c("Control", "COPD"), c("Control", "IPF"),
                    c("COPD", "IPF"));
    xlab    <- "Disease"
    outfile <- paste0(outdir,"/normalized_expression_boxplot_celltype_", xlab)
    bxp = dt |> group_by_at(c("Patient",xlab,"Cell_Annotation","Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, add="jitter",
                  facet.by=c("Gene_Name", "Cell_Annotation"), ylim=c(0, 4.25),
                  legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=6) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp, 36, 6)
}

generate_boxplots <- function(dt, outdir){
                                        # Statistics
    plot_box_stats(dt, "Disease", outdir)
    plot_box_stats_celltype(dt, outdir)
    plot_box_stats_celltype_all(dt, outdir)
}

### Main script section
                                        # Plot venn diagrams by age group
venn_diagrams()
                                        # All cells
outdir <- "all_cells"
dir.create(outdir)
plot_dotplot(outdir, FALSE)
df0    <- prepare_data()
generate_boxplots(df0, outdir)
df0 |> data.table::fwrite("normalized_expression.tsv", sep='\t')
                                        # Filtered cells
outdir <- "filter_cells"
dir.create(outdir)
plot_dotplot(outdir, TRUE)
df     <- prepare_filtered_data()
generate_boxplots(df, outdir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
