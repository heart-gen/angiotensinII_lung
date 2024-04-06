## Analyze HLCA version 2 scRNA-seq data for angiotensin II.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggvenn)
    library(Seurat)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    fn  <- here::here("inputs/hlca/_m/hlca_core.rds")
    sce <- zellkonverter::readH5AD(fn)
    names(assays(sce)) <- "counts"
    sce <- scuttle::logNormCounts(sce)
    sce <- scran::computeSumFactors(sce)
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_coarse_for_GWAS_and_modeling
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type
    return(sce)
}
memSC <- memoise::memoise(load_data)

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

extract_angiotensinII <- function(){
                                        # Load data
    sce <- memSC()
                                        # Add QC
    sce <- add_qc(sce)
                                        # Filter QC
    sce <- filter_qc(sce)
                                        # Update rownames to gene name
    rownames(sce) <- rowData(sce)[, 2]
                                        # Subset for angiotensin receptors
    angiotensin2  <- sce[c("AGTR1", "AGTR2"),]
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
    df.seurat$cell_type   <- forcats::fct_rev(factor(df.seurat$cell_type))
    df.seurat$compartment <- forcats::fct_rev(factor(df.seurat$compartment))
    df.seurat$subclusters <- forcats::fct_rev(factor(df.seurat$subclusters))
    df.seurat$clusters    <- forcats::fct_rev(factor(df.seurat$clusters))
    df.seurat$active.ident <- df.seurat$cell_type
    Idents(df.seurat) <- "cell_type"
    if(FILTER){
        df.seurat <- subset(x=df.seurat, subset=AGTR1 > 0 | AGTR2 > 0)
        h = 6; w = 6
    } else {
        h = 10; w = 7.5
    }
    pp = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2")) +
        RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_cell_annotation"), pp, w, h)
    rr = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="compartment") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_compartment"), rr, 4, 5)
    qq = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="subclusters") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_subcluster"), qq, w, h+2.5)
    gg = DotPlot(object=df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="clusters") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_clusters"), gg, w, h)
}

prepare_data <- function(){
    angiotensin2 <- memAGTR()
    df <- logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Subcluster=colData(angiotensin2)$subclusters,
               Cluster=colData(angiotensin2)$clusters,
               Compartment=colData(angiotensin2)$compartment) |>
        tidyr::pivot_longer(-c(Cell_Annotation, Patient, Compartment,
                               Subcluster, Cluster),
                            names_to="Gene_Name",
                            values_to="Normalized Expression")
    return(df)
}

plot_bar <- function(dt, xlab, r, w, h, outdir){
    outfile = paste0(outdir,"/proportion_plot_", xlab)
    bar = ggbarplot(dt, x=xlab, y="Proportion", fill="Gene_Name",
                    lab.pos="out", color="Gene_Name", palette="Paired",
                    label=FALSE, xlab="",
                    ggtheme=theme_pubr(base_size=15)) +
        rotate_x_text(r)
    save_ggplots(tolower(outfile), bar, w, h)
}

calculate_proportion <- function(df, outdir, FILTER=FALSE){
    ## Print proportions to screen
    dt = df |> group_by(Gene_Name) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    print(dt)
    dt = df |> group_by(Gene_Name, Compartment) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    print(dt); plot_bar(dt, "Compartment", 45, 3, 5, outdir)
    dt = df |> group_by(Gene_Name, Subcluster) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    if(FILTER){
        print(arrange(dt, desc(Total_Cells)) |> head(10));
    } else {
        print(arrange(dt, desc(Proportion)) |> head(10));
    }
    plot_bar(dt, "Subcluster", 90, 12, 8, outdir)
    outfile <- paste0(outdir,"/lung_angiotensinII_proportions.subclusters.tsv")
    dt |> data.table::fwrite(outfile, sep="\t")
    ## Save free annotation proportion
    dt = df |> group_by(Gene_Name, Cell_Annotation) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    if(FILTER){
        print(arrange(dt, desc(Total_Cells)) |> head(10));
    } else {
        print(arrange(dt, desc(Proportion)) |> head(10));
    }
    plot_bar(dt, "Cell_Annotation", 90, 12, 8, outdir)
    outfile <- paste0(outdir,"/lung_angiotensinII_proportions.tsv")
    dt |> data.table::fwrite(outfile, sep="\t")
}

### Main script section
venn_diagrams()
                                        # All cells
outdir = "all_cells"
dir.create(outdir)
plot_dotplot(outdir, FALSE)
df0 <- prepare_data()
df0 |> data.table::fwrite("normalized_expression.txt.gz", sep='\t', 
                          compress="auto")
calculate_proportion(df0, outdir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
