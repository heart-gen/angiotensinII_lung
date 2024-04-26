## Analyze Lung SmartSeq2 scRNA-seq data for angiotensin II.

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggpubr)
    library(ggvenn)
    library(Seurat)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    fn  <- here("inputs/copd/_m/COPD_dataset.RData")
    load(fn, verbose=TRUE)
    sce <- UpdateSeuratObject(lungsoup) |>
        as.SingleCellExperiment()
    colData(sce)$cell_type = colData(sce)$cell2
    colData(sce)$disease = colData(sce)$disease.ident
    colData(sce)$patient = colData(sce)$subject.ident
    colLabels(sce) <- colData(sce)$subject_celltype2
    return(sce)
}
memSC <- memoise::memoise(load_data)

extract_angiotensinII <- function(){
                                        # Load data
    sce <- memSC()
    genes = c("AGTR1", "AGTR2")
    angiotensin2 <- sce[genes,]
    return(angiotensin2)
}
memAGTR <- memoise::memoise(extract_angiotensinII)

venn_diagrams <- function(){
    outfile = "angiotensinII_venn_diagram"
    x = list(
        AGTR1 = counts(memAGTR()) |> t() |> as.data.frame() |>
            tibble::rownames_to_column("Cell_ID") |>
            filter(AGTR1 > 0) |> select(Cell_ID) |> unlist(),
        AGTR2 = counts(memAGTR()) |> t() |> as.data.frame() |>
            tibble::rownames_to_column("Cell_ID") |>
            filter(AGTR2 > 0) |> select(Cell_ID) |> unlist()
    )
    vv <- ggvenn(x, fill_color=get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_ggplots(tolower(outfile), vv, 5, 5)
}

plot_dotplot <- function(outdir, FILTER=FALSE){
    df.seurat <- as.Seurat(memAGTR(), counts = "counts", data = "logcounts")
    df.seurat$cell_type = forcats::fct_rev(factor(df.seurat$cell_type))
    df.seurat$disease = forcats::fct_rev(factor(df.seurat$disease))
    df.seurat$active.ident <- df.seurat$cell_type
    Idents(df.seurat) <- "cell_type"
    if(FILTER){
        df.seurat <- subset(x=df.seurat, subset=AGTR1 >0 | AGTR2 > 0)
        h = 6; w = 6
    } else {
        h = 10; w = 6
    }
    pp = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2")) +
        RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_cell_annotation"), pp, w, h)
    qq = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="disease") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_disease"), qq, 5, 5)
}

prepare_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Disease=colData(angiotensin2)$disease) |>
        mutate_if(is.character, as.factor) |>
        tidyr::pivot_longer(-c(Cell_Annotation, Patient, Disease),
                            names_to="Gene_Name",
                            values_to="Normalized Expression")
    return(df)
}

prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Disease=colData(angiotensin2)$disease) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Cell_Annotation, Patient, Disease),
                            names_to="Gene_Name",
                            values_to="Normalized Expression")
    return(df)
}

plot_box_stats <- function(dt, xlab, outdir){
    mycomps <- list(c("Control", "COPD"))
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
    mycomps <- list(c("Control", "COPD")); xlab = "Disease"
    celltypes1 <- c("AT2 A", "AT2 B", "Pericytes", "SMC")
    celltypes2 <- c("AT2 B", "Fibroblast Adventitial",
                   "Fibroblast Alveolar", "Pericytes")
    outfile = paste0(outdir,"/normalized_expression_boxplot_celltype_by_", xlab)
    bxp1 <- dt |> group_by_at(c("Patient", xlab, "Cell_Annotation",
                                 "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        filter(Cell_Annotation %in% celltypes1) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, add="jitter",
                  facet.by=c("Gene_Name", "Cell_Annotation"), ylim=c(0, 2.75),
                  legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=4) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(paste0(tolower(outfile), ".v2"), bxp1, 12, 8)
    bxp2 <- dt |> group_by_at(c("Patient", xlab, "Cell_Annotation", "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        filter(Cell_Annotation %in% celltypes2) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, add="jitter",
                  facet.by=c("Gene_Name", "Cell_Annotation"), ylim=c(0, 2.75),
                  legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=4) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp2, 12, 8)
}

plot_box_stats_celltype_all <- function(dt, outdir){
    mycomps <- list(c("Control", "COPD")); xlab = "Disease"
    outfile = paste0(outdir,"/normalized_expression_boxplot_celltype_", xlab)
    bxp = dt |> group_by_at(c("Patient", xlab, "Cell_Annotation", "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, add="jitter",
                  facet.by=c("Gene_Name", "Cell_Annotation"), ylim=c(0, 2.75),
                  legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", palette="npg",
                  panel.labs.font=list(face='bold'), add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=6) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp, 36, 6)
}

generate_boxplots <- function(df, outdir){
                                        # Statistics
    plot_box_stats(df, "Disease", outdir)
    plot_box_stats_celltype(df, outdir)
    ##plot_box_stats_celltype_all(dt, outdir)
}

### Main script section
venn_diagrams()
                                        # All cells
outdir = "all_cells"
dir.create(outdir)
plot_dotplot(outdir, FALSE)
df0 <- prepare_data()
df0 |> data.table::fwrite("normalized_expression.tsv", sep='\t')
generate_boxplots(df0, outdir)

                                        # Filtered cells
outdir = "filter_cells"
dir.create(outdir)
plot_dotplot(outdir, TRUE)
df <- prepare_filtered_data()
generate_boxplots(df, outdir)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
