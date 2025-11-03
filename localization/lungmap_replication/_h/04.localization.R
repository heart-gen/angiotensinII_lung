## Analyze Lung 10X scRNA-seq data for angiotensin II.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(ggvenn)
    library(Seurat)
    library(SingleCellExperiment)
})

#### Functions
save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    sce <- zellkonverter::readH5AD("lungmap_transferred.h5ad")
    return(sce)
}
memSC <- memoise::memoise(load_data)

extract_angiotensinII <- function(){
    sce <- memSC()
    genes = c("AGTR1", "AGTR2")
    angiotensin2 <- sce[genes,]
    return(angiotensin2)
}
memAGTR <- memoise::memoise(extract_angiotensinII)

prepare_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Donor=colData(angiotensin2)$donor,
               Celltype=colData(angiotensin2)$celltype,
               Lineage=colData(angiotensin2)$lineage) |>
        tidyr::pivot_longer(-c(Donor, Celltype, Lineage),
                            names_to="Gene_Name",
                            values_to="Normalized Expression") |>
        mutate_if(is.character, as.factor)
    return(df)
}

prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Donor=colData(angiotensin2)$donor,
               Celltype=colData(angiotensin2)$celltype,
               Lineage=colData(angiotensin2)$lineage) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Donor, Celltype, Lineage),
                            names_to="Gene_Name",
                            values_to="Normalized Expression") |>
        mutate_if(is.character, as.factor)
    return(df)
}

                                        # Plotting functions
venn_diagrams <- function(){
    outfile <- "angiotensinII_venn_diagram"
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

plot_dotplot <- function(outdir="all_cells"){
    df.seurat <- as.Seurat(memAGTR(), counts = "counts", data = "logcounts")
    df.seurat$celltype <- forcats::fct_rev(factor(df.seurat$celltype))
    df.seurat$lineage  <- forcats::fct_rev(factor(df.seurat$lineage))
    df.seurat$active.ident <- df.seurat$celltype
    Idents(df.seurat)  <- "celltype"
    pp = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2")) +
        RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_celltype"), pp, 6, 7)
    qq = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="lineage") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_lineage"), qq, 6, 7)
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

calculate_proportion <- function(df, outdir){
    dt = df |> group_by(Gene_Name) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    print(dt)
    dt = df |> group_by(Gene_Name, Lineage) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    print(dt); plot_bar(dt, "Lineage", 45, 3, 5, outdir)
    dt = df |> group_by(Gene_Name, Celltype) |>
        summarise(Proportion=sum(`Normalized Expression` > 0)/n(),
                  Total_Cells=sum(`Normalized Expression` > 0))
    print(arrange(dt, desc(Proportion)) |> head(10));

    plot_bar(dt, "Celltype", 90, 12, 8, outdir)
    outfile <- paste0(outdir,"/lung_angiotensinII_proportions.tsv")
    dt |> data.table::fwrite(outfile, sep="\t")
}

### Main script section
                                        # Plot venn diagram
venn_diagrams()
                                        # All cells
outdir = "all_cells"; dir.create(outdir)

plot_dotplot()
df1 <- prepare_data()
df1 |> data.table::fwrite("normalized_expression.txt.gz", sep='\t', 
                          compress="auto")
calculate_proportion(df1, outdir)

                                        # Filtered
outdir = "filtered_cells"; dir.create(outdir)
df2 <- prepare_filtered_data()
calculate_proportion(df2, outdir)


#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
