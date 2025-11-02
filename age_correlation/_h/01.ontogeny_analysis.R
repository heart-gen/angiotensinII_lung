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
        mutate_if(is.character, as.factor)
    return(df)
}

prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Donor=colData(angiotensin2)$donor,
               Celltype=colData(angiotensin2)$celltype,
               Lineage=colData(angiotensin2)$lineage,
               Age=colData(angiotensin2)$age) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Donor, Celltype, Lineage, Age),
                            names_to="Gene_Name",
                            values_to="Normalized Expression") |>
        mutate_if(is.character, as.factor)
    return(df)
}

                                        # Plotting functions
venn_diagrams <- function(age_var){
    outfile = paste0("angiotensinII_venn_diagram_", age_var)
    x = list(
        AGTR1 = counts(memAGTR()) |> t() |> as.data.frame() |>
            mutate(Age=colData(memAGTR())$age) |>
            tibble::rownames_to_column("Cell_ID") |>
            filter(AGTR1 > 0, Age == age_var) |>
            select(Cell_ID) |> unlist(),
        AGTR2 = counts(memAGTR()) |> t() |> as.data.frame() |>
            mutate(Age=colData(memAGTR())$age) |>
            tibble::rownames_to_column("Cell_ID") |>
            filter(AGTR2 > 0, Age == age_var) |>
            select(Cell_ID) |> unlist()
    )
    vv <- ggvenn(x, fill_color=get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_ggplots(tolower(outfile), vv, 5, 5)
}

plot_dotplot <- function(){
    outdir = "all_cells"
    dir.create(outdir)
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
    rr = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="age") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_age"), rr, 6, 7)
}

plot_line_age <- function(df1, df2, outdir){
    line1 <- df1 |> group_by_at(c("Gene_Name", "Donor")) |>
        summarize(tot_ct = sum(Counts), Age=Age) |> distinct() |>
        ggline(x="Age", y="tot_ct", add=c("mean_se", "jitter"),
               linetype="Gene_Name", ylab="Counts",
               add.params=list(alpha=0.8),
               ggtheme=theme_pubr(base_size=15, border=TRUE))
    line2 <- df2 |> group_by_at(c("Gene_Name", "Donor")) |>
        summarize(mean_exp = mean(`Normalized Expression`, na.rm=TRUE),
                  Age=Age) |> distinct() |>
        ggline(x="Age", y="mean_exp", add=c("mean_se", "jitter"),
               linetype="Gene_Name", ylab="Normalized Expression", 
               add.params=list(alpha=0.8),
               ggtheme=theme_pubr(base_size=15, border=TRUE))
    save_ggplots(paste0(outdir, "/line_graph_counts"), line1, 5, 5)
    save_ggplots(paste0(outdir, "/line_graph_expr"), line2, 5, 5)
}

plot_line_variable <- function(df1, df2, xlab, cols, w, outdir){
    line1 <- df1 |> group_by_at(c("Gene_Name", xlab, "Donor")) |>
        summarize(tot_ct = sum(Counts), Age=Age) |> distinct() |>
        ggline(x="Age", y="tot_ct", add=c("mean_se", "jitter"),
               linetype="Gene_Name", facet.by=xlab, ylab="Counts", 
               scales="free_y", add.params=list(alpha=0.8), 
               panel.labs.font=list(face='bold'),
               ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=cols)
    line2 <- df2 |> group_by_at(c("Gene_Name", xlab, "Donor")) |>
        summarize(mean_exp = mean(`Normalized Expression`, na.rm=TRUE),
                  Age=Age) |> distinct() |>
        ggline(x="Age", y="mean_exp", add=c("mean_se", "jitter"),
               linetype="Gene_Name",facet.by=xlab, ylab="Normalized Expression",
               add.params=list(alpha=0.8), panel.labs.font=list(face='bold'),
               ggtheme=theme_pubr(base_size=15, border=TRUE), 
               ncol=cols, scales="free_y")
    save_ggplots(paste0(outdir, "/line_graph_counts_", tolower(xlab)), 
                 line1, w, 12)
    save_ggplots(paste0(outdir, "/line_graph_expr_", tolower(xlab)), 
                 line2, w, 12)
}

generate_lineplots <- function(df1, df2, outdir){
    plot_line_age(df1, df2, outdir)
    plot_line_variable(df1, df2, "Celltype", 5, 15, outdir)
    plot_line_variable(df1, df2, "Lineage", 1, 3, outdir)
}

                                        # Statistical modeling
statistics_byGene <- function(df, xlab){
    model = paste("`Normalized Expression` ~", xlab)
    est_fit <- df |> group_by(Gene_Name, Age) |>
        do(fitEST=broom::tidy(lm(model, data=.))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf),
               term = gsub(xlab, "", term))
    ##print(est_fit |> dplyr::count(Gene_Name, Age, p.bonf.cat))
    return(est_fit)
}

statistics_byGene_byPatient <- function(df, xlab){
    model = paste("mean_var ~", xlab)
    est_fit <- df |> group_by_at(c("Donor", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE), Age=Age) |>
        group_by(Gene_Name, Age) |>
        do(fitEST=broom::tidy(lm(model, data=.))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf),
               term = gsub(xlab, "", term))
    ##print(est_fit |> dplyr::count(Gene_Name, Age, p.bonf.cat))
    return(est_fit)
}

statistics_byGene_byPatient_acrossAge <- function(df, xlab){
    model = paste("mean_var ~ Age")
    est_fit <- df |> group_by_at(c("Donor", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE), Age=Age) |>
        group_by_at(c("Gene_Name", xlab)) |>
        do(fitEST=broom::tidy(lm(model, data=.))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf),
               term = gsub(xlab, "", term))
    return(est_fit)
}

statistics_within <- function(df, xlab){
    model = "`Normalized Expression` ~ Gene_Name"
    est_fit <- df |> group_by_at(c(xlab, "Age")) |>
        do(fitEST=broom::tidy(lm(model, data=.))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf),
               term = gsub("Gene_Name", "", term)) |>
        data.table::setnames(old=xlab, new="Feature") |>
        mutate(Feature=gsub(" ", "_", Feature))
    ##print(est_fit |> dplyr::count(Feature, Age, p.bonf.cat))
    return(est_fit)
}

statistics_within_byPatient <- function(df, xlab){
    model = "mean_var ~ Gene_Name"
    est_fit <- df |> group_by_at(c("Donor", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE), Age=Age) |>
        group_by_at(c(xlab, "Age")) |>
        do(fitEST=broom::tidy(lm(model, data=.))) |>
        tidyr::unnest(fitEST) |> filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf,
                                breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf),
               term = gsub("Gene_Name", "", term)) |>
        data.table::setnames(old=xlab, new="Feature") |>
        mutate(Feature=gsub(" ", "_", Feature))
    ##print(est_fit |> dplyr::count(Feature, Age, p.bonf.cat))
    return(est_fit)
}

run_statistical_model <- function(df, label, fnc){
    outfile = paste0("linear_model_comparison_byAge_", label, ".tsv")
    datalist = list()
    for(label in c("Lineage", "Celltype")){
        datalist[[label]] <- fnc(df, label)
    }
    bind_rows(datalist) |> data.table::fwrite(outfile, sep='\t')
}

### Main script section
                                        # Plot venn diagrams by age group
for(age in c("31wk", "3yr", "31yr")){ venn_diagrams(age) }
                                        # All cells
plot_dotplot()
                                        # Filtered cells
outdir = "filter_cells"
dir.create(outdir)
df1 <- prepare_counts()
df2 <- prepare_filtered_data()
generate_lineplots(df1, df2, outdir)
                                        # Statistical model
run_statistical_model(df2, "byGene", statistics_byGene)
run_statistical_model(df2, "byGene_byPatient", statistics_byGene_byPatient)

run_statistical_model(df2, "within", statistics_within)
run_statistical_model(df2, "within_byPatient", statistics_within_byPatient)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
