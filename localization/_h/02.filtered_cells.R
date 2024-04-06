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
    sce <- as.SingleCellExperiment(readRDS(fn))
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
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_cell_annotation"),
                 pp, w, h)
    rr = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="compartment") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_compartment"),
                 rr, 4, 5)
    qq = DotPlot(object = df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="subclusters") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_subcluster"),
                 qq, w, h+2.5)
    gg = DotPlot(object=df.seurat, features = c("AGTR1", "AGTR2"),
                 group.by="clusters") + RotatedAxis()
    save_ggplots(paste0(outdir,"/dotplot_angiotensinII_clusters"),
                 gg, w, h)
}

prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df <- logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Patient=colData(angiotensin2)$patient,
               Cell_Annotation=colData(angiotensin2)$cell_type,
               Subcluster=colData(angiotensin2)$subclusters,
               Cluster=colData(angiotensin2)$clusters,
               Compartment=colData(angiotensin2)$compartment) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
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

plot_box_stats <- function(dt, xlab, outdir, mycomps){
    outfile = paste0(outdir,"/normalized_expression_boxplot_", xlab)
    bar = ggboxplot(dt, x=xlab, y="Normalized Expression", fill=xlab,
                    facet.by="Gene_Name", add="jitter", legend="bottom",
                    outlier.shape=NA, xlab="",
                    panel.labs.font=list(face='bold'),
                    add.params=list(alpha=0.5),
                    ggtheme=theme_pubclean(base_size=15),
                    ncol=2) + rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bar, 7, 6)
}

plot_box_noStats <- function(dt, xlab, outdir, r, w, h){
    outfile = paste0(outdir,"/normalized_expression_boxplot_noStat_", xlab)
    bar = dt |> mutate_if(is.character, as.factor) |>
        ggboxplot(x=xlab, y="Normalized Expression", fill=xlab,
                  facet.by="Gene_Name", add="jitter", legend="bottom",
                  outlier.shape=NA, xlab="",
                  panel.labs.font=list(face='bold'),
                  add.params=list(alpha=0.5),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=2) +
        rotate_x_text(r)
    save_ggplots(tolower(outfile), bar, w, h)
}

plot_box_byPatient <- function(dt, xlab, outdir, r, w, h){
    outfile = paste0(outdir,"/normalized_expression_boxplot_byPatient_", xlab)
    bxp = dt |> mutate_if(is.character, as.factor) |>
        group_by_at(c("Patient", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        ggboxplot(x=xlab, y="mean_var", fill=xlab, facet.by="Gene_Name",
                  add="jitter", legend="bottom", outlier.shape=NA, xlab="",
                  ylab="Normalized Expression", #palette="npg",
                  panel.labs.font=list(face='bold'),
                  add.params=list(alpha=0.8),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), ncol=2) +
        rotate_x_text(r)
    save_ggplots(tolower(outfile), bxp, w, h)
}

generate_boxplots <- function(df, outdir){
    compartment <- list(
        c("Endothelial", "Epithelial"), c("Epithelial", "Immune"),
        c("Stroma", "Immune"), c("Endothelial", "Immune"),
        c("Stroma", "Epithelial"), c("Endothelial", "Stroma")
    )
                                        # No statistics
    plot_box_noStats(df, "Compartment", outdir, 45, 6, 5)
    plot_box_noStats(df, "Cell_Annotation", outdir, 90, 9, 7)
    plot_box_noStats(df, "Subcluster", outdir, 90, 9, 7)
                                        # by Patient
    plot_box_byPatient(df, "Compartment", outdir, 45, 6, 5)
    plot_box_byPatient(df, "Cell_Annotation", outdir, 90, 9, 7)
    plot_box_byPatient(df, "Subcluster", outdir, 90, 9, 7)
}

statistics_byGene <- function(df, xlab){
    model = paste("`Normalized Expression` ~", xlab)
    est_fit <- df |> group_by(Gene_Name) |>
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
    #print(est_fit)
    return(est_fit)
}

statistics_byGene_byPatient <- function(df, xlab){
    model = paste("mean_var ~", xlab)
    est_fit <- df |> group_by_at(c("Patient", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        group_by(Gene_Name) |>
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
    #print(est_fit)
    return(est_fit)
}

statistics_within <- function(df, xlab){
    model = "`Normalized Expression` ~ Gene_Name"
    est_fit <- df |> group_by_at(xlab) |>
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
    #print(est_fit)
    return(est_fit)
}

statistics_within_byPatient <- function(df, xlab){
    model = "`mean_var` ~ Gene_Name"
    est_fit <- df |> group_by_at(c("Patient", xlab, "Gene_Name")) |>
        summarize(mean_var = mean(`Normalized Expression`, na.rm=TRUE)) |>
        group_by_at(xlab) |>
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
    #print(est_fit)
    return(est_fit)
}

run_statistical_model <- function(df, label, fnc){
    outfile = paste0("linear_model_comparison_", label, ".tsv")
    datalist = list()
    for(label in c("Compartment", "Cell_Annotation", "Subcluster")){
        datalist[[label]] <- fnc(df, label)
    }
    bind_rows(datalist) |> data.table::fwrite(outfile, sep='\t')
}

### Main script section
                                        # Filtered cells
outdir = "filter_cells"
dir.create(outdir)
plot_dotplot(outdir, TRUE)
df <- prepare_filtered_data()
calculate_proportion(df, outdir, TRUE)
generate_boxplots(df, outdir)
                                        # Statistical model
run_statistical_model(df, "byGene", statistics_byGene)
run_statistical_model(df, "byGene_byPatient", statistics_byGene_byPatient)
run_statistical_model(df, "within", statistics_within)
run_statistical_model(df, "within_byPatient", statistics_within_byPatient)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
