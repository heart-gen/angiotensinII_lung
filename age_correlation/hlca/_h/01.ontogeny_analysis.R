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
    fn  <- here::here("inputs/hlca/_m/hlca_core.h5ad")
    sce <- zellkonverter::readH5AD(fn)
    ## Generate RSE object
    names(assays(sce)) <- "counts"
    sce <- scuttle::logNormCounts(sce)
    sce <- scran::computeSumFactors(sce)
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_coarse_for_GWAS_and_modeling
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$patient     <- sce$donor_id
    colData(sce)$age         <- sce$age_or_mean_of_age_range
    colLabels(sce) <- sce$cell_type
    return(sce)
}
memSC <- memoise::memoise(load_data)

add_qc <- function(sce){
    is_mito <- grep("MT-", rowData(sce)$feature_name)
    sce     <- scuttle::addPerCellQCMetrics(sce, subsets=list(mito=is_mito))
    sce     <- scuttle::addPerFeatureQCMetrics(sce)
    return(sce)
}

filter_qc <- function(sce){
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
                                        # Select angiotensin II genes
    rownames(sce) <- rowData(sce)[, 2]
    genes = c("AGTR1", "AGTR2")
    angiotensin2 <- sce[genes,]
    return(angiotensin2)
}
memAGTR <- memoise::memoise(extract_angiotensinII)

                                        # Prepare data functions
prepare_filtered_data <- function(){
    angiotensin2 <- memAGTR()
    df = logcounts(angiotensin2) |> t() |> as.data.frame() |>
        mutate(Donor=colData(angiotensin2)$patient,
               Celltype=colData(angiotensin2)$cell_type,
               Compartment=colData(angiotensin2)$compartment,
               Age=colData(angiotensin2)$age) |>
        filter(AGTR1 > 0 | AGTR2 > 0) |> droplevels() |>
        tidyr::pivot_longer(-c(Donor, Celltype, Compartment, Age),
                            names_to="Gene_Name",
                            values_to="Normalized Expression") |>
        mutate_if(is.character, as.factor)
    return(df)
}

                                        # Plotting functions
plot_age <- function(df2, outdir){
    sca2 <- df2 |> group_by_at(c("Gene_Name", "Donor")) |>
        reframe(mean_exp = mean(`Normalized Expression`, na.rm=TRUE),
                  Age=Age) |> distinct() |>
        ggscatter(x="Age", y="mean_exp", facet.by="Gene_Name",
                  panel.labs.font=list(face='bold'),
                  add="reg.line", conf.int=TRUE, cor.coef=TRUE,
                  ylab="Normalized Expression",
                  add.params=list(color="blue", fill="lightgray"),
                  ggtheme=theme_pubr(base_size=15, border=TRUE))
    save_ggplots(paste0(outdir, "/scatter_expr"), sca2, 6, 4)
}

plot_age_variable <- function(df2, gene, xlab, cols, w, outdir){
    sca2 <- df2 |> group_by_at(c("Gene_Name", xlab, "Donor")) |>
        reframe(mean_exp = mean(`Normalized Expression`, na.rm=TRUE),
                Age=Age) |> distinct() |> filter("Gene_Name" == gene) |>
        ggscatter(x="Age", y="mean_exp", facet.by=xlab, add="reg.line", 
                  conf.int=TRUE, cor.coef=TRUE, ylab="Normalized Expression", 
                  panel.labs.font=list(face='bold'),
                  add.params=list(color="blue", fill="lightgray"),
                  ggtheme=theme_pubr(base_size=15, border=TRUE), 
                  ncol=cols, scales="free_y")
    save_ggplots(paste0(outdir, "/scatter_expr.", gene, ".", 
                        tolower(xlab)), sca2, w, 12)
}

generate_plots <- function(df2, outdir){
    plot_age(df2, outdir)
    for(gene in c("AGTR1", "AGTR2")){
        plot_age_variable(df2, gene, "Celltype", 5, 15, outdir)
        plot_age_variable(df2, gene, "Compartment", 1, 3, outdir)
    }
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
