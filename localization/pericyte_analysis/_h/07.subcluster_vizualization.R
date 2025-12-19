## Analyze pericyte subclusters

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(SingleCellExperiment)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_data <- function(){
    sce <- zellkonverter::readH5AD("results/airspace_pericyte_states.h5ad")
    return(sce)
}

prepare_data <- function(sce){
    cols_to_keep <- c(
        "donor_id", "leiden_pericytes", "airspace_score",
        "capillary_sig_z", "arteriolar_sig_z", "airspace_dist_z",
        "pericyte_combined_score", "pericyte_state", "age",
        "sex", "disease", "self_reported_ethnicity"
    )
    sub_sce <- sce[c("AGTR1"),]

    expr_long <- Matrix::summary(logcounts(sub_sce)) |>
        transmute(
            Gene_Name = rownames(sub_sce)[i],
            cell_id   = colnames(sub_sce)[j],
            Normalized_Expression = x
        )

    cell_meta <- as.data.frame(colData(sub_sce)) |>
        filter(subclusters == "Pericytes") |>
        mutate(age = age_or_mean_of_age_range) |>
        select(all_of(cols_to_keep)) |>
        tibble::rownames_to_column("cell_id")

    df <- inner_join(expr_long, cell_meta, by="cell_id")
    donor_df <- group_by(df, donor_id) |>
        summarise(
            Normalized_Expression = mean(Normalized_Expression, na.rm = TRUE),
            .groups = "drop"
        )
    return(df)
}

boxplot_subcluster <- function(df, label="Normalized_Expression") {
    bxp <- ggboxplot(df, x="leiden_pericytes", y=label)
}

correlation_subcluster <- function(df, label) {
    sca <- ggscatter(df, x="Normalized_Expression", y=label,
                     facet.by="leiden_pericytes")
}

### Main script section
sce <- load_data()
print(sce)

df  <- prepare_data(sce)
print(dim(df))
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
