## This script examines differences in expression
suppressPackageStartupMessages({
    library("here")
    library("dplyr")
    library("ggpubr")
    library("Seurat")
    library("speckle")
    library("viridisLite")
    library("DeconvoBuddies")
    library("SingleCellExperiment")
})

local_setup <- function() {
    system_readline <- "/usr/lib/libreadline.so.8"
    Sys.setenv(LD_PRELOAD = system_readline)
}

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_data <- function() {
    input_file <- "pericyte.hlca_core.subclustered.h5ad"
    sce  <- zellkonverter::readH5AD(input_file)
    rowData(sce)$ensembl_id <- rownames(sce)
    rownames(sce) <- rowData(sce)$feature_name
    return(sce)
}

plot_overlap <- function(sce) {
    seur <- as.Seurat(sce, counts = "counts", data = "logcounts")
    seur$leiden       <- forcats::fct_rev(factor(seur$leiden))
    seur$active.ident <- seur$leiden
    Idents(seur) <- "leiden"
    v_cols <- viridis(10)
    col_1  <- v_cols[1]; col_2 <- v_cols[7]
    pp <- FeaturePlot(seur, features = c("ACTA2", "AGTR1"), blend = TRUE,
                      reduction = "X_umap", pt.size = 0.1,
                      cols = c("lightgrey", col_1, col_2))
    ggplot2::ggsave("colocalization_feature_plot.pdf", plot = pp,
                    width = 12, height = 4)
    ggplot2::ggsave("colocalization_feature_plot.png", plot = pp,
                    width = 12, height = 4, dpi = 300)
}

generate_marker_genes <- function(sce) {
    ratios   <- get_mean_ratio(sce, cellType_col = "leiden",
                               assay_name = "logcounts",
                               gene_ensembl = "ensembl_id",
                               gene_name = "feature_name")
    write.csv(ratios, file = "marker_stats_genes.csv",
              row.names = FALSE)
                                        # Plot the top ten markers
    plot_marker_express_ALL(sce, ratios, n_genes = 2,
                            pdf_fn = "top2-marker-genes.pdf",
                            cellType_col = "leiden",
                            color_pal = NULL, plot_points = FALSE)
    return(ratios)
}

props_analysis <- function(sce, gene_name, transform = "asin",
                           DIR = "cell_props") {
                                        # Generate mask
    expr_matrix <- counts(sce)[gene_name, , drop = FALSE]
    is_FEATURE  <- as.integer(expr_matrix > 0)
    colData(sce)$is_FEATURE <- as.vector(is_FEATURE)
                                        # Calculate change in proportion
    prop_df <- propeller(clusters = sce$leiden, sample = sce$donor_id,
                         group = sce$is_FEATURE, transform = transform)
    print(prop_df)
    out_file1 <- paste0(DIR, "/", tolower(gene_name),
                        ".propeller_proportion.t-test.tsv")
    data.table::fwrite(prop_df, file=out_file1,
                       sep="\t", row.names=TRUE)

                                        # Transform proportions
    pheno_df <- colData(sce) |> as.data.frame() |>
        select("donor_id", "is_FEATURE") |> distinct()
    props    <- getTransformedProps(sce$leiden, sce$donor_id,
                                    transform = transform)
    df  <- as.data.frame(props$TransformedProps) |>
        inner_join(pheno_df, by=c("sample"="donor_id"),
                   relationship = "many-to-many") |>
        tidyr::drop_na()
    df$clusters <- factor(df$clusters)
    df$is_FEATURE <- factor(df$is_FEATURE)
    levels(df$is_FEATURE) <- c("Absent", "Present")
                                        # Plot proportions
    bar <- ggbarplot(df, x="clusters", y="Freq", fill="is_FEATURE",
                     add="mean_se", palette="npg",
                     position=position_dodge(),
                     xlab="Subclusters", ylab="Arcsin(Proportion)",
                     ggtheme=theme_pubr(base_size=15)) +
        rotate_x_text(45)
    out_file2 <- paste0(DIR, "/", tolower(gene_name), ".barplot")
    save_plot(bar, out_file2, 6, 6)
}
#### Main analysis ####
                                        # For local run
local_setup()

                                        # Load data
sce <- load_data()

                                        # Plot overlap
plot_overlap()

                                        # Marker genes
ratios <- generate_marker_genes(sce)

                                        # Cell proportions
dir.create("cell_props")
for(gene in c("AGTR1", "ACTA2")) {
    props_analysis(sce, gene, DIR="cell_props")
}

#### Reproducibility information ####
cat("Reproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
reticulate::py_config()
sessioninfo::session_info()
