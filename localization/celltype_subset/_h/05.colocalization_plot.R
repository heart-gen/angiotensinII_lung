## This script examines differences in expression
suppressPackageStartupMessages({
    library("here")
    library("dplyr")
    library("Seurat")
    library("viridisLite")
    library("DeconvoBuddies")
    library("SingleCellExperiment")
})

local_setup <- function() {
    system_readline <- "/usr/lib/libreadline.so.8"
    Sys.setenv(LD_PRELOAD = system_readline)
}

load_data <- function() {
    input_file <- "pericyte.hlca_core.subclustered.h5ad"
    sce  <- zellkonverter::readH5AD(input_file)
    rowData(sce)$ensembl_id <- rownames(sce)
    rownames(sce) <- rowData(sce)$feature_name
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
#### Main analysis ####
                                        # For local run
local_setup()

                                        # Load data
sce <- load_data()

                                        # Plot overlap
plot_overlap()

#### Reproducibility information ####
cat("Reproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
reticulate::py_config()
sessioninfo::session_info()
