## This script examines differences in expression
suppressPackageStartupMessages({
    library("here")
    library("dplyr")
    library("scran")
    library("edgeR")
    library("ggpubr")
    library("Seurat")
    library("speckle")
    library("scuttle")
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

load_data <- function(model_type = c("core", "full")) {
    model_type <- match.arg(model_type)
    input_file <- switch(
        model_type,
        "core" = "pericyte.hlca_core.subclustered.h5ad",
        "full" = "pericyte.hlca_full.subclustered.h5ad"
    )
    message("Loading model: ", model_type, " (", input_file, ")")
    sce  <- zellkonverter::readH5AD(input_file)
    rowData(sce)$ensembl_id <- rownames(sce)
    rownames(sce) <- rowData(sce)$feature_name
    return(sce)
}

compute_cluster_avg <- function(sce, cluster_key = "cell_types") {
    clusters <- colData(sce)[[cluster_key]]
    avg <- aggregateAcrossCells(sce, ids = clusters)
    avg <- logNormCounts(avg)
    cluster_avg_expr <- as.data.frame(assay(avg, "logcounts"))
    colnames(cluster_avg_expr) <- levels(factor(clusters))
    return(cluster_avg_expr)
}

calc_specificity_ratio <- function(cluster_avg_expr_df) {
    specificity <- as.data.frame(matrix(NA, 
                                        nrow = nrow(cluster_avg_expr_df),
                                        ncol = ncol(cluster_avg_expr_df)))
    rownames(specificity) <- rownames(cluster_avg_expr_df)
    colnames(specificity) <- colnames(cluster_avg_expr_df)
  
    for (cluster in colnames(cluster_avg_expr_df)) {
        cluster_expr <- cluster_avg_expr_df[[cluster]]
        other_expr <- rowSums(
            cluster_avg_expr_df[ , setdiff(colnames(cluster_avg_expr_df),
                                           cluster), drop=FALSE]
        )
        specificity[[cluster]] <- cluster_expr / (other_expr + 1e-9)
    }
    colnames(specificity) <- paste0("cluster_", colnames(specificity))
    return(specificity)
}

plot_overlap <- function(sce, model_tag, dirname = "figures") {
    outdir   <- file.path(dirname, model_tag)
    dir.create(outdir, showWarnings=FALSE)
    seur <- as.Seurat(sce, counts = "counts", data = "logcounts")
    seur$leiden       <- forcats::fct_rev(factor(seur$leiden))
    seur$active.ident <- seur$leiden
    Idents(seur) <- "leiden"

    v_cols <- viridis(10)
    col_1  <- v_cols[1]; col_2 <- v_cols[7]

    pp <- FeaturePlot(
        seur, features = c("ACTA2", "AGTR1"), blend = TRUE,
        reduction = "X_umap", pt.size = 0.1,
        cols = c("lightgrey", col_1, col_2)
    )
    ggplot2::ggsave(file.path(outdir, "colocalization_feature_plot.pdf"),
                    plot = pp, width = 12, height = 4)
    ggplot2::ggsave(file.path(outdir, "colocalization_feature_plot.png"),
                    plot = pp, width = 12, height = 4, dpi = 300)
}

generate_marker_genes <- function(sce, model_tag) {
    outdir   <- file.path("marker_genes", model_tag)
    ratios   <- get_mean_ratio(sce, cellType_col = "leiden",
                               assay_name = "logcounts",
                               gene_ensembl = "ensembl_id",
                               gene_name = "feature_name")
    write.csv(ratios,
              file = file.path(outdir, "marker_stats_genes.csv"),
              row.names = FALSE)
                                        # Plot the top ten markers
    outfile <- file.path(outdir, "top2-marker-genes.pdf")
    plot_marker_express_ALL(
        sce, ratios, n_genes = 2, pdf_fn = outfile,
        cellType_col = "leiden", color_pal = NULL,
        plot_points = FALSE
    )
    return(ratios)
}

run_de_analysis <- function(sce, cluster_key = "leiden") {
    clusters <- factor(colData(sce)[[cluster_key]])
    design <- model.matrix(~0 + clusters)
    colnames(design) <- paste0("cluster_", levels(clusters))
    
    y <- edgeR::DGEList(counts = counts(sce), genes = rowData(sce))
    y <- edgeR::calcNormFactors(y)
    
    v <- limma::voom(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
  
                                        # Make contrasts: each cluster vs rest
    results_list <- list()
    cluster_names <- paste0("cluster_", levels(clusters))
    for (cl in cluster_names) {
        contrast_vec <- rep(-1 / (length(cluster_names) - 1),
                            length(cluster_names))
        names(contrast_vec) <- cluster_names
        contrast_vec[cl]    <- 1
        
        fit2 <- limma::contrasts.fit(fit, contrasts = contrast_vec)
        fit2 <- limma::eBayes(fit2)
        tab  <- limma::topTable(fit2, number = Inf, sort.by = "P")
        results_list[[cl]] <- tab
    }
    return(results_list)
}

props_analysis <- function(sce, gene_name, model_tag,
                           transform = "asin", dirname = "cell_props") {
    outdir <- file.path(dirname, model_tag)
    dir.create(outdir, showWarnings = FALSE)

                                        # Generate mask
    expr_matrix <- counts(sce)[gene_name, , drop = FALSE]
    is_FEATURE  <- as.integer(expr_matrix > 0)
    colData(sce)$is_FEATURE <- as.vector(is_FEATURE)

                                        # Calculate change in proportion
    prop_df <- propeller(
        clusters = sce$leiden, sample = sce$donor_id,
        group = sce$is_FEATURE, transform = transform
    )
    print(prop_df)
    out_file1 <- file.path(
        outdir, paste0(tolower(gene_name),
                        ".propeller_proportion.t-test.tsv")
    )
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
    bar <- ggbarplot(
        df, x="clusters", y="Freq", fill="is_FEATURE",
        add="mean_se", palette="npg", position=position_dodge(),
        xlab="Subclusters", ylab="Arcsin(Proportion)",
        ggtheme=theme_pubr(base_size=15)
    ) + rotate_x_text(45)
    out_file2 <- file.path(
        outdir, paste0(tolower(gene_name), ".barplot")
    )
    save_plot(bar, out_file2, 6, 6)
    
}

select_markers <- function(cluster_avg_expr_df, specificity_ratio,
                           de_results, specificity_threshold = 4.0,
                           expression_threshold = 1.0, fdr_threshold = 0.05) {
    markers <- list()
    for (cl in names(de_results)) {
        de_df <- de_results[[cl]]
        common_genes <- intersect(rownames(de_df),
                                  rownames(cluster_avg_expr_df))
        common_genes <- intersect(common_genes,
                                  rownames(specificity_ratio))
        if (length(common_genes) == 0) next
        
        spec <- specificity_ratio[common_genes, cl]
        avg <- cluster_avg_expr_df[common_genes, cl]
        de_sub <- de_df[common_genes, , drop = FALSE]
    
        candidates <- common_genes[
            spec > specificity_threshold &
            avg > expression_threshold &
            de_sub$adj.P.Val < fdr_threshold
        ]
    
        if (length(candidates) > 0) {
            markers[[cl]] <- data.frame(
                cluster = cl,
                gene = candidates,
                specificity_ratio = spec[candidates],
                mean_expression = avg[candidates],
                pval_adj = de_sub[candidates, "adj.P.Val"],
                logfc = de_sub[candidates, "logFC"],
                stringsAsFactors = FALSE
            )
        }
    }
    return(dplyr::bind_rows(markers) |>
           arrange(cluster, desc(specificity_ratio)))
}

plot_marker_expression <- function(sce, gene, cluster_key = "leiden",
                                   dirname = "marker_genes", model = "core") {
    outdir <- file.path(dirname, model)
    dir.create(outdir, showWarnings = FALSE)
    df <- data.frame(
        expression = assay(sce, "logcounts")[gene, ],
        cluster = colData(sce)[[cluster_key]]
    )
    p <- ggplot(df, aes(x = cluster, y = expression)) +
        geom_violin(trim = TRUE, fill = "lightblue") +
        geom_boxplot(width = 0.1, outlier.shape = NA) +
        ggtitle(paste("Expression of", gene, "across clusters (logcounts)")) +
        theme_bw()
  
    ggsave(file.path(outdir, paste0(tolower(gene), "_expression.png")),
           p, width = 6, height = 4, dpi = 300)
    ggsave(file.path(outdir, paste0(tolower(gene), "_expression.pdf")),
           p, width = 6, height = 4)
}

calculate_specificity_markers <- function(sce, model, cluster_key,
                                          dirname = "marker_genes") {
    outdir <- file.path(dirname, model)
    cluster_avg_expr_df <- compute_cluster_avg(sce, cluster_key)
    specificity_ratio   <- calc_specificity_ratio(cluster_avg_expr_df)
    de_results <- run_de_analysis(sce, cluster_key)
    markers_df <- select_markers(cluster_avg_expr_df, specificity_ratio,
                                 de_results,
                                 specificity_threshold = 1.25)
    write.table(specificity_ratio, file.path(outdir, "specificity_data.tsv"),
                sep = "\t", quote = FALSE)
    write.table(markers_df, file.path(outdir, "candidate_markers_for_qpcr.tsv"),
                sep = "\t", quote = FALSE)
    message("Saved ", nrow(markers_df), " candidate markers.")
}

#### Main analysis ####
for (model in c("core", "full")) {   
                                        # Load data
    sce <- load_data(model)

                                        # Plot overlap
    plot_overlap(sce, model)

                                        # Marker genes
    ratios <- generate_marker_genes(sce, model)
    calculate_specificity_markers(sce, model, "leiden")

                                        # Cell proportions
    for(gene in c("AGTR1", "ACTA2")) {
        props_analysis(sce, gene, model_tag = model,
                       dirname="cell_props")
    }
}

#### Reproducibility information ####
cat("Reproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
reticulate::py_config()
sessioninfo::session_info()
