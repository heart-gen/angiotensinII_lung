#' Get Mean Ratio for Each Gene x Cell Type

library(dplyr)
library(purrr)
library(matrixStats)

get_mean_ratio2 <- function(sce, cellType_col = "cellType", assay_name = "logcounts", add_symbol = TRUE) {
    # RCMD fix
    cellType.target <- NULL
    cellType <- NULL
    ratio <- NULL

    cell_types <- unique(sce[[cellType_col]])
    names(cell_types) <- cell_types

    sce_assay <- as.matrix(SummarizedExperiment::assays(sce)[[assay_name]])

    ## Get mean expression for each gene for each cellType
    cell_means <- map(cell_types, ~ as.data.frame(base::rowMeans(sce_assay[, sce[[cellType_col]] == .x])))

    cell_means <- do.call("rbind", cell_means)
    colnames(cell_means) <- "mean"
    ## Define columns
    cell_means$cellType <- rep(cell_types, each = nrow(sce))
    cell_means$gene <- rep(rownames(sce), length(cell_types))
    # print(head(cell_means))

    ## Filter and calculate ratio for each celltype
    ratio_tables <- map(cell_types, ~ .get_ratio_table(.x, sce, sce_assay, cellType_col, cell_means))

    ratio_tables <- do.call("rbind", ratio_tables)

    # max_digits <- nchar(max(ratio_tables$ratio_tables))

    if (add_symbol) {
        ratio_tables$Symbol <- SummarizedExperiment::rowData(sce)[ratio_tables$gene, ]$Symbol
        # ratio_tables <- ratio_tables %>%
        #    mutate(ratio_anno = paste0(stringr::str_pad(rank_ratio, max_digits, "left"),": ",Symbol))
    }
    ratio_tables <- ratio_tables %>%
        mutate(anno_ratio = paste0(cellType.target, "/", cellType, ": ", base::round(ratio, 3)))

    return(ratio_tables)
}


.get_ratio_table <- function(x, sce, sce_assay, cellType_col, cell_means) {
    # RCMD Fix
    mean.target <- NULL
    gene <- NULL
    ratio <- NULL
    cellType.target <- NULL
    cellType <- NULL

    # filter target median != 0
    median_index <- matrixStats::rowMedians(sce_assay[, sce[[cellType_col]] == x]) != 0
    # message("Median == 0: ", sum(!median_index))
    # filter for target means
    target_mean <- cell_means[cell_means$cellType == x, ]
    target_mean <- target_mean[median_index, ]
    colnames(target_mean) <- c("mean.target", "cellType.target", "gene")

    nontarget_mean <- cell_means[cell_means$cellType != x, ]

    ratio_table <- dplyr::left_join(target_mean, nontarget_mean, by = "gene") %>%
        mutate(ratio = mean.target / mean) %>%
        dplyr::group_by(gene) %>%
        arrange(ratio) %>%
        dplyr::slice(1) %>%
        dplyr::select(gene, cellType.target, mean.target, cellType, mean, ratio) %>%
        arrange(-ratio) %>%
        dplyr::ungroup() %>%
        mutate(rank_ratio = dplyr::row_number())

    return(ratio_table)
}
