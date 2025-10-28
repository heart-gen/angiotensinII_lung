#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(memoise)
  library(ggfittext)
  library(ggpubr)
  library(scales)
  library(sessioninfo)
})

#----------------------------#
# Fisher's exact test helpers
#----------------------------#
fisher_exact_generic <- function(df, group_col, in_group_label, in_group_name,
                                 out_group_name) {
                                        # Compute mean expression per patient x group
    dx <- df |> 
        group_by(.data$Donor, .data[[group_col]]) |>
        summarize(
            mean_expr = mean(`Normalized Expression`, na.rm = TRUE),
            .groups = "drop"
        )

                                        # Present / not present counts
    group_yes <- df |>
        filter(.data[[group_col]] == in_group_label, `Normalized Expression` > 0) |>
        nrow()

    group_no <- df |>
        filter(.data[[group_col]] == in_group_label, `Normalized Expression` == 0) |>
        nrow()

    other_yes <- df |>
        filter(.data[[group_col]] != in_group_label, `Normalized Expression` > 0) |>
        nrow()

    other_no <- df |>
        filter(.data[[group_col]] != in_group_label, `Normalized Expression` == 0) |>
        nrow()

    contingency <- data.frame(
        Present      = c(group_yes, other_yes),
        `Not Present` = c(group_no, other_no),
        row.names    = c(in_group_name, out_group_name)
    )
    return(fisher.test(contingency))
}

fisher_exact_compartment <- function(df, xlab) {
    fisher_exact_generic(
        df = df, group_col = "Lineage", in_group_label = xlab,
        in_group_name = "In Compartment", out_group_name = "Not in Compartment"
    )
}

fisher_exact_annotation <- function(df, xlab) {
    fisher_exact_generic(
        df = df, group_col = "Celltype", in_group_label = xlab,
        in_group_name = "In Cell type", out_group_name = "Not in Cell type"
    )
}

#--------------------------------------#
# Run enrichment for a set of locations
#--------------------------------------#

enrichment_loop <- function(df, fisher_fun, locations, annotation_label,
                            genes = c("AGTR1", "AGTR2")) {
    results <- lapply(genes, function(gene) {
        gene_df <- df |> filter(Gene_Name == gene)

        loc_stats <- lapply(locations, function(loc) {
            ft <- fisher_fun(gene_df, loc)

            data.frame(
                Annotation       = annotation_label,
                Annotation_Name  = loc,
                Gene             = gene,
                OR               = unname(ft$estimate),
                P                = ft$p.value,
                stringsAsFactors = FALSE
            )
        })

        bind_rows(loc_stats)
    })


    return(bind_rows(results) |> mutate(FDR = p.adjust(P, method = "bonferroni")))
}

#----------------------#
# Full enrichment stage
#----------------------#
run_enrichment <- function(expr_path) {
    df <- data.table::fread(expr_path)

                                        # Compartment
    comp_locs <- unique(df$Lineage)
    dat_comp <- enrichment_loop(
        df = df,
        fisher_fun = fisher_exact_compartment,
        locations = comp_locs,
        annotation_label = "Compartment"
    )

                                        # Cell type
    celltype_locs <- unique(df$Celltype)
    dat_celltype <- enrichment_loop(
        df = df,
        fisher_fun = fisher_exact_annotation,
        locations = celltype_locs,
        annotation_label = "Cell type"
    )

    bind_rows(dat_comp, dat_celltype)
}

#----------------------------#
# Plotting
#----------------------------#
save_plot <- function(p, file_stub, w, h) {
    for (ext in c(".pdf", ".png")) {
        ggsave(
            filename = paste0(file_stub, ext), plot = p, width = w, height = h
        )
    }
}

load_enrichment <- function(tsv_path) {
    return(data.table::fread(tsv_path))
}
mem_enrich <- memoise::memoise(load_enrichment)

prepare_plot_data <- function(tsv_path) {
    err  <- 1e-7
    err2 <- 1e-100

    mem_enrich(tsv_path) |>
        mutate_if(is.character, as.factor) |>
        mutate(
            `-log10(FDR)`   = ifelse(FDR != 0, -log10(FDR), -log10(err)),
            `OR Percentile` = OR / (1 + OR),
            p_fdr_sig       = FDR < 0.05,
            `log2(OR)`      = log2(OR + err2),
            p_fdr_cat       = cut(
                FDR,
                breaks = c(1, 0.05, 0.01, 0.005, 0),
                labels = c("<= 0.005", "<= 0.01", "<= 0.05", "> 0.05"),
                include.lowest = TRUE
            )
        )
}
mem_plot_df <- memoise::memoise(prepare_plot_data)

plot_tile <- function(tsv_path, annotation_label, w, h,
                      out_stub_prefix = "tileplot_enrichment_") {
    dt <- mem_plot_df(tsv_path) |>
        filter(Annotation == annotation_label)

    y0 <- min(dt$`log2(OR)`, na.rm = TRUE) - 0.1
    y1 <- max(dt$`log2(OR)`, na.rm = TRUE) + 0.1

    p <- ggplot(
        dt,
        aes(x = Gene, y = Annotation_Name, fill = `log2(OR)`,
            label = ifelse(
                p_fdr_sig,
                format(round(`-log10(FDR)`, 1), nsmall = 1),
                "")
            )) +
        geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(
            colors = c("blue", "white", "red"),
            values = scales::rescale(c(y0, 0, y1)),
            limits = c(y0, y1)) +
        labs(x = NULL, y = NULL) +
        ggpubr::theme_pubr(base_size = 20, border = FALSE) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            axis.title = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold"),
            strip.text = element_text(face = "bold")
        )

    outfile <- paste0(out_stub_prefix,
                      gsub(" ", "_", tolower(annotation_label)))

    save_plot(p, outfile, w, h)
}

#----------------------------#
# Main pipeline
#----------------------------#
main <- function() {
    expression_file <- "../_m/normalized_expression.txt.gz"
    enrichment_tsv  <- "cell_annotation_enrichment_analysis.tsv"

                                        # Enrichment analysis + write TSV
    enrichment_results <- run_enrichment(expression_file)
    data.table::fwrite(enrichment_results, enrichment_tsv, sep = "\t")

                                        # Plots from TSV
    plot_tile(tsv_path = enrichment_tsv, annotation_label = "Compartment",
              w = 5.4, h = 4)
    plot_tile(tsv_path = enrichment_tsv, annotation_label = "Cell type",
              w = 8, h = 13)

                                        # Reproducibility info
    cat("Reproducibility information:\n")
    print(Sys.time())
    print(proc.time())
    options(width = 120)
    print(sessioninfo::session_info())
}

                                        # Execute when run as script
main()
