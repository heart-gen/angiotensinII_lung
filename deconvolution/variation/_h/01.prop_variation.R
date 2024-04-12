################################################################################
## This script examines variance within cell type populations in the GTEx v8  ##
## lung data as well as plots proportions.                                    ##
################################################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(plotly)
    library(ggpubr)
    library(BayesPrism)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

get_cellprop <- function(FULL){
    if(FULL){
        load("../../_m/prism_results.03.lung_GTEx.RDat")
    } else {
        load("../../_m/prism_results.02.lung_GTEx.RDat")
    }
    theta <- get.fraction(bp=bp.res,
                          which.theta="final",
                          state.or.type="type") |>
        as.data.frame() |>
        tibble::rownames_to_column("Sample_ID") |>
        tidyr::pivot_longer(!Sample_ID,
                            names_to="Cell_Type",
                            values_to="Proportion") |>
        janitor::clean_names()
    return(theta)
}

get_compartment <- function(FULL){
    load("../../_m/prism_results.01.lung_GTEx.RDat")
    theta <- get.fraction(bp=bp.res,
                          which.theta="final",
                          state.or.type="type") |>
        as.data.frame() |>
        tibble::rownames_to_column("Sample_ID") |>
        tidyr::pivot_longer(!Sample_ID,
                            names_to="Cell_Type",
                            values_to="Proportion") |>
        janitor::clean_names()
    return(theta)
}

plotNsave_bar <- function(fnc, label, w, h, FULL=FALSE){
    bar <- fnc(FULL) |> group_by(cell_type) |>
        summarize(Proportion=round(mean(proportion), 2)) |>
        mutate(Tissue = "Lung") |>
        ggbarplot(x="Tissue", y="Proportion", fill="cell_type",
                  color="cell_type", palette="npg", label=TRUE, xlab="",
                  ylab="Mean Proportion", lab.pos="in", legend="left",
                  ggtheme=theme_pubr(base_size=15,border=TRUE))
    save_ggplots(paste0(label, "_barplot"), bar, w, h)
}

plotNsave_pie <- function(FULL){
    outfile <- ifelse(FULL, "celltype_piechart.full",
                      "celltype_piechart.sig")
    df   <- get_cellprop(FULL) |> group_by(cell_type) |>
        summarize(Proportion=round(mean(proportion), 3)) |>
        mutate(Tissue = "Lung")
    labs <- paste0(df$cell_type, " (", df$Proportion * 100, "%)")
    fig <- plot_ly(df, labels=~cell_type, values=~Proportion, type="pie",
                   textposition="inside", textinfo="label+percent",
                   insidetextfont=list(color="#FFFFFF"), showlegend=FALSE)
    pie  <- ggpie(df, "Proportion", label=labs, fill="cell_type",
                  color="white", repel=TRUE, legend="",
                  ggtheme=theme_pubr(base_size=15, border=TRUE))
    save_ggplots(outfile, pie, 10, 10)
}

plotNsave_prop <- function(fnc, label, w, h, FULL=FALSE){
    bxp <- ggboxplot(fnc(FULL), x="cell_type", y="proportion", #add="jitter",
                     ylab="Proportion", #fill="cell_type", #outlier.shape=NA,
                     xlab="", panel.labs.font=list(face='bold'), legend="",
                     add.params=list(alpha=0.5),
                     ggtheme=theme_pubr(base_size=15,border=TRUE)) +
        rotate_x_text(45) + font("y.title", face="bold")
    save_ggplots(paste0(label, "_prop"), bxp, w, h)
}

#### Main
                                        # Sample variation
lf <- logr::log_open("summary.log", logdir=FALSE, autolog=FALSE)
logr::log_print("Summary of proportion variation:")
logr::log_print("Compartment:")
sample_variance <- get_compartment(FALSE) |>
    group_by(cell_type) |>
    summarize(Mean=mean(proportion),
              SD=sd(proportion),
              Variance=var(proportion))
logr::log_print(sample_variance |> as.data.frame())
logr::log_print("Cell types: Marker Genes")
sample_variance <- get_cellprop(FALSE) |>
    group_by(cell_type) |>
    summarize(Mean=mean(proportion),
              SD=sd(proportion),
              Variance=var(proportion))
logr::log_print(sample_variance |> as.data.frame())
logr::log_print("Cell types: All Genes")
sample_variance <- get_cellprop(TRUE) |>
    group_by(cell_type) |>
    summarize(Mean=mean(proportion),
              SD=sd(proportion),
              Variance=var(proportion))
logr::log_print(sample_variance |> as.data.frame())
logr::log_close()

                                        # Proportion plotting
plotNsave_prop(get_compartment, "compartment", 5, 4)
plotNsave_bar(get_compartment, "compartment", 4, 6)
plotNsave_prop(get_cellprop, "celltype.marker_genes", 14, 8, FALSE)
plotNsave_pie(FALSE)
plotNsave_prop(get_cellprop, "celltype.all_genes", 14, 8, TRUE)
plotNsave_pie(TRUE)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
