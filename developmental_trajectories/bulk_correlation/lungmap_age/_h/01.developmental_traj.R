## This script examines angiotensin II in bulk lung tissues
## across different developmental time points.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

source("../_h/age_plotter.R")
save_ggplot <- function(fn, pp, w, h){
    for(ext in c(".svg", ".pdf")){
        ggsave(paste0(fn, ext), plot=pp, width=w, height=h)
    }
}

get_phenotypes <- function(){
    fn = here::here("inputs/lungmap/_m/rnaseq",
                    "LMEX0000003691_sample_metadata.xlsx")
    df <- readxl::read_excel(fn) |>
        mutate(Age=sapply(gsub(" Days","/365",
                               gsub(" Years","/1",
                                    gsub(" Months", "/12", `Age Of Donor`))),
                          function(x) eval(parse(text=x)))) |>
        select(c("Donor Id", "Sample Inventory ID", "Age", "Donor Sex",
                 "Donor Race", "Age Cohort", "Derivative Type")) |>
        janitor::clean_names()
    return(df)
}
memPHENO <- memoise::memoise(get_phenotypes)

get_norm_angiotensin <- function(){
    fn <- here::here("inputs/lungmap/_m/rnaseq",
                     "LMEX0000003691_normalized_counts.txt.gz")
    return(data.table::fread(fn) |>
           filter(V1 %in% c("AGTR1", "AGTR2")) |>
           tibble::column_to_rownames("V1") |>
           t() |> as.data.frame() |>
           tibble::rownames_to_column("sample_inventory_id"))
}
memCOUNTS <- memoise::memoise(get_norm_angiotensin)

merge_data <- function(){
    return(memPHENO() |>
           inner_join(memCOUNTS(), by="sample_inventory_id"))
}
memDF <- memoise::memoise(merge_data)

plot_boxplot <- function(){
    dx <- memDF() |> tibble::column_to_rownames("sample_inventory_id") |>
        select(age, age_cohort, derivative_type, AGTR1, AGTR2) |>
        tidyr::pivot_longer(-c(age, age_cohort, derivative_type),
                            names_to="Gene_Name",
                            values_to="Normalized_Expression") |>
        mutate_if(is.character, as.factor) |>
        mutate(Expression=log10(Normalized_Expression+0.5))
    ## By derivative type
    bxp <- ggboxplot(dx, x="age_cohort", y="Expression", color="Gene_Name",
                     add="jitter", facet.by="derivative_type", xlab="",
                     palette="npg", ylab="Normalized Expression",
                     ggtheme=theme_pubr(base_size=15)) + rotate_x_text(90)
    save_ggplot("boxplot_angiotensinII_derivativeType", bxp, 8, 8)
    ## Grouped by age cohort
    bxp <- ggboxplot(dx, x="age_cohort", y="Expression", fill="age_cohort",
                     add="jitter", facet.by="Gene_Name", xlab="", legend="",
                     palette="npg", ylab="Normalized Expression",
                     ggtheme=theme_pubr(base_size=15)) + rotate_x_text(90)
    save_ggplot("boxplot_angiotensinII", bxp, 7, 6)
}

plot_line_graph <- function(){
    dx <- memDF() |> tibble::column_to_rownames("sample_inventory_id") |>
        select(age, age_cohort, derivative_type, AGTR1, AGTR2) |>
        tidyr::pivot_longer(-c(age, age_cohort, derivative_type),
                            names_to="Gene_Name",
                            values_to="Normalized_Expression") |>
        mutate_if(is.character, as.factor) |>
        mutate(Expression=log10(Normalized_Expression+0.5))
    ## By derivative type
    line <- ggline(dx, x="age_cohort", y="Expression", linetype="Gene_Name",
                   shape="Gene_Name", add=c("mean_se"), xlab="",
                   panel.labs.font=list(face='bold'),
                   facet.by="derivative_type", ylab="Normalized Expression",
                   ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        rotate_x_text(90)
    save_ggplot("linegraph_angiotensinII_derivativeType", line, 8, 8)
    ## Grouped by age cohort
    line <- ggline(dx, x="age_cohort", y="Expression", linetype="Gene_Name",
                   shape="Gene_Name", add=c("mean_se"),
                   xlab="", ylab="Normalized Expression",
                   ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        rotate_x_text(90)
    save_ggplot("linegraph_angiotensinII", line, 3, 6)
}

plot_scatter_breaks <- function(){
    dx = memDF() |> tibble::column_to_rownames("sample_inventory_id") |>
        mutate(AGTR1=log10(AGTR1+0.5), AGTR2=log10(AGTR2+0.5))
    for(gene_name in c("AGTR1", "AGTR2")){
        R.devices::devEval(
                       c("pdf", "png"), name="agePlotter_normalized_lung", path=".",
                       scale=1.25, tags=gene_name, sep="_", aspectRatio=0.6, {
                           agePlotter_v3(dx[, gene_name], dx[, "age"],
                                         mainText=gene_name,
                                         ylab="Normalized Expression",
                                         fetal_present=FALSE,
                                         ageBreaks=c(0, 1, 2, 50),
                                         mod=model.matrix(~age+donor_sex+donor_race, dx))
                       }
                   )
    }
}

main <- function(){
    lf <- logr::log_open("summary.log", logdir=FALSE, autolog=FALSE)
    logr::log_print("Summary of correlation between normalized expression and age:")
    for(gene_name in c("AGTR1", "AGTR2")){
        model <- paste0(gene_name, "~age+derivative_type+donor_sex+donor_race")
        fit0  <- lm(model, data=memDF())
        os    <- segmented::selgmented(fit0, seg.Z=~age,
                                       return.fit=FALSE,
                                       type='bic')
        if("npsi" %in% names(os)){
            npsi  <- os$npsi
        } else {
            npsi  <- os$selection.psi$npsi
        }
        if(npsi == 0){
            fit <- fit0
            dir <- coefficients(fit)[["age"]] |> sign()
            names(dir) <- "slope1"
        } else {
            fit <- segmented::segmented(fit0, seg.Z=~age, npsi=npsi)
            dir <- segmented::slope(fit)$age[, "Est."] |> sign()
        }
        logr::log_print(summary(fit))
    }
    logr::log_close()
    plot_scatter_breaks()
    plot_boxplot()
    plot_line_graph()
}

main()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
