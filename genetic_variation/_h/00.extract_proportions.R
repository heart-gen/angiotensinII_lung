## This script extracts the cell type proportions from R variable
library(dplyr)
library(BayesPrism)

get_cellprop <- function(){
    fn <- here::here("deconvolution/_m",
                     "prism_results.02.lung_GTEx.RDat")
    load(fn)
    return(theta <- get.fraction(bp=bp.res,
                                 which.theta="final",
                                 state.or.type="type") |>
               as.data.frame() |>
               tibble::rownames_to_column("Sample_ID") |>
               janitor::clean_names())
}

get_pheno <- function(){
    fn <- here::here("inputs/gtex/_m/gtex_v8_sample_data.tsv")
    return(data.table::fread(fn) |> filter(SMTS == "Lung") |>
           janitor::clean_names() |>
           inner_join(get_cellprop(),
                      by=c("sampid"="sample_id")))
}
memPHENO <- memoise::memoise(get_pheno)

#### Main
memPHENO() |>
    data.table::fwrite("gtex_phenotypes_lung.tsv", sep='\t')

get_cellprop() |>
    data.table::fwrite("gtex_cell_proportions.tsv", sep="\t")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
