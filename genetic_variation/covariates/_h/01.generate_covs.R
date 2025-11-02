## Prepare covariates for FastQTL / tensorQTL analysis
library(here)
library(dplyr)
                                        # Load data
sample_df <- data.table::fread("../../_m/sample_id_to_subject_id.tsv")
covs_file <- here("inputs/gtex/covariates/_m",
                  "GTEx_Analysis_v8_eQTL_covariates",
                  "Lung.v8.covariates.txt")

                                        # Extract covariates
covs <- data.table::fread(covs_file) |>
    select(ID, any_of(sample_df$subjid)) |>
    tibble::column_to_rownames("ID") |> t() |>
    as.data.frame() |> tibble::rownames_to_column("subjid") |>
    inner_join(sample_df, by=c("subjid")) |>
    select(-"sampid") |> rename("ID"="subjid") |>
    tibble::column_to_rownames("ID") |> t() |>
    as.data.frame() |> tibble::rownames_to_column("ID")

                                        # Save file
data.table::fwrite(covs, "genes.combined_covariates.txt", sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
