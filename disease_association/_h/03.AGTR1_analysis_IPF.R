## This script examines AGTR1-positive cell
## distribution across disease state.

library(dplyr)

load_data <- function(){
    return(data.table::fread("normalized_expression.tsv.gz") |>
           filter(Gene_Name == "AGTR1") |>
           group_by(Patient, Disease) |>
           reframe(Mean_Expr = mean(`Normalized Expression`)) |>
           mutate(Present = ifelse(Mean_Expr > 0.001, 1, 0)) |>
           as.data.frame() |> mutate_if(is.character, as.factor))
}

#### Main
df <- load_data()
group_by(df, Disease) |>
    summarize(Positive= sum(Present),
              N = n(),
              Prop = Positive / N) |>
    as.data.frame()
res_aov <- aov(Mean_Expr ~ Disease, data=df)
summary(res_aov)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
