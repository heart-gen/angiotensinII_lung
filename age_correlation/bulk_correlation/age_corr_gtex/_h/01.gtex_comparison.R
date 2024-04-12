## This script compares TPM expression across GTEx v8 tissues

library(here)
library(dplyr)
library(ggpubr)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

get_tpm <- function(){
    fn = here("inputs/gtex/_m",
              "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
    return(data.table::fread(fn))
}
memTPM <- memoise::memoise(get_tpm)

get_pheno <- function(){
    fn = here("inputs/gtex/_m/gtex_v8_sample_data.tsv")
    return(data.table::fread(fn))
}
memPHENO <- memoise::memoise(get_pheno)

get_angiotensin <- function(){
    return(filter(memTPM(), Description %in% c("AGTR1", "AGTR2")) |>
           tibble::column_to_rownames("Description") |>
           select(-"Name") |> t() |> as.data.frame())
}
memAGTR <- memoise::memoise(get_angiotensin)

merge_data <- function(){
    return(memAGTR() |> tibble::rownames_to_column("SAMPID") |>
           tidyr::pivot_longer(!SAMPID, names_to="Symbol", values_to="TPM") |>
           mutate(log_TPM=log10(TPM+1)) |>
           inner_join(memPHENO(), by="SAMPID"))
}
memDF <- memoise::memoise(merge_data)

#### Main
                                        # Select Lung tissues
dt <- memDF() |> filter(SMTS == "Lung")
                                        # Correlation with age
lf <- logr::log_open("summary.log", logdir=FALSE, autolog=FALSE)
logr::log_print("Summary of correlation between log2 TPM expression and age:")
est_fit <- dt |> group_by(Symbol) |>
    do(fit=broom::tidy(lm(log_TPM ~ AGE, data=.))) |> tidyr::unnest(fit) |>
    filter(term != "(Intercept)") |>
    mutate(p.bonf = p.adjust(p.value, "bonf"), p.bonf.sig = p.bonf < 0.05,
           p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                            labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                            include.lowest = TRUE),
           p.fdr = p.adjust(p.value, "fdr"),
           log.p.bonf = -log10(p.bonf))
logr::log_print(est_fit |> count(p.bonf.cat))
logr::log_print(est_fit |> as.data.frame())
logr::log_close()
                                        # Plot correlation
sca <- ggscatter(dt, x="AGE", y="log_TPM", facet.by="Symbol", add="reg.line",
                 add.params=list(color="blue", fill="lightgray"),
                 conf.int=TRUE, panel.labs.font=list(face="bold"),
                 cor.coef=TRUE, ylab="log10(TPM + 1)",
                 cor.coeff.args=list(method="pearson",label.sep='\n'),
                 ggtheme=theme_pubr(base_size=15, border=TRUE)) +
    font("xy.title", face="bold", size=18)
save_ggplots("gtex_v8_lung_age_corr", sca, 8, 4)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
