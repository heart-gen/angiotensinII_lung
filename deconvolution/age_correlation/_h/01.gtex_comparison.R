## This script compares TPM expression across GTEx v8 tissues
##
## Normalize expression with count data
#############################################################

library(here)
library(dplyr)
library(ggpubr)
library(BayesPrism)

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
        janitor::clean_names()
    return(theta)
}

get_pheno <- function(FULL){
    fn = here("inputs/gtex/_m/gtex_v8_sample_data.tsv")
    return(data.table::fread(fn) |> filter(SMTS == "Lung") |>
           janitor::clean_names() |>
           inner_join(get_cellprop(FULL), by=c("sampid"="sample_id")))
}
memPHENO <- memoise::memoise(get_pheno)

get_counts <- function(){
    fn = here("inputs/gtex/_m/genes_gtex_v8_counts.txt.gz")
    return(data.table::fread(fn))
}
memCOUNTS <- memoise::memoise(get_counts)

select_lung_data <- function(FULL){
                                        # Clean data
    counts <- tibble::column_to_rownames(memCOUNTS(), "Name") |>
        select(any_of(memPHENO(FULL)$sampid))
    genes  <- memCOUNTS() |> select(Name, Description)
    pheno  <- memPHENO(FULL) |> filter(sampid %in% colnames(counts))
                                        # Filter low expression
    x      <- edgeR::DGEList(counts=counts, genes=genes, samples=pheno)
    mod    <- model.matrix(~age, data=x$samples)
    keep.x <- edgeR::filterByExpr(x, design=mod)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x      <- x[keep.x, , keep.lib.sizes=FALSE]
                                        # Normalize counts
    x      <- edgeR::calcNormFactors(x, method="TMM")
                                        # Build model
    ct_model  <- paste(colnames(x$samples[, 22:60]), collapse=" + ")
    new_model <- paste("~age + sex + race", ct_model, sep=" + ")
    mod       <- model.matrix(formula(new_model), data=x$samples)
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept",
                          colnames(mod))
                                        # Voom normalize
    v          <- limma::voom(x, mod)
                                        ## # Residualize expression
    ## null_model <- v$design |> as.data.frame() |>
    ##     select(-c("age")) |> as.matrix()
    ## fit_res    <- limma::lmFit(v, design=null_model)
    ## res        <- v$E - ( fit_res$coefficients %*% t(null_model) )
    ## res_sd     <- apply(res, 1, sd); res_mean = apply(res, 1, mean)
    ## res_norm   <- (res - res_mean) / res_sd
                                        # Annotate genes
    ## res_norm   <- as.data.frame(res_norm) |>
    ##     tibble::rownames_to_column("Name")
    expr  <- as.data.frame(v$E) |>
        tibble::rownames_to_column("Name")
    return(inner_join(genes, expr, by="Name"))
}
memRES <- memoise::memoise(select_lung_data)

get_angiotensin <- function(FULL){
    return(filter(memRES(FULL), Description %in% c("AGTR1", "AGTR2")) |>
           tibble::column_to_rownames("Description") |>
           select(-"Name") |> t() |> as.data.frame())
}
memAGTR <- memoise::memoise(get_angiotensin)

merge_data <- function(FULL){
    return(tibble::rownames_to_column(memAGTR(FULL), "SAMPID") |>
           tidyr::pivot_longer(!SAMPID,
                               names_to="Symbol",
                               values_to="Expression") |>
           janitor::clean_names() |>
           inner_join(memPHENO(FULL), by="sampid"))
}
memDF <- memoise::memoise(merge_data)

plotting_age_corr <- function(dt, yvar, ylabel, mlabel){
    fn  <- paste("gtex_v8.lung.age_corr", yvar, mlabel, sep=".")
    sca <- ggscatter(dt, x="age", y=yvar, facet.by="symbol", add="reg.line",
                     add.params=list(color="blue", fill="lightgray"),
                     conf.int=TRUE, panel.labs.font=list(face="bold"),
                     cor.coef=TRUE, xlab="Age", ylab=ylabel,
                     cor.coeff.args=list(method="pearson",label.sep='\n'),
                     ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=18)
    save_ggplots(fn, sca, 8, 4)   
}

plotting_expression_corr <- function(dt, xvar, xlabel, mlabel){
    fn  <- paste("gtex_v8.lung", xvar, "corr", mlabel, sep=".")
    sca <- ggscatter(dt, x=xvar, y="expression",
                     facet.by="symbol", add="reg.line",
                     add.params=list(color="blue", fill="lightgray"),
                     conf.int=TRUE, panel.labs.font=list(face="bold"),
                     cor.coef=TRUE, xlab=xlabel,
                     ylab="Normalized Expression",
                     cor.coeff.args=list(method="pearson",label.sep='\n'),
                     ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=18)
    save_ggplots(fn, sca, 8, 4)
}

#### Main
dt1 <- memDF(FALSE)
dt2 <- memDF(TRUE)
                                        # Correlation with age
model <- "expression ~ age"
lf    <- logr::log_open("summary.log", logdir=FALSE, autolog=FALSE)
logr::log_print("Summary of correlation between residualized expression and age:")
logr::log_print("Marker genes")
est_fit <- dt1 |> group_by(symbol) |>
    do(fit=broom::tidy(lm(model, data=.))) |>
    tidyr::unnest(fit) |>
    filter(term != "(Intercept)") |>
    mutate(p.bonf = p.adjust(p.value, "bonf"),
           p.bonf.sig = p.bonf < 0.05,
           p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                            labels = c("<= 0.005","<= 0.01",
                                       "<= 0.05","> 0.05"),
                            include.lowest = TRUE),
           p.fdr = p.adjust(p.value, "fdr"),
           log.p.bonf = -log10(p.bonf))
logr::log_print(est_fit |> count(p.bonf.cat))
logr::log_print(est_fit |> as.data.frame())

logr::log_print("All genes")
est_fit <- dt2 |> group_by(symbol) |>
    do(fit=broom::tidy(lm(model, data=.))) |>
    tidyr::unnest(fit) |>
    filter(term != "(Intercept)") |>
    mutate(p.bonf = p.adjust(p.value, "bonf"),
           p.bonf.sig = p.bonf < 0.05,
           p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                            labels = c("<= 0.005","<= 0.01",
                                       "<= 0.05","> 0.05"),
                            include.lowest = TRUE),
           p.fdr = p.adjust(p.value, "fdr"),
           log.p.bonf = -log10(p.bonf))
logr::log_print(est_fit |> count(p.bonf.cat))
logr::log_print(est_fit |> as.data.frame())
logr::log_close()
                                        # Plot correlation
vars <- c("expression", "Pericytes", "AT2",
          "Peribronchial_Fibroblasts",
          "Smooth_Muscle", "Alveolar_Fibroblasts",
          "Adventitial_Fibroblasts",
          "Alveolar_Macrophages")
for(var in vars){
    if(var == "expression"){
        ylab <- "Normalized Expression"
        ## Age Correlation
        plotting_age_corr(dt1, tolower(var), ylab, "sig")
        plotting_age_corr(dt2, tolower(var), ylab, "all")
    } else {
        ylab <- paste("Proportion of", gsub("_", " ", var))
        ## Age Correlation
        plotting_age_corr(dt1, tolower(var), ylab, "sig")
        plotting_age_corr(dt2, tolower(var), ylab, "all")
        ## Proportion Correlation
        plotting_expression_corr(dt1, tolower(var), ylab, "sig")
        plotting_expression_corr(dt2, tolower(var), ylab, "all")
    }    
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
