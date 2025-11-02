## Plotting PHATE results
suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

plot_box_stats <- function(){
    dt <- data.table::fread("at2_phate.normalized_expression.tsv.gz",
                            sep="\t")
    mycomps <- list(c("Control", "COPD"), c("Control", "IPF"),
                    c("COPD", "IPF"))
    outfile <- "at2_phate.normalized_expression.boxplot_disease"
    bxp <- dt |> filter(AGTR2 > 0) |> 
        ggboxplot(x="disease", y="AGTR2", fill="disease",
                  facet.by="PHATE", add="jitter", legend="bottom",
                  outlier.shape=NA, xlab="", palette="npg",
                  ylab="Normalized Expression (AGTR2)",
                  panel.labs.font=list(face='bold'),
                  add.params=list(alpha=0.8), ncol=1,
                  ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        rotate_x_text(45) +
        stat_compare_means(comparisons=mycomps, method="t.test")
    save_ggplots(tolower(outfile), bxp, 3, 6)
}

### Main script section
plot_box_stats()

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
