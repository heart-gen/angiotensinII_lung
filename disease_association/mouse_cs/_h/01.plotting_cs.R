## This script plots the qPCR data from mouse CS model

library(ggpubr)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

#### MAIN
                                        # Load data
df <- data.table::fread("../_h/Mouse-receptor-data-with-smoke-exposure.csv")
                                        # Format for plotting
df <- tidyr::pivot_longer(df, !V1, names_to="samples",
                          values_to="relative") |>
    as.data.frame() |> tidyr::drop_na() |> 
    dplyr::mutate(receptor=ifelse(V1=="AT1","AGTR1","AGTR2"),
                  samples=ifelse(samples=="CS","Cigarette Smoke","Room Air"))
                                        # Generate plot
bxp <- ggboxplot(df, x="samples", y="relative", fill="samples",
                 facet.by="receptor", palette="npg",
                 outlier.shape=NA, add="jitter",
                 ylab="Relative fold change", xlab="", ylim=c(0,2),
                 add.params=list(alpha=0.8), legend="none",
                 panel.labs.font=list(face='bold', size=16), 
                 ggtheme=theme_pubr(base_size=15, border=TRUE)) +
    font("ylab", size=18, face="bold") + rotate_x_text(45) +
    stat_compare_means(method="t.test", label.y=1.9)
                                        # Save plot
save_ggplots("mouse_receptor.smoke_exposure", bxp, 5, 6)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
