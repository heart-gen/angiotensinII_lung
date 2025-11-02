## Plotting enrichment
suppressPackageStartupMessages({
    library(ggplot2)
    library(tidyverse)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png")){
        ggsave(filename=paste0(fn,ext), plot=p,
               width=w, height=h)
    }
}

load_enrichment <- function(){
    return(data.table::fread("cell_annotation_enrichment_analysis.tsv"))
}
memENRICH <- memoise::memoise(load_enrichment)

gen_data <- function(){
    err = 0.0000001; err2 = 1e-100
    dt <- memENRICH() |> mutate_if(is.character, as.factor) |>
        mutate(`-log10(FDR)`= ifelse(FDR != 0, -log10(FDR), -log10(err)),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE))
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(label, w, h){
    dt <- memDF() |> filter(Annotation == label)
    y0 <- min(dt$`log2(OR)`)-0.1
    y1 <- max(dt$`log2(OR)`)+0.1
    tile_plot <- ggplot(dt, aes(x = Gene, y = Annotation_Name,
                                fill = `log2(OR)`,
                                label = ifelse(p.fdr.sig,
                                               format(round(`-log10(FDR)`,1),
                                                      nsmall=1), ""))) +
        ylab('') + xlab("") + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=FALSE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              legend.position="right",
              axis.title=element_text(face="bold"),
              axis.text.y=element_text(face="bold"),
              strip.text = element_text(face="bold"))
    outfile <- paste0("tileplot_enrichment_",
                      gsub(" ", "_",tolower(label)))
    save_plot(tile_plot, outfile, w, h)
}

## Run script
plot_tile("Compartment", 5.4, 4)
plot_tile("Cell type", 8, 13)
plot_tile("Subcluster", 8, 17)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
