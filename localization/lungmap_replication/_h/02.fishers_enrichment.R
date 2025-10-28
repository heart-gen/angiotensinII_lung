## Test for enrichment
suppressPackageStartupMessages({
    library(dplyr)
})

fishers_exact_compartment <- function(dt, xlab){
    dx        <- dt |> group_by(Patient, Compartment) |>
        summarize(mean_expr=mean(`Normalized Expression`,na.rm=TRUE))
    group_yes <- dim(filter(dx,Compartment == xlab,mean_expr > 0))[1]
    group_no  <- dim(filter(dx,Compartment == xlab,mean_expr ==0))[1]
    loc_yes   <- dim(filter(dx,Compartment != xlab,mean_expr > 0))[1]
    loc_no    <- dim(filter(dx,Compartment != xlab,mean_expr ==0))[1]
    table <- data.frame("present"=c(group_yes, loc_yes),
                        "not_present"=c(group_no, loc_no),
                        row.names=c("In Compartment", "Not in Compartment"))
    colnames(table) <- c("Present", "Not Present")
    return(fisher.test(table))
}

fishers_exact_annotation <- function(dt, xlab){
    dx        <- dt |> group_by(Patient, Cell_Annotation) |>
        summarize(mean_expr=mean(`Normalized Expression`,na.rm=TRUE))
    group_yes <- dim(filter(dx,Cell_Annotation == xlab,mean_expr > 0))[1]
    group_no  <- dim(filter(dx,Cell_Annotation == xlab,mean_expr ==0))[1]
    loc_yes   <- dim(filter(dx,Cell_Annotation != xlab,mean_expr > 0))[1]
    loc_no    <- dim(filter(dx,Cell_Annotation != xlab,mean_expr ==0))[1]
    table <- data.frame("present"=c(group_yes, loc_yes),
                        "not_present"=c(group_no, loc_no),
                        row.names=c("In Cell type", "Not in Cell type"))
    colnames(table) <- c("Present", "Not Present")
    return(fisher.test(table))
}

fishers_exact_subcluster <- function(dt, xlab){
    dx        <- dt |> group_by(Patient, Subcluster) |>
        summarize(mean_expr=mean(`Normalized Expression`,na.rm=TRUE))
    group_yes <- dim(filter(dx,Subcluster == xlab,mean_expr > 0))[1]
    group_no  <- dim(filter(dx,Subcluster == xlab,mean_expr ==0))[1]
    loc_yes   <- dim(filter(dx,Subcluster != xlab,mean_expr > 0))[1]
    loc_no    <- dim(filter(dx,Subcluster != xlab,mean_expr ==0))[1]
    table <- data.frame("present"=c(group_yes, loc_yes),
                        "not_present"=c(group_no, loc_no),
                        row.names=c("In Cell type", "Not in Cell type"))
    colnames(table) <- c("Present", "Not Present")
    return(fisher.test(table))
}

enrichment_loop <- function(df, fnc, locs, label){
    datalist  <- list(); locations <- c();
    gnames <- c(); pvalues <- c(); oddratio <- c();
    for(gene in c("AGTR1","AGTR2")){
        dt <- df |> filter(Gene_Name == gene)
        for(location in locs){
            gnames    <- c(gnames, gene)
            locations <- c(locations, location)
            pvalues   <- c(pvalues,fnc(dt,location)$p.value)
            oddratio  <- c(oddratio,fnc(dt,location)$estimate)
        }
    }
    fdr <- p.adjust(pvalues, method="fdr")
    return(data.frame("Annotation"=label,
                      "Annotation_Name"=locations,
                      "Gene"=gnames, "OR"=oddratio,
                      "P"=pvalues, "FDR"=fdr))
}

enrichment_analysis <- function(){
    df <- data.table::fread("../../_m/normalized_expression.txt.gz")
                                        # Compartment
    locs <- df$Compartment |> unique()
    dat1 <- enrichment_loop(df, fishers_exact_compartment,
                            locs, "Compartment")
                                        # Cell type
    locs <- df$Cell_Annotation |> unique()
    dat2 <- enrichment_loop(df, fishers_exact_annotation,
                            locs, "Cell type")
                                        # Subclusters
    locs <- df$Subcluster |> unique()
    dat3 <- enrichment_loop(df, fishers_exact_subcluster,
                            locs, "Subcluster")
    return(bind_rows(dat1, dat2, dat3))
}

#### MAIN
enrichment_analysis() |>
    data.table::fwrite("cell_annotation_enrichment_analysis.tsv",
                       sep='\t')

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
