## Generation of new priors from Lung Atlas

library(dplyr)
library(ggpubr)
library(Seurat)
library(scPower)

round_up <- function(, to = 1000) {
                                        # Round data
    return(ceiling(x / to) * to)
}

get_cell_fraction <- function(ct, meta_df) {
                                        # Calculate the cell type frequency
    num <- dim(filter(meta_df, cell_type == ct))[1]
    tot <- dim(meta_df)[1]
    return(num / tot)
}

read_data <- function(fn) {
                                        # Load data
    col_names <- c("gene_id", "snp_id", "tss", "p_value", "beta")
    dt <- data.table::fread(fn, header=FALSE, col.names=col_names)
                                        # Calculate Rsq
    n  <- 24 + 40 + 83 # Three datasets Bryois, et al.
    df <- n - 2; var_x <- 1; var_epsilon <- 1
    return(mutate(dt, var_hat_y = beta^2 * var_x,
                  var_y = var_hat_y + var_epsilon,
                  Rsq = var_hat_y / var_y,
                  name = "Bryois"))
}

print_expressed_genes <- function(annot, counts.matrix) {
                                        # Filter the annotation data
    annot       <- annot[colnames(counts.matrix),]
                                        # Reformat to pseudobulk matrix
    pseudo.bulk <- create.pseudobulk(counts.matrix, annot, colName="cell_type")
                                        # Calculate expressed genes
    expressed.genes <- calculate.gene.counts(pseudo.bulk, min.counts=3,
                                             perc.indiv=0.5)
    df <- data.frame(num.cells=ncol(counts.matrix),
                     expressed.genes=nrow(expressed.genes))
    print(df)
}

subset_data <- function(seur, meta_df, ct) {
    seur  <- subset(x=seur, subset = cell_type == ct)
    coutns.matrix <- seur[["RNA"]]$counts |> as.matrix()
    annot <- select(meta_df, cell_type, donor_id) |>
        filter(cell_type == ct) |>
        tibble::rownames_to_column("cell") |>
        rename("individual" = "donor_id")
    rownames(annot) <- annot$cell
                                        # Make sure annotation data is in
                                        # same order as count matrix
    all(colnames(counts.matrix) == annot$cell)
    return(list("matrix"=as.matrix(seur[["RNA"]]$counts), "annotation"=annot))
}

extract_parameters <- function(counts.matrix, label) {
    ##                                     # Remove non-expressing genes
    ## counts.matrix <- counts.matrix[rowSums(counts.matrix != 0) > 0, ]
                                        # Estimate negative binomial
    temp <- nbinom.estimation(counts.matrix)
    norm.mean.values <- temp[[1]]
    disp.param <- temp[[3]]
                                        # Estimate gamma mixed distribution
                                        # Number of cells per cell type as
                                        # censoring point
    censoredPoint  <- 1 / ncol(counts.matrix)
    num.genes.kept <- round_up(sum(norm.mean.values$mean != 0), 10000)
    gamma.fit <- mixed.gamma.estimation(norm.mean.values$mean,
                                        num.genes.kept = num.genes.kept,
                                        censoredPoint = censoredPoint)
                                        # Visualize gamma
    g1 <- visualize.gamma.fits(norm.mean.values$mean, gamma.fit,
                               nGenes=num.genes.kept)
    ct_name <- gsub(" ", "_", tolower(label))
    outfile <- paste0("gamma_visualize.", ct_name, ".pdf")
    ggsave(outfile, plot=g1, width=9, height=6)
                                        # Parameterization of parameters
    umi.values <- data.frame("mean.umi"=meanUMI.calculation(counts.matrix))
    gamma.fit  <- convert.gamma.parameters(cbind(gamma.fit, umi.values))
    gamma.fit$ct <- "New_ct"
    return(list("gamma"=gamma.fit, "disp"=disp.fun.general.new,
                "nGenes"=num.genes.kept))
}

#### Main
                                        # Load data
lung <- readRDS(here::here("inputs/hlca/_m/hlca_core.rds"))
mat  <- lung[["RNA"]]$counts
mat  <- mat[rowSums(mat)>0, ] # remove all 0 genes
mat  <- t(t(mat) / colSums(mat)) # Normalize by count per cell
