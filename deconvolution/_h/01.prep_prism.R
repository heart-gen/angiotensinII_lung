## Prepare data for deconvolution using HLCA version 2.

suppressPackageStartupMessages({
    library(dplyr)
    library(BayesPrism)
    library(SingleCellExperiment)
})

source("../_h/get_mean_ratio2.R")

get_counts <- function(){
    fn = here::here("inputs/gtex/_m/",
                    "genes_gtex_v8_counts.txt.gz")
    return(data.table::fread(fn))
}
memCOUNTS <- memoise::memoise(get_counts)

get_pheno <- function(){
    fn = here::here("inputs/gtex/_m/gtex_v8_sample_data.tsv")
    return(data.table::fread(fn) |> filter(SMTS == "Lung"))
}
memPHENO <- memoise::memoise(get_pheno)

select_lung_data <- function(){
                                        # Clean data
    counts <- tibble::column_to_rownames(memCOUNTS(), "Name") |>
        select(any_of(memPHENO()$SAMPID))
    genes  <- memCOUNTS() |> select(Name, Description)
    pheno  <- memPHENO() |> filter(SAMPID %in% colnames(counts))
                                        # Filter low expression
    x <- edgeR::DGEList(counts=counts, genes=genes, samples=pheno)
    keep.x <- edgeR::filterByExpr(x)
    x <- x[keep.x, , keep.lib.sizes=FALSE]
                                        # Update rownames
    bulk_counts <- x$counts |> as.data.frame() |>
        tibble::rownames_to_column("Name") |>
        inner_join(genes, by="Name") |>
        distinct(Description, .keep_all=TRUE) |>
        tibble::column_to_rownames("Description") |>
        select(-Name)
    return(bulk_counts)
}
memDF <- memoise::memoise(select_lung_data)

#### Main
                                        # Prepare GTEx lung bulk data
bk.dat <- memDF() |> t()

                                        # Load single cell data
load("../_m/scRNA_HLCA_version2.RData")
rownames(sce) <- rowData(sce)[, "feature_name"]

                                        # Generate meta data
meta <- colData(sce) |> as.data.frame() |>
    select(patient, cell_type, compartment, subclusters)

                                        # Prepare single cell reference
sc.dat        <- assays(sce)$counts |> t()
plot.cor.phi(input=sc.dat, pdf.prefix="lung.cor",
             input.labels=meta$cell_type,
             title="Cell type correlation",
             cexRow=0.2, cexCol=0.2, 
             margins=c(2,2))

                                        # Select marker genes
marker_stats <- get_mean_ratio2(sce, cellType_col="cell_type")
rm(sce); gc()

                                        # Quality control plots
plot.scRNA.outlier(
    input=sc.dat, 
    pdf.prefix="lung.sc_stat",
    cell.type.labels=meta$cell_type,
    species="hs", return.raw=FALSE
)
plot.bulk.outlier(
    bulk.input=bk.dat, sc.input=sc.dat,
    cell.type.labels=meta$cell_type,
    species="hs", return.raw=FALSE,
    pdf.prefix="lung.bk_stats"
)

                                        # Cleanup genes
sc.filtered <- cleanup.genes(input=sc.dat, input.type="count.matrix",
                             species="hs", exp.cells=5,
                             gene.group=c("Rb", "Mrp", "other_Rb",
                                          "chrM", "chrX", "chrY", "MALAT1"))

plot.bulk.vs.sc(sc.input=sc.filtered,
                bulk.input=bk.dat,
                pdf.prefix="lung.bk_vs_sc")
rm(sc.dat); gc()

                                        # Subset genes to improve 
                                        # computational efficiency
marker_genes <- marker_stats |>
    filter(rank_ratio <= 100, gene %in% colnames(bk.dat),
           gene %in% colnames(sc.filtered))
print(table(marker_genes$cellType.target))
sc.filtered.sig <- sc.filtered[, pull(marker_genes, gene)]

                                        # Run BayesPrism
lungPrism <- new.prism(reference=sc.filtered.sig,
                       mixture=bk.dat, input.type="count.matrix",
                       cell.type.labels=meta$cell_type,
                       cell.state.labels=meta$cell_type,
                       key=NULL, outlier.cut=0.01,
                       outlier.fraction=0.1)

lungPrism0 <- new.prism(reference=sc.filtered.sig,
                        mixture=bk.dat, input.type="count.matrix",
                        cell.type.labels=meta$compartment,
                        cell.state.labels=meta$cell_type,
                        key=NULL, outlier.cut=0.01,
                        outlier.fraction=0.1)

lungPrismX <- new.prism(reference=sc.filtered,
                        mixture=bk.dat, input.type="count.matrix",
                        cell.type.labels=meta$cell_type,
                        cell.state.labels=meta$cell_type,
                        key=NULL, outlier.cut=0.01,
                        outlier.fraction=0.1)

                                        # Save results
save(lungPrism, lungPrism0, lungPrismX,
     file="lung_prop_celltypes.RDat")

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
