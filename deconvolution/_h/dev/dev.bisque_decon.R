## This script runs cellular deconvolution using BisqueRNA

suppressPackageStartupMessages({
    library(dplyr)
    library(BisqueRNA)
    library(SingleCellExperiment)
})

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

#### Prepare input ####

#### Bulk data     ####
                                        # Prepare GTEx lung bulk data
counts <- tibble::column_to_rownames(memCOUNTS(), "Name") |>
    select(any_of(memPHENO()$SAMPID))
genes  <- memCOUNTS() |> select(Name, Description)
pheno  <- memPHENO() |> filter(SAMPID %in% colnames(counts))
                                        # Filter low expression
x <- edgeR::DGEList(counts=counts, genes=genes, samples=pheno)
keep.x <- edgeR::filterByExpr(x)
x <- x[keep.x, , keep.lib.sizes=FALSE]
                                        # Update rownames
X <- x$counts |> as.data.frame() |>
    tibble::rownames_to_column("Name") |>
    inner_join(genes, by="Name") |>
    distinct(Description, .keep_all=TRUE) |>
    tibble::column_to_rownames("Description") |>
    select(-Name) |> as.matrix()
bulk.eset <- Biobase::ExpressionSet(assayData=X)
                                        # Clean variables
rm(x, pheno, genes, counts, X)
gc()

#### Single cell data ####
                                        # Prepare SCE
load("scRNA_HLCA_version2.RData", verbose=TRUE)
rownames(sce) <- rowData(sce)[, "feature_name"]
                                        # Generate reference
ref.data      <- t(counts(sce))
                                        # Generate meta data
meta    <- colData(sce) |> as.data.frame() |>
    select(patient, cell_type, compartment, subclusters) |>
    Biobase::AnnotatedDataFrame()
rm(sce); gc()
                                        # Cleanup genes
ref.filtered  <- BayesPrism::cleanup.genes(
    input=ref.data, input.type="count.matrix",
    species="hs", exp.cells=5,
    gene.group=c("Rb", "Mrp", "other_Rb", "chrM", "chrX", 
                 "chrY", "MALAT1")) |>
    t() |> as.matrix()
rm(ref.data); gc()
                                        # Make SCE variable
sc.eset <- Biobase::ExpressionSet(assayData=ref.filtered,
                                  phenoData=meta)
                                        # Clean variables
rm(ref.filtered, meta)
gc()

#### Run Bisque ####
run_bisque <- function(celltype){
    est_prop <- ReferenceBasedDecomposition(bulk.eset = bulk.eset,
                                            sc.eset = sc.eset,
                                            cell.types = celltype,
                                            subject.names = "patient",
                                            use.overlap = FALSE)
    est_prop$bulk.props <- t(est_prop$bulk.props)
    print(round(colMeans(est_prop$bulk.props),3))
                                        # Add long data
    est_prop$Est.prop.long <- est_prop$bulk.props |> as.data.frame() |>
        tibble::rownames_to_column("SAMPID") |>
        tidyr::pivot_longer(-SAMPID, names_to="Cell_Type",
                            values_to="Proprotion")
                                        # Isometric log ratio transformation
    est_prop$ilr <- compositions::ilr(est_prop$bulk.props)
    colnames(est_prop$ilr) <- paste0("ilr_",1:ncol(est_prop$ilr))
                                        # Save R variable
    outfile <- paste0(celltype,"_est_prop_Bisque.Rdata")
    save(est_prop, file = outfile)
}

run_bisque("compartment")
run_bisque("cell_type")
run_bisque("subclusters")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
