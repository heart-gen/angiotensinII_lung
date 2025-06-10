library(dplyr)
library(ggpubr)
library(Seurat)
library(scPower)
library(SingleCellExperiment)

#### FUNCTIONS
load_data <- function() {
                                        # Load data
    lung <- readRDS("../_h/hlca_core.rds")
    sce  <- as.SingleCellExperiment(lung)
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$cell_type   <- sce$ann_coarse_for_GWAS_and_modeling
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$individual  <- sce$donor_id
    colLabels(sce)           <- sce$cell_type
    return(sce)
}

filter_stroma_data <- function(sce) {
                                        # Filter for mesenchymal cells
    stroma <- sce[, colData(sce)$compartment == "Stroma"]
    colData(stroma)$cell_type <- droplevels(colData(stroma)$cell_type)
    return(as.data.frame(colData(stroma)))
}

get_cell_fraction <- function(ct, meta_df) {
                                        # Calculate the cell type frequency
    num <- dim(filter(meta_df, cell_type == ct))[1]
    tot <- dim(meta_df)[1]
    return(num / tot)
}

subset_data <- function(sce, ct) {
    sce        <- sce[, colData(sce)$cell_type == ct]
    counts.mat <- assays(sce)$counts |> as.matrix()
    annot <- colData(sce) |> as.data.frame() |>
        select(cell_type, individual, clusters, subclusters) |>
        tibble::rownames_to_column("cell")
    rownames(annot) <- annot$cell
                                        # Make sure annotation data is in
                                        # same order as count matrix
    all(colnames(counts.mat) == annot$cell)
    return(list("matrix"=counts.mat, "annotation"=annot))
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

round_up <- function(x, to = 1000) { ceiling(x / to) * to }

extract_parameters <- function(counts.matrix, label) {
                                        # Estimate negative binomial
    temp <- nbinom.estimation(counts.matrix, sizeFactorMethod = "poscounts")
    norm.mean.values <- temp[[1]]
    disp.param       <- temp[[3]]
                                        # Estimate gamma mixed distribution
                                        # Number of cells per cell type as
                                        # censoring point
    censoredPoint  <- 1 / ncol(counts.matrix)
    num.genes.kept <- round_up(sum(norm.mean.values$mean != 0), 1000)
    gamma.fit <- mixed.gamma.estimation(norm.mean.values$mean,
                                        num.genes.kept = num.genes.kept,
                                        censoredPoint = censoredPoint)
                                        # Visualize gamma
    g1      <- visualize.gamma.fits(norm.mean.values$mean, gamma.fit,
                                    nGenes=num.genes.kept)
    outfile <- paste0("gamma_visualize.", tolower(label), ".pdf")
    ggsave(outfile, plot=g1, width=9, height=6)
                                        # Parameterization of parameters
    umi.values <- data.frame("mean.umi"=meanUMI.calculation(counts.matrix))
    gamma.fit  <- convert.gamma.parameters(cbind(gamma.fit, umi.values))
    gamma.fit$ct <- "New_ct"
                                        # Estimate median dispersion
    disp.fun.general.new <- dispersion.function.estimation(disp.param)
    disp.fun.general.new$ct <- "New_ct"
    return(list("gamma"=gamma.fit, "disp"=disp.fun.general.new,
                "nGenes"=num.genes.kept))
}

#### MAIN
                                        # Load data
sce     <- load_data()
meta_df <- filter_stroma_data(sce)

                                        # Power estimation to detection
ct      <- "Pericytes"
frac    <- get_cell_fraction(ct, meta_df)
ns <- c(); cells <- c(); cell_per <- c()
for(n_samples in c(10, 15, 30, 60)) {
    for(min_cell in seq(25, 150, by=25)) {
        est_frac <- number.cells.detect.celltype(prob.cutoff = 0.95,
                                                 min.num.cells=min_cell,
                                                 cell.type.frac = frac,
                                                 nSamples=n_samples)
        ns       <- c(ns, n_samples)
        cells    <- c(cells, min_cell)
        cell_per <- c(cell_per, est_frac)
    }
}
power_df <- data.frame(N=ns, Min=cells, Cells=cell_per) |>
    mutate(N = as.factor(as.character(N)))

xlabel <- "Minimal number of cells from pericytes per individual (95% probability)"
gline  <- ggline(power_df, x="Min", y="Cells", linetype="N", shape="N",
                 size=0.75, ylab="Cells per Individual", xlab=xlabel,
                 ggtheme=theme_pubr(base_size=16)) +
    font("xy.title", size=18, face="bold")
pwrfile <- paste0("power_detection.", tolower(ct), ".pdf")
ggsave(pwrfile, plot=gline, width=12, height=5)

                                        # Dispersion priors
data <- subset_data(sce, ct)

                                        # Number of expressed genes
print_expressed_genes(data[["annotation"]], data[["matrix"]])

                                        # Dispersion and gamma parameters
model_priors <- extract_parameters(data[["matrix"]], ct)

                                        # Simulate effect size in DE
ranks      <- uniform.ranks.interval(start=5001, end=10000, numGenes=200)
foldChange <- effectSize.DE.simulation(mean=2, sd=0.5, numGenes=200)
simulated.de.genes <- data.frame(ranks=ranks, FoldChange=foldChange,
                                 name="Simulated")

detect_pwr <- c()
for(n_samples in c(10, 15, 30, 60)) {
    for(ncells in seq(500, 4000, by=500)) {
        power <- power.sameReadDepth.withDoublets(
            nSamples = n_samples, nCells = ncells, ct.freq = frac,
            type = "de", ref.study = simulated.de.genes,
            ref.study.name = "Simulated", samplesPerLane = 4,
            gamma.parameters = model_priors[["gamma"]], ct="New_ct",
            disp.fun.param = model_priors[["disp"]],
            mappingEfficiency = 0.8, min.UMI.counts = 3,
            perc.indiv.expr = 0.5, nGenes = model_priors[["nGenes"]],
            sign.threshold = 0.05, MTmethod = "Bonferroni")
        detect_pwr <- rbind(detect_pwr, power)
    }
}
detect_df <- data.frame(detect_pwr)
detect_df |> data.table::fwrite(paste0("power_analysis.", tolower(ct), ".tsv"),
                                sep="\t")

dotchart <- ggdotchart(detect_df, x="sampleSize", y="totalCells",
                       color="powerDetect", ylab="Cells per Individual",
                       size=8, xlab="Sample Size",
                       ggtheme=theme_pubr(base_size=16, legend="right")) +
    font("xy.title", size=18, face="bold") +
    gradient_color(viridisLite::viridis(10)) +
    rotate_x_text(angle=0, hjust=0.5) +
    labs(color="Detection Power")

detfile <- paste0("power_sampleSize.", tolower(ct), ".pdf")
ggsave(detfile, plot=dotchart, width=6, height=5)


#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
