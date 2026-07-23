## Data-driven pericyte expression programs via CoGAPS (Bayesian NMF).
##
## Runs CoGAPS on the HLCA pericyte genes x cells matrix from
## 00.prepare_cogaps_input.py over a sweep of nPatterns x random seeds, so
## 02.select_rank.R can pick the rank *de novo* (cross-seed pattern
## robustness + reconstruction error), rather than fixing nP a priori to the
## curated program count. Uses single-cell distributed CoGAPS (cells partitioned
## across sets) because there are many cells.
##
## CoGAPS lives in the project .Rlib (Bioconductor; not in shared R_env).
##
## Outputs per (nPatterns N, seed tag TAG):
##   cogaps_result_npN{TAG}.rds        (full CogapsResult)
##   patterns_npN{TAG}.tsv.gz          (cell x pattern weights; sampleFactors)
##   feature_loadings_npN{TAG}.tsv.gz  (gene x pattern; featureLoadings)
##   pattern_markers_npN{TAG}.tsv.gz   (PatternMarkers genes per pattern)
##   cogaps_meta_npN{TAG}.tsv          (one row: np, seed, meanChiSq -- for
##                                      de-novo rank selection; per-run file so
##                                      concurrent array tasks don't race)
suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(CoGAPS)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INDIR     <- parse_arg("--indir", "../_m")
OUTDIR    <- parse_arg("--outdir", "../_m")
NPATTERNS <- as.integer(strsplit(parse_arg("--npatterns", "5,7,9"), ",")[[1]])
NITER     <- as.integer(parse_arg("--niterations", "5000"))
NSETS     <- as.integer(parse_arg("--nsets", "8"))
NTHREADS  <- as.integer(parse_arg("--nthreads", "1"))
SEED      <- as.integer(parse_arg("--seed", "13"))
## optional filename suffix so re-runs at other seeds don't clobber the canonical
## (seed-13) outputs; "" keeps the existing names (cogaps_result_np5.rds, ...).
TAG       <- parse_arg("--tag", "")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
fn <- function(stem, np, ext) file.path(
    OUTDIR, sprintf("%s_np%d%s.%s", stem, np, TAG, ext))

## ---- Load genes x cells matrix --------------------------------------------
mat <- readMM(file.path(INDIR, "cogaps_input.mtx"))
genes <- readLines(file.path(INDIR, "cogaps_genes.tsv"))
cells <- readLines(file.path(INDIR, "cogaps_barcodes.tsv"))
mat <- as.matrix(mat)
rownames(mat) <- genes
colnames(mat) <- cells
cat(sprintf("matrix: %d genes x %d cells\n", nrow(mat), ncol(mat)))

run_one <- function(np) {
    cat(sprintf("\n==== CoGAPS nPatterns=%d (nIter=%d, nSets=%d) ====\n",
                np, NITER, NSETS))
    params <- CogapsParams(nPatterns = np, nIterations = NITER, seed = SEED,
                           sparseOptimization = TRUE)
    ## single-cell distributed: partition cells across sets for scalability.
    params <- setDistributedParams(params, nSets = NSETS)
    res <- tryCatch(
        CoGAPS(mat, params, distributed = "single-cell",
               nThreads = NTHREADS, messages = FALSE),
        error = function(e) {
            cat("distributed CoGAPS failed (", conditionMessage(e),
                "); falling back to standard\n")
            CoGAPS(mat, CogapsParams(nPatterns = np, nIterations = NITER,
                                     seed = SEED, sparseOptimization = TRUE),
                   nThreads = NTHREADS, messages = FALSE)
        })

    saveRDS(res, fn("cogaps_result", np, "rds"))

    P <- res@sampleFactors            # cells x patterns
    A <- res@featureLoadings          # genes x patterns
    colnames(P) <- colnames(A) <- paste0("Pattern_", seq_len(ncol(P)))
    data.table::fwrite(
        data.frame(barcode = rownames(P), P),
        fn("patterns", np, "tsv.gz"), sep = "\t")
    data.table::fwrite(
        data.frame(gene = rownames(A), A),
        fn("feature_loadings", np, "tsv.gz"), sep = "\t")

    pm <- tryCatch(patternMarkers(res, threshold = "all"),
                   error = function(e) NULL)
    if (!is.null(pm)) {
        long <- do.call(rbind, lapply(names(pm$PatternMarkers), function(p)
            data.frame(pattern = p, gene = pm$PatternMarkers[[p]],
                       stringsAsFactors = FALSE)))
        data.table::fwrite(
            long, fn("pattern_markers", np, "tsv.gz"), sep = "\t")
    }
    ## one-row metadata for de-novo rank selection (read by 04, no rds reload).
    data.table::fwrite(
        data.frame(np = np,
                   seed = SEED,
                   tag = if (TAG == "") "seed13" else sub("^_", "", TAG),
                   niter = NITER, nsets = NSETS,
                   meanChiSq = res@metadata$meanChiSq,
                   n_genes = nrow(mat), n_cells = ncol(mat)),
        fn("cogaps_meta", np, "tsv"), sep = "\t")
    cat(sprintf("nP=%d done: meanChiSq=%.1f\n", np, res@metadata$meanChiSq))
}

for (np in NPATTERNS) run_one(np)

cat("\n---- sessionInfo ----\n")
sessioninfo::session_info()
