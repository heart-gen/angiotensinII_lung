## Subset lung single-cell data
                                        # Load required libraries
suppressPackageStartupMessages({
    library(here)
    library(scater)
    library(scuttle)
    library(harmony)
    library(SingleCellExperiment)
})

                                        # Set zellkonverter parameters
Sys.setenv(ZELLKONVERTER_USE_BASILISK = "FALSE")
reticulate::use_python("/ocean/projects/bio250020p/shared/opt/env/R_env/bin/python",
                       required = TRUE)

                                        # Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript script.R <model>")
}

model <- args[1]
cat("Model selected:", model, "\n")

                                        # Functions
subset_data <- function(input_file, COMPARTMENT = FALSE) {
					# Load data
    sce <- zellkonverter::readH5AD(input_file)
    if ("soupX" %in% names(assays(sce))) {
        names(assays(sce)) <- c("counts", "soupX")
    } else {
        names(assays(sce)) <- c("counts")
    }

                                        # Remove missing or unknown
                                        # annotations (level 4)
    sce <- sce[, !is.na(colData(sce)$ann_level_4)]
    sce <- sce[, !colData(sce)$ann_level_4 %in% c("None", "Unknown")]

                                        # Update annotation columns
    colData(sce)$subclusters <- sce$ann_finest_level
    colData(sce)$cell_type   <- sce$ann_level_4
    colData(sce)$compartment <- sce$ann_level_1
    colData(sce)$clusters    <- sce$cell_type
    colData(sce)$patient     <- sce$donor_id
    colLabels(sce) <- sce$cell_type

                                        # Remove studies with fewer
                                        # than 20 pericytes
    if ("study" %in% colnames(colData(sce))) {
        study_counts  <- table(colData(sce)$study)
        valid_studies <- names(study_counts[study_counts >= 20])
        sce <- sce[, colData(sce)$study %in% valid_studies]
    } else {
        warning("Variable 'study' not found in colData; skipping study filter.")
    }

    # Subset data
    if (COMPARTMENT) {
        sce_sub <- sce[, colData(sce)$compartment %in% c("Stroma")]
    } else {
        sce_sub <- sce[, colData(sce)$subclusters %in% c("Pericytes")]
    }
    return(sce_sub)
}

preprocess_data <- function(sce) {    
                                        # Check for batch variables
    batch_vars <- c(
        "data", "assay", "tissue_sampling_method", "sequencing_platform",
        "development_stage", "tissue", "subject_type", "study",
        "lung_condition", "sex", "self_reported_ethnicity", "age_or_mean_of_age_range"
    )
    present_vars <- batch_vars[batch_vars %in% colnames(colData(sce))]
    missing_vars <- setdiff(batch_vars, present_vars)
    if (length(missing_vars) > 0) {
        warning(paste("Missing batch variables:",
                      paste(missing_vars, collapse = ", ")))
    }

                                        # Preprocessing
    sce <- scuttle::logNormCounts(sce)

                                        # Run Harmony batch correction
    if (length(present_vars) > 0) {
        sce <- scater::runPCA(sce)
        harmony_df <- as.data.frame(colData(sce)[, present_vars, drop = FALSE])
        harmony_embeddings <- HarmonyMatrix(
            reducedDim(sce, "PCA"),
            harmony_df,
            vars_use = present_vars
        )
        reducedDim(sce, "HARMONY") <- harmony_embeddings
    } else {
        warning("No batch variables found for Harmony correction.")
    }
    return(sce)
}

#### Main execution loop

for (COMPARTMENT in c(FALSE, TRUE)) {
    label    <- ifelse(COMPARTMENT, "stroma", "pericyte")
    out_file <- paste0(label, ".hlca_", model, ".dataset.h5ad")
    in_file  <- here("inputs/hlca/_m", paste0("hlca_", model, ".h5ad"))
                                        # Call function
    sce <- subset_data(in_file, COMPARTMENT)
    sce <- preprocess_data(sce)
                                        # Write as H5AD
    zellkonverter::writeH5AD(sce, file = out_file)
}
    
#### Reproducibility information ####
cat("Reproducibility information:\n")
Sys.time()
proc.time()
options(width = 120)
reticulate::py_config()
sessioninfo::session_info()
