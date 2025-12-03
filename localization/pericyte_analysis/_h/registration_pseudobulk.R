suppressPackageStartupMessages({
    library(edgeR)
    library(scuttle)
    library(SingleCellExperiment)
})

registration_pseudobulk <-
    function(sce,
             var_registration,
             var_sample_id,
             covars = NULL,
             assay_type = "counts",
             min_ncells = 10,
             pseudobulk_rds_file = NULL) {
        ## Check that inputs are correct
        stopifnot(is(sce, "SingleCellExperiment"))
        stopifnot(var_registration %in% colnames(colData(sce)))
        stopifnot(var_sample_id %in% colnames(colData(sce)))
        stopifnot(all(
            !c("registration_sample_id", "registration_variable") %in%
            colnames(colData(sce))
        ))

        ## Avoid any incorrect inputs that are otherwise hard to detect
        stopifnot(!var_registration %in% covars)
        stopifnot(!var_sample_id %in% covars)
        stopifnot(var_registration != var_sample_id)

        ## Check that the values in the registration variable are ok
        uniq_var_regis <- unique(sce[[var_registration]])
        if (any(grepl("\\+|\\-", uniq_var_regis))) {
            stop(
                "Remove the + and - signs in colData(sce)[, '",
                var_registration,
                "'] to avoid downstream issues.",
                call. = FALSE
            )
        }

        ## Pseudo-bulk for our current BayesSpace cluster results
        message(Sys.time(), " make pseudobulk object")
        ## I think this needs counts assay
        sce_pseudo <- scuttle::aggregateAcrossCells(
            sce,
            DataFrame(
                registration_variable = sce[[var_registration]],
                registration_sample_id = sce[[var_sample_id]]
            ),
            use.assay.type = assay_type
        )
        colnames(sce_pseudo) <-
            paste0(
                sce_pseudo$registration_sample_id,
                "_",
                sce_pseudo$registration_variable
            )

        ## Check that the covariates are present
        if (!is.null(covars)) {
            for (covariate_i in covars) {
                if (sum(is.na(sce_pseudo[[covariate_i]])) == ncol(sce_pseudo)) {
                    stop(
                        "Covariate '",
                        covariate_i,
                        "' has all NAs after pseudo-bulking. Might be due to not being a sample-level covariate.",
                        call. = FALSE
                    )
                }
            }
        }

        ## Drop pseudo-bulked samples that had low initial contribution
        ## of raw-samples. That is, pseudo-bulked samples that are not
        ## benefiting from the pseudo-bulking process to obtain higher counts.
        if (!is.null(min_ncells)) {
            message(
                Sys.time(),
                " dropping ",
                sum(sce_pseudo$ncells < min_ncells),
                " pseudo-bulked samples that are below 'min_ncells'."
            )
            sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= min_ncells]
        }

        if (is.factor(sce_pseudo$registration_variable)) {
            ## Drop unused var_registration levels if we had to drop some due
            ## to min_ncells
            sce_pseudo$registration_variable <- droplevels(sce_pseudo$registration_variable)
        }

        ## Drop lowly-expressed genes
        message(Sys.time(), " drop lowly expressed genes")
        keep_expr <-
            edgeR::filterByExpr(sce_pseudo, group = sce_pseudo$registration_variable)
        sce_pseudo <- sce_pseudo[which(keep_expr), ]

        ## Compute the logcounts
        message(Sys.time(), " normalize expression")
        logcounts(sce_pseudo) <-
            edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
                log = TRUE,
                prior.count = 1
            )

        if (is(sce_pseudo, "SpatialExperiment")) {
            ## Drop things we don't need
            spatialCoords(sce_pseudo) <- NULL
            imgData(sce_pseudo) <- NULL
        }
        if (!is.null(pseudobulk_rds_file)) {
            message(Sys.time(), " saving sce_pseudo to ", pseudobulk_rds_file)
            saveRDS(sce_pseudo, file = pseudobulk_rds_file)
        }

        ## Done!
        return(sce_pseudo)
    }
