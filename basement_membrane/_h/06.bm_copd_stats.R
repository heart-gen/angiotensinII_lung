## COPD basement-membrane contrast in GSE136831 (Adams/Kaminski), bounded honestly.
##
## Literature reports LAMB1 and LAMA4 dysregulation in COPD. LAMB1 and LAMA4 are
## therefore PRE-SPECIFIED as the primary hypotheses; the remaining BM genes are
## exploratory and labelled as such in every output.
##
## THE POWER CEILING, stated first because it bounds every claim below: at a
## 5-cell-per-donor floor this dataset has 6 COPD but only ONE Control donor with
## >=5 pericytes. A pericyte-specific COPD-vs-Control test is NOT ESTIMABLE. The
## pericyte compartment is therefore reported descriptively, with a formal
## minimum detectable effect, and NO p-value. Compartments carrying the
## confirmatory tests are those with >=5 donors on both arms (Endothelial,
## Fibroblast, Myofibroblast; Mural is marginal and flagged exploratory).
##
## IPF-vs-Control is reported alongside as an internal positive control: IPF is
## well powered here, so if the pipeline detects nothing in IPF either, the assay
## rather than the disease is in question.
##
## Single study, so no (1|study) term -- but `ever_smoker` IS available here and
## is adjusted for, which HLCA does not permit.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(data.table)
    library(emmeans)
})

PRIMARY_GENES <- c("LAMB1", "LAMA4")
POWERED_MIN_DONORS <- 5L

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--panels", type = "character"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(x, file)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)

pb <- fread(opt$pseudobulk)
panels <- fread(opt$panels)
bm_genes <- panels[panel == "basement_membrane", unique(gene)]
bm_genes <- intersect(bm_genes,
                      sub("__expr$", "", grep("__expr$", names(pb), value = TRUE)))

pb <- pb[n_cells >= opt$min_cells]
pb[, disease := factor(disease, levels = c("Control", "COPD", "IPF"))]
pb[, ever_smoker := factor(ever_smoker)]
pb[, sex := factor(sex)]
pb[, age := suppressWarnings(as.numeric(age))]

## ------------------------------------------------------- power inventory ----
power <- pb[, .(n_donors = uniqueN(donor_id), n_cells = sum(n_cells)),
            by = .(compartment, disease)]
power_wide <- dcast(power, compartment ~ disease, value.var = "n_donors",
                    fill = 0L)
write_tsv_safe(power_wide, file.path(opt$outdir, "bm_copd_power_inventory.tsv"))
message("Donors per compartment (>= ", opt$min_cells, " cells):")
print(power_wide)

powered <- power_wide[Control >= POWERED_MIN_DONORS & COPD >= POWERED_MIN_DONORS,
                      compartment]
message("Powered compartments (>=", POWERED_MIN_DONORS,
        " donors on both arms): ", paste(powered, collapse = ", "))
if (!length(powered))
    warning("no compartment is powered for COPD vs Control")

## ------------------------------------------------------------- modelling ----
fit_one <- function(dt, gene, comp) {
    col <- paste0(gene, "__expr")
    d <- dt[compartment == comp & !is.na(get(col))]
    d[, y := get(col)]
    d <- d[!is.na(disease)]
    d[, disease := droplevels(disease)]
    if (uniqueN(d$disease) < 2) return(NULL)
    ## Adjust for depth and demographics where available; drop any covariate
    ## that is constant or fully missing rather than failing the whole fit.
    covs <- c("mean_log10_counts", "sex", "age", "ever_smoker")
    covs <- covs[vapply(covs, function(cc)
        cc %in% names(d) && sum(!is.na(d[[cc]])) == nrow(d) &&
            length(unique(d[[cc]])) > 1, logical(1))]
    fit <- try(lm(reformulate(c("disease", covs), "y"), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    emm <- emmeans(fit, specs = "disease")
    ctr <- as.data.frame(pairs(emm, adjust = "none"))
    ctr <- ctr[grepl("Control", ctr$contrast), , drop = FALSE]
    if (!nrow(ctr)) return(NULL)
    data.frame(compartment = comp, gene = gene, ctr,
               n_donors = nrow(d),
               covariates = paste(covs, collapse = "+"),
               row.names = NULL)
}

run_family <- function(genes, comps, label) {
    res <- rbindlist(
        lapply(comps, function(cm)
            rbindlist(lapply(genes, function(g) fit_one(pb, g, cm)), fill = TRUE)),
        fill = TRUE)
    if (!nrow(res)) return(res)
    ## BH within the family, and only over the COPD-vs-Control tests: the
    ## IPF contrast is a positive control, not part of the hypothesis family.
    res[, family := label]
    res[, is_copd_contrast := grepl("COPD", contrast)]
    res[is_copd_contrast == TRUE,
        p_BH := p.adjust(p.value, method = "BH")]
    res[]
}

primary <- run_family(intersect(PRIMARY_GENES, bm_genes), powered, "primary")
exploratory <- run_family(setdiff(bm_genes, PRIMARY_GENES), powered, "exploratory")

write_tsv_safe(primary, file.path(opt$outdir, "bm_copd_primary.tsv"))
write_tsv_safe(exploratory, file.path(opt$outdir, "bm_copd_exploratory.tsv"))
message("Primary family: ", nrow(primary[is_copd_contrast == TRUE]),
        " COPD-vs-Control tests across ", length(powered), " compartments")

## ------------------------------ pericytes: descriptive + minimum detectable ----
## No test is run here. The MDE converts "we found nothing" into the quantitative
## statement "this design could not have detected anything smaller than X".
mde <- function(n1, n2, sd_pooled, power = 0.80, alpha = 0.05) {
    if (any(is.na(c(n1, n2, sd_pooled))) || n1 < 2 || n2 < 2 || !is.finite(sd_pooled))
        return(NA_real_)
    df <- n1 + n2 - 2
    (qt(1 - alpha / 2, df) + qt(power, df)) * sd_pooled * sqrt(1 / n1 + 1 / n2)
}

peri_rows <- lapply(bm_genes, function(g) {
    col <- paste0(g, "__expr")
    d <- pb[compartment == "Pericyte" & !is.na(get(col))]
    if (!nrow(d)) return(NULL)
    n_ctrl <- d[disease == "Control", uniqueN(donor_id)]
    n_copd <- d[disease == "COPD", uniqueN(donor_id)]
    ## Between-donor SD is taken from the best-populated arm (IPF), because the
    ## COPD/Control arms are too small to estimate it themselves.
    sd_ref <- d[disease == "IPF", sd(get(col), na.rm = TRUE)]
    data.frame(
        gene = g, n_donors_control = n_ctrl, n_donors_copd = n_copd,
        n_donors_ipf = d[disease == "IPF", uniqueN(donor_id)],
        mean_control = d[disease == "Control", mean(get(col), na.rm = TRUE)],
        mean_copd = d[disease == "COPD", mean(get(col), na.rm = TRUE)],
        mean_ipf = d[disease == "IPF", mean(get(col), na.rm = TRUE)],
        sd_between_donor_ipf = sd_ref,
        mde_log1p_cp10k = mde(n_ctrl, n_copd, sd_ref),
        estimable = n_ctrl >= 2 && n_copd >= 2,
        row.names = NULL)
})
peri <- rbindlist(Filter(Negate(is.null), peri_rows), fill = TRUE)
write_tsv_safe(peri, file.path(opt$outdir, "bm_pericyte_power.tsv"))

readme <- c(
    "COPD basement-membrane contrast (GSE136831) -- generated summary",
    sprintf("Per-donor cell floor: %d", opt$min_cells),
    "",
    "Donors per compartment x disease:",
    paste(utils::capture.output(print(power_wide)), collapse = "\n"),
    "",
    sprintf("Powered compartments (>=%d donors both arms): %s",
            POWERED_MIN_DONORS, paste(powered, collapse = ", ")),
    sprintf("Primary family: %s x %d compartments = %d COPD-vs-Control tests",
            paste(intersect(PRIMARY_GENES, bm_genes), collapse = "+"),
            length(powered), nrow(primary[is_copd_contrast == TRUE])),
    "",
    "PERICYTES: NOT ESTIMABLE for COPD vs Control -- descriptive only, no p-value.",
    paste(utils::capture.output(print(
        peri[, .(gene, n_donors_control, n_donors_copd, mde_log1p_cp10k, estimable)])),
        collapse = "\n"))
writeLines(readme, file.path(opt$outdir, "bm_copd_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
