## =============================================================================
## Independent evaluation of the AGTR1 finding in GSE136831 (Adams/Kaminski).
##
## WHAT THIS IS FOR. The HLCA analysis (disease_association/_h/05) reports that
## disease-associated AGTR1 variation is carried by FIBROBLAST populations rather
## than by pericytes. GSE136831 is the only dataset in this project with a real
## COPD arm alongside Control and IPF, and it is fully independent of HLCA, so it
## is the natural place to ask whether that fibroblast pattern reproduces.
##
## WHAT IT CANNOT DO, stated first because it bounds every claim below. At a
## 5-cell-per-donor floor GSE136831 has exactly ONE Control donor with >= 5
## pericytes (against 6 COPD and 15 IPF). A pericyte-specific contrast is
## therefore NOT ESTIMABLE in either direction. This dataset can corroborate or
## fail to corroborate the FIBROBLAST half of the claim; it cannot establish or
## refute pericyte-specific AGTR1 regulation, and no output here should be read
## as doing so. That asymmetry is written into the outputs rather than left to
## the caption.
##
## Design:
##   UNIT      : donor x compartment pseudobulk (sum raw counts within the unit,
##               CP10K, log1p) from basement_membrane/_h/05.bm_copd.py run
##               against the RAS gene panel.
##   PRIMARY   : AGTR1, Control-vs-COPD and Control-vs-IPF, in the FIBROBLAST
##               and MYOFIBROBLAST compartments -- the fibroblast-lineage
##               compartments named by the HLCA result.
##
##               REPORTED AT NOMINAL ALPHA, NOT BH-CORRECTED (changed
##               2026-07-28). This is a directional REPLICATION test of a
##               hypothesis that was fixed before these data were touched -- the
##               compartments come from the HLCA analysis in
##               disease_association/_h/05, and GSE136831 is independent of HLCA.
##               A multiplicity penalty models a search over many candidate
##               answers; there is no search here. Worse, the family it would be
##               applied over is the wrong one: the two myofibroblast tests are
##               the SAME compartment evaluated in two disease arms that share
##               the same Control donors, so they are strongly correlated and BH
##               treats them as if they were independent looks. `p_BH_reference`
##               is still emitted so a reader who wants the corrected value can
##               see it, but it is a reference column, NOT the reporting gate.
##               The protection against fishing is the pre-specification plus the
##               `family` column, which keeps every non-pre-specified compartment
##               and gene labelled exploratory.
##   EXPLORATORY: the same model in every other powered compartment, and the
##               remaining RAS genes (AGTR2/AGT/ACE/ACE2), labelled exploratory
##               in the output rather than mixed into the primary family.
##   DESCRIPTIVE: pericytes -- group means, donor counts and a minimum
##               detectable effect, no p-value.
##   MODEL     : lm(y ~ disease + mean_log10_counts + sex + age + ever_smoker).
##               Single study, so no random effect -- but `ever_smoker` IS
##               recorded here, which HLCA does not permit, and it matters for a
##               COPD contrast.
## =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(data.table); library(dplyr); library(emmeans)
})

PRIMARY_GENE        <- "AGTR1"
PRIMARY_COMPARTMENTS <- c("Fibroblast", "Myofibroblast")
POWERED_MIN_DONORS  <- 5L

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--outdir", type = "character", default = "stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) fwrite(as.data.table(x), file.path(opt$outdir, f), sep = "\t")

pb <- fread(opt$pseudobulk)
pb <- pb[n_cells >= opt$min_cells]
pb[, disease := factor(disease, levels = c("Control", "COPD", "IPF"))]
pb[, `:=`(ever_smoker = factor(ever_smoker), sex = factor(sex),
          age = suppressWarnings(as.numeric(age)))]
genes <- sub("__expr$", "", grep("__expr$", names(pb), value = TRUE))
cat("genes in pseudobulk:", paste(genes, collapse = ", "), "\n")
stopifnot(PRIMARY_GENE %in% genes)

## ------------------------------------------------------- power inventory ----
power <- dcast(pb[, .(n_donors = uniqueN(donor_id)), by = .(compartment, disease)],
               compartment ~ disease, value.var = "n_donors", fill = 0L)
wt(power, "agtr1_copd_power_inventory.tsv")
cat("\n== donors per compartment (>=", opt$min_cells, " cells) ==\n"); print(power)

powered_copd <- power[Control >= POWERED_MIN_DONORS & COPD >= POWERED_MIN_DONORS, compartment]
powered_ipf  <- power[Control >= POWERED_MIN_DONORS & IPF  >= POWERED_MIN_DONORS, compartment]
cat("\npowered for COPD vs Control:", paste(powered_copd, collapse = ", "), "\n")
cat("powered for IPF  vs Control:", paste(powered_ipf,  collapse = ", "), "\n")
if (!"Pericyte" %in% powered_copd)
    cat("\nNOTE: Pericyte is NOT powered here -- descriptive only, no p-value.\n")

## ------------------------------------------------------------- modelling ----
fit_one <- function(gene, comp) {
    col <- paste0(gene, "__expr")
    d <- pb[compartment == comp & !is.na(get(col)) & !is.na(disease)]
    d[, y := get(col)]
    d[, disease := droplevels(disease)]
    if (uniqueN(d$disease) < 2) return(NULL)
    ## Drop any covariate that is constant or partly missing rather than failing
    ## the whole compartment.
    covs <- c("mean_log10_counts", "sex", "age", "ever_smoker")
    covs <- covs[vapply(covs, function(cc)
        cc %in% names(d) && sum(!is.na(d[[cc]])) == nrow(d) &&
            uniqueN(d[[cc]]) > 1, logical(1))]
    fit <- try(lm(reformulate(c("disease", covs), "y"), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    emm <- emmeans(fit, specs = "disease")
    ctr <- as.data.table(as.data.frame(pairs(emm, adjust = "none")))
    ctr <- ctr[grepl("Control", contrast)]
    setnames(ctr, "SE", "se")
    ## Sign flip so a POSITIVE estimate means "higher in disease". emmeans emits
    ## `Control - COPD`; leaving that convention in place would invert the
    ## direction of every statement in the figure caption.
    ctr[, `:=`(contrast = sub("^Control - ", "", contrast),
               estimate = -estimate)]
    ctr[, `:=`(gene = gene, compartment = comp, n_donors = nrow(d),
               ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se,
               covariates = paste(covs, collapse = "+"))]
    ns <- d[, .N, by = disease]
    ctr[, n_control := ns[disease == "Control", N][1]]
    ctr[]
}

comps_all <- power$compartment
grid <- CJ(gene = genes, compartment = comps_all, sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i)
    fit_one(grid$gene[i], grid$compartment[i])), fill = TRUE)

## A contrast is only reported as a TEST where both arms clear the donor floor.
res <- merge(res, melt(power, id.vars = "compartment", variable.name = "disease",
                       value.name = "n_donors_arm")[disease != "Control"],
             by.x = c("compartment", "contrast"), by.y = c("compartment", "disease"),
             all.x = TRUE)
res <- merge(res, power[, .(compartment, n_control_donors = Control)],
             by = "compartment", all.x = TRUE)
res[, estimable := n_donors_arm >= POWERED_MIN_DONORS & n_control_donors >= POWERED_MIN_DONORS]

res[, family := fifelse(gene == PRIMARY_GENE & compartment %in% PRIMARY_COMPARTMENTS &
                        estimable, "primary", "exploratory")]
## Reference only -- see the header. The primary family is reported at nominal
## alpha because this is a pre-specified replication, not a screen.
res[family == "primary", p_BH_reference := p.adjust(p.value, method = "BH")]

setcolorder(res, c("gene", "compartment", "contrast", "family", "estimable",
                   "estimate", "se", "ci_lo", "ci_hi", "df", "t.ratio", "p.value",
                   "p_BH_reference"))
res <- res[order(family, gene != PRIMARY_GENE, p.value)]
wt(res, "agtr1_copd_all_contrasts.tsv")
wt(res[family == "primary"], "agtr1_copd_primary.tsv")

cat("\n== PRIMARY family: AGTR1 in fibroblast-lineage compartments (nominal alpha) ==\n")
print(res[family == "primary",
          .(compartment, contrast, estimate, ci_lo, ci_hi, p.value, p_BH_reference,
            n_donors)])
## Concordance across the two disease arms is the replication evidence that a
## per-test p-value does not capture: an effect present in COPD AND IPF, in the
## same compartment and the same direction, is a stronger signal than either
## test alone.
conc <- dcast(res[family == "primary"], compartment ~ contrast, value.var = "estimate")
if (all(c("COPD", "IPF") %in% names(conc)))
    conc[, same_direction := sign(COPD) == sign(IPF)]
cat("\n== direction concordance across disease arms ==\n"); print(conc)
cat("\n== AGTR1, every compartment (exploratory beyond the primary family) ==\n")
print(res[gene == PRIMARY_GENE,
          .(compartment, contrast, estimable, estimate, ci_lo, ci_hi, p.value, n_donors)])

## ------------------------------------------- descriptive pericyte report ----
## No p-value. Group means, donor counts, and the effect this dataset COULD have
## detected at 80% power given its actual n and residual SD, so "we saw nothing"
## is separable from "we could not have seen anything".
desc <- pb[, .(n_donors = uniqueN(donor_id), n_cells = sum(n_cells),
               mean_AGTR1 = mean(AGTR1__expr, na.rm = TRUE),
               sd_AGTR1 = sd(AGTR1__expr, na.rm = TRUE),
               detect_AGTR1 = mean(AGTR1__detect, na.rm = TRUE)),
           by = .(compartment, disease)][order(compartment, disease)]
wt(desc, "agtr1_copd_descriptive.tsv")
cat("\n== descriptive AGTR1 by compartment x disease ==\n"); print(desc)

mde <- rbindlist(lapply(comps_all, function(comp) {
    d <- pb[compartment == comp & !is.na(AGTR1__expr)]
    s <- d[, sd(AGTR1__expr, na.rm = TRUE)]
    rbindlist(lapply(c("COPD", "IPF"), function(arm) {
        n1 <- d[disease == "Control", .N]; n2 <- d[disease == arm, .N]
        ok <- n1 >= 2 && n2 >= 2 && is.finite(s) && s > 0
        data.table(compartment = comp, contrast = arm, n_control = n1, n_disease = n2,
                   sd_log1p_cp10k = s,
                   ## two-sample, alpha 0.05 two-sided, 80% power
                   mde_80pct = if (ok) 2.802 * s * sqrt(1 / n1 + 1 / n2) else NA_real_,
                   estimable = ok && n1 >= POWERED_MIN_DONORS && n2 >= POWERED_MIN_DONORS)
    }))
}))
wt(mde, "agtr1_copd_mde.tsv")
cat("\n== minimum detectable AGTR1 effect (80% power) ==\n"); print(mde)

readme <- c(
    "GSE136831 (Adams/Kaminski) independent AGTR1 evaluation -- generated summary",
    sprintf("Per-donor cell floor: %d", opt$min_cells),
    "",
    "PRIMARY family: AGTR1 x {Fibroblast, Myofibroblast} x {COPD, IPF} vs Control.",
    "Reported at NOMINAL alpha: these compartments were fixed by the independent",
    "HLCA analysis (disease_association/_h/05) before these data were touched, so",
    "this is a directional replication, not a screen. p_BH_reference is provided",
    "for reference only and is NOT the reporting gate. Everything else is",
    "exploratory and labelled so in the `family` column.",
    "",
    "Estimates are signed so that POSITIVE = higher in disease than in Control.",
    "",
    "PERICYTES ARE NOT ESTIMABLE in this dataset: see agtr1_copd_power_inventory.tsv",
    "and agtr1_copd_mde.tsv. Only one Control donor clears the cell floor, so no",
    "pericyte contrast is reported and none should be inferred. This dataset",
    "evaluates the FIBROBLAST half of the HLCA result only.",
    "",
    "Powered for COPD vs Control: " , paste(powered_copd, collapse = ", "),
    "Powered for IPF vs Control: "  , paste(powered_ipf,  collapse = ", "))
writeLines(readme, file.path(opt$outdir, "agtr1_copd_README.txt"))

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
