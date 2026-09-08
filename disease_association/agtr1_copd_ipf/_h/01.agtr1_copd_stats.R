## =============================================================================
## Independent evaluation of the AGTR1 finding in GSE136831 (Adams/Kaminski).
##
## WHAT THIS IS FOR. The HLCA analysis (disease_association/_h/05) reports that
## disease-associated AGTR1 variation is carried by FIBROBLAST populations rather
## than by pericytes. GSE136831 is the only dataset in this project with a real
## COPD arm alongside Control and IPF, and it is fully independent of HLCA, so it
## is the natural place to ask whether that fibroblast pattern reproduces.
##
## This is an INDEPENDENT EVALUATION, not a replication, and the word matters.
## Nothing in the HLCA cell-type analysis is significant (min BH 0.199 on the
## contrasts, 0.332 on the omnibus), so there is no directional finding to
## replicate: the correct statement of a disagreement is "no directional finding
## to replicate", never "the replication failed". The cross-cohort sign table
## emitted below reports signs only, with the scale difference stated, and never
## differences the two estimates.
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
##               and MYOFIBROBLAST compartments.
##
##               ** HOW THE FAMILY WAS ACTUALLY CHOSEN (corrected 2026-09-07,
##               defect P1-19). ** This header used to say the two compartments
##               were "named by the HLCA result". That is true of FIBROBLAST and
##               false of MYOFIBROBLAST, and the distinction matters because
##               myofibroblast is the compartment carrying both significant
##               results. Timeline: the family was fixed here on 2026-07-28; the
##               HLCA analysis (disease_association/_h/05) was rebuilt on
##               2026-07-30, when the `age > 20` gate was fixed and myofibroblasts
##               were admitted to it for the first time (1 -> 10 fibrotic donors).
##               So when this family was pre-specified, HLCA had NO myofibroblast
##               estimate to name.
##
##               Myofibroblast was chosen on LINEAGE grounds -- a perfectly good
##               pre-specification, and the nominal-alpha argument below still
##               holds, because the family was still fixed before these data were
##               touched. But it is a different pre-specification from the one
##               previously claimed, and it must be described as the one it is.
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
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    ## The HLCA side of the comparison. Optional: if absent the cross-cohort block
    ## is skipped with a message rather than failing the job.
    make_option("--hlca-effects", type = "character", dest = "hlca",
                default = "../../_m/mean_expr/agtr1_celltype_disease_effects.tsv")
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

## THESE ARE DONOR-COUNT GATES, NOT POWER CALCULATIONS (renamed 2026-09-07, P3-14).
## They were called `powered_*`, and the log and generated README announced
## "Powered for COPD vs Control: ...". That word is not earned: the criterion is
## `>= POWERED_MIN_DONORS` donors in each arm and nothing else. In four of the
## compartments it admits, AGTR1 is IDENTICALLY ZERO in a whole arm
## (agtr1_copd_descriptive.tsv: ATI COPD 0 and IPF 0; ATII Control 0;
## Endothelial COPD 0; SMC Control 0), with minimum detectable effects of
## 0.007-0.015 log units in a gene whose fibroblast group means are 0.17-0.24.
## Real power lives in agtr1_copd_mde.tsv; this is a floor.
has_donor_floor_copd <- power[Control >= POWERED_MIN_DONORS & COPD >= POWERED_MIN_DONORS, compartment]
has_donor_floor_ipf  <- power[Control >= POWERED_MIN_DONORS & IPF  >= POWERED_MIN_DONORS, compartment]
cat("\nmeets the >=", POWERED_MIN_DONORS, "-donor floor, COPD vs Control:",
    paste(has_donor_floor_copd, collapse = ", "), "\n")
cat("meets the >=", POWERED_MIN_DONORS, "-donor floor, IPF  vs Control:",
    paste(has_donor_floor_ipf,  collapse = ", "), "\n")
cat("  (a donor-count gate, NOT a power calculation -- see agtr1_copd_mde.tsv)\n")
if (!"Pericyte" %in% has_donor_floor_copd)
    cat("\nNOTE: Pericyte does NOT meet the donor floor -- descriptive only, no p-value.\n")

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
    ##
    ## The test statistic has to turn with the estimate. Flipping `estimate`
    ## alone left `t.ratio` describing the opposite comparison in every row of
    ## agtr1_copd_all_contrasts.tsv (68) and agtr1_copd_primary.tsv (4) --
    ## the same defect fixed in basement_membrane and agt_axis. Two-sided
    ## p-values, `se` and the CI bounds (recomputed from the flipped estimate
    ## below) were never affected.
    tstat <- intersect(c("t.ratio", "z.ratio"), names(ctr))
    if (length(tstat)) ctr[, (tstat) := lapply(.SD, function(x) -x),
                           .SDcols = tstat]
    ctr[, `:=`(contrast = sub("^Control - ", "", contrast),
               estimate = -estimate)]
    ## `n_model_donors`, not `n_donors` (renamed 2026-09-07, P3-16): this is the
    ## FULL model's donor count -- all three arms -- sitting beside a two-group
    ## contrast label like "IPF". The per-arm counts are `n_donors_arm` and
    ## `n_control_donors`, added below. Not a wrong number; a naming hazard, in
    ## the table most likely to be read directly.
    ctr[, `:=`(gene = gene, compartment = comp, n_model_donors = nrow(d),
               ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se,
               covariates = paste(covs, collapse = "+"))]
    ns <- d[, .N, by = disease]
    ctr[, n_control := ns[disease == "Control", N][1]]
    ctr[]
}

## ---- COVARIATE BALANCE + UNADJUSTED CONTRAST (P2-27, added 2026-09-07) ----
## The header presents `ever_smoker` as this module's advantage over HLCA, and
## having it IS better than not having it. But it is close to a perfect separator
## of the COPD arm, so the COPD coefficient is not "COPD adjusted for smoking" --
## it is COPD identified through a smoking slope borrowed from the IPF arm, under
## an untested no-interaction assumption. Two things are needed to let a reader
## see that, and neither existed: the balance itself, and the same contrast
## without the covariate.
balance_one <- function(comp) {
    d <- pb[compartment == comp & !is.na(disease)]
    if (!nrow(d)) return(NULL)
    d[, disease := droplevels(disease)]
    d[, .(n_donors = .N,
          n_ever_smoker = sum(tolower(as.character(ever_smoker)) %in%
                                  c("1", "true", "y", "yes"), na.rm = TRUE),
          pct_ever_smoker = 100 * mean(tolower(as.character(ever_smoker)) %in%
                                           c("1", "true", "y", "yes"), na.rm = TRUE),
          n_female = sum(tolower(as.character(sex)) %in% c("f", "female"), na.rm = TRUE),
          age_median = suppressWarnings(median(as.numeric(as.character(age)), na.rm = TRUE)),
          mean_log10_counts = mean(mean_log10_counts, na.rm = TRUE)),
      by = disease][, compartment := comp][]
}
balance <- rbindlist(lapply(unique(pb$compartment), balance_one), fill = TRUE)
if (nrow(balance)) {
    ## A covariate that takes one value inside an arm cannot be adjusted for
    ## within it. Flag it explicitly rather than leaving it to be noticed.
    balance[, ever_smoker_constant_in_arm := pct_ever_smoker %in% c(0, 100)]
    wt(balance[order(compartment, disease)], "agtr1_copd_covariate_balance.tsv")
    cat("\n== covariate balance by compartment and arm (P2-27) ==\n")
    print(balance[order(compartment, disease)])
}

## Same contrast under three covariate sets, so the adjustment's effect is
## visible instead of assumed. `adjusted` is the shipped primary; the other two
## are diagnostics and are written to their own file.
fit_arms <- function(gene, comp) {
    col <- paste0(gene, "__expr")
    d <- pb[compartment == comp & !is.na(get(col)) & !is.na(disease)]
    d[, y := get(col)]; d[, disease := droplevels(disease)]
    if (uniqueN(d$disease) < 2) return(NULL)
    keep <- function(cc) cc %in% names(d) && sum(!is.na(d[[cc]])) == nrow(d) &&
        uniqueN(d[[cc]]) > 1
    full <- Filter(keep, c("mean_log10_counts", "sex", "age", "ever_smoker"))
    sets <- list(adjusted   = full,
                 no_smoking = setdiff(full, "ever_smoker"),
                 unadjusted = character(0))
    rbindlist(lapply(names(sets), function(nm) {
        f <- try(lm(reformulate(c("disease", sets[[nm]]), "y"), data = d), silent = TRUE)
        if (inherits(f, "try-error")) return(NULL)
        ct <- as.data.table(as.data.frame(pairs(emmeans(f, specs = "disease"),
                                                adjust = "none")))
        ct <- ct[grepl("Control", contrast)]
        setnames(ct, "SE", "se")
        tstat <- intersect(c("t.ratio", "z.ratio"), names(ct))
        if (length(tstat)) ct[, (tstat) := lapply(.SD, function(x) -x), .SDcols = tstat]
        ct[, `:=`(contrast = sub("^Control - ", "", contrast), estimate = -estimate)]
        ct[, `:=`(gene = gene, compartment = comp, model_arm = nm,
                  covariates = if (length(sets[[nm]])) paste(sets[[nm]], collapse = "+") else "none",
                  n_model_donors = nrow(d),
                  ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)]
        ct[]
    }), fill = TRUE)
}

comps_all <- power$compartment
grid <- CJ(gene = genes, compartment = comps_all, sorted = FALSE)
res <- rbindlist(lapply(seq_len(nrow(grid)), function(i)
    fit_one(grid$gene[i], grid$compartment[i])), fill = TRUE)

arms <- rbindlist(lapply(seq_len(nrow(grid)), function(i)
    fit_arms(grid$gene[i], grid$compartment[i])), fill = TRUE)
if (nrow(arms)) {
    wt(arms[order(gene, compartment, contrast, model_arm)],
       "agtr1_copd_model_arms.tsv")
    cat("\n== same contrast, three covariate sets (P2-27) ==\n")
    print(dcast(arms[gene == "AGTR1"], compartment + contrast ~ model_arm,
                value.var = "estimate"))
}

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
            n_model_donors)])
## Concordance across the two disease arms is the replication evidence that a
## per-test p-value does not capture: an effect present in COPD AND IPF, in the
## same compartment and the same direction, is a stronger signal than either
## test alone.
conc <- dcast(res[family == "primary"], compartment ~ contrast, value.var = "estimate")
if (all(c("COPD", "IPF") %in% names(conc)))
    conc[, same_direction := sign(COPD) == sign(IPF)]
cat("\n== direction concordance across disease arms ==\n"); print(conc)
## ------------------------------------------------- cross-cohort sign table ----
## The computation this module was missing (defect P1-19). The in-script
## concordance check above compares GSE136831's COPD arm to its own IPF arm and
## never touches HLCA -- so a module described as a replication had no comparison
## to the thing it was replicating.
##
## THREE RULES, enforced here rather than left to a caption:
##  1. SIGNS ONLY. HLCA estimates are in within-cell-type SD units; these are in
##     log1p CP10K. The two are never differenced, and no ratio is formed.
##  2. HLCA HAS NO COPD ARM. `05.agtr1_celltype_disease.R` excludes COPD outright
##     (its TRI_LEVELS are Healthy / Fibrotic_ILD / Other), and its "Other" group
##     is a COVID-dominated grab-bag, NOT a COPD arm. So the GSE136831 COPD
##     contrast has nothing to compare against, and the row says so instead of
##     silently borrowing "Other".
##  3. NEITHER HLCA ESTIMATE IS SIGNIFICANT. `agreement` therefore reports sign
##     agreement against a NULL reference, and `hlca_is_null` marks it.
##
## Population correspondence is ASSERTED, NOT ESTABLISHED: HLCA
## `ann_finest_level` "Myofibroblasts" and GSE136831 `Manuscript_Identity`
## "Myofibroblast" are not mapped to each other anywhere in this repository, and
## HLCA splits fibroblasts into three subtypes where GSE136831 has one coarse
## class. `mapping_confidence` records that.
XCOHORT_MAP <- list(
    list(gse = "Myofibroblast", hlca = "Myofibroblasts",
         conf = "name-matched only; populations not established as equivalent"),
    list(gse = "Fibroblast",
         hlca = c("Alveolar fibroblasts", "Adventitial fibroblasts",
                  "Peribronchial fibroblasts"),
         conf = "one GSE136831 class vs three HLCA subtypes; range reported, not a mean"))

if (file.exists(opt$hlca)) {
    hl <- fread(opt$hlca)
    hl <- hl[response == "z_AGTR1" & contrast == "Fibrotic_ILD - Healthy"]
    xc <- rbindlist(lapply(XCOHORT_MAP, function(m) {
        rbindlist(lapply(c("IPF", "COPD"), function(arm) {
            g <- res[family == "primary" & compartment == m$gse & contrast == arm]
            if (nrow(g) != 1) return(NULL)
            h <- hl[cell_type %in% m$hlca]
            has_h <- arm == "IPF" && nrow(h) > 0
            data.table(
                gse_compartment  = m$gse,
                gse_contrast     = arm,
                gse_estimate_log1p_cp10k = g$estimate,
                gse_p            = g$p.value,
                gse_sign         = if (g$estimate > 0) "+" else "-",
                hlca_cell_types  = if (has_h) paste(h$cell_type, collapse = "; ") else NA_character_,
                hlca_contrast    = if (has_h) "Fibrotic_ILD - Healthy" else NA_character_,
                hlca_estimate_sd_units = if (has_h) paste(sprintf("%+.3f", h$estimate), collapse = "; ") else NA_character_,
                hlca_p           = if (has_h) paste(sprintf("%.3f", h$p.value), collapse = "; ") else NA_character_,
                hlca_sign        = if (!has_h) NA_character_
                                   else if (all(h$estimate > 0)) "+"
                                   else if (all(h$estimate < 0)) "-" else "mixed",
                hlca_is_null     = if (has_h) all(h$p.value >= 0.05) else NA,
                agreement = if (!has_h) {
                        "NO HLCA COMPARATOR -- HLCA/05 excludes COPD; its 'Other' arm is COVID-dominated, not COPD"
                    } else {
                        hs <- if (all(h$estimate > 0)) "+" else if (all(h$estimate < 0)) "-" else "mixed"
                        gs <- if (g$estimate > 0) "+" else "-"
                        ## The qualifier tracks `hlca_is_null` rather than being
                        ## asserted: one HLCA fibroblast subtype is nominally
                        ## significant (peribronchial, p = 0.030) even though
                        ## nothing in that analysis survives BH (min 0.199).
                        q <- if (all(h$p.value >= 0.05))
                                 "against a NON-SIGNIFICANT HLCA estimate -- no directional finding to replicate"
                             else
                                 "against an HLCA estimate significant at nominal alpha but NOT after BH (min BH in that analysis is 0.199)"
                        if (hs == "mixed")
                            paste0("HLCA subtypes disagree in sign among themselves (", 
                                   paste(sprintf("%+.3f", h$estimate), collapse = ", "), ")")
                        else if (hs == gs) paste0("same sign, ", q)
                        else paste0("OPPOSITE sign, ", q)
                    },
                mapping_confidence = m$conf,
                scale_note = "NOT COMPARABLE IN MAGNITUDE: HLCA is within-cell-type SD, GSE136831 is log1p CP10K. Signs only.")
        }))
    }))
    wt(xc, "agtr1_cross_cohort_signs.tsv")
    cat("\n== cross-cohort SIGN comparison vs HLCA (magnitudes are NOT comparable) ==\n")
    print(xc[, .(gse_compartment, gse_contrast, gse_sign, hlca_sign, hlca_is_null, agreement)])
} else {
    cat("\n== cross-cohort sign table SKIPPED: no HLCA effects file at ", opt$hlca,
        " ==\n", sep = "")
}

cat("\n== AGTR1, every compartment (exploratory beyond the primary family) ==\n")
print(res[gene == PRIMARY_GENE,
          .(compartment, contrast, estimable, estimate, ci_lo, ci_hi, p.value,
            n_model_donors)])

## ------------------------------------------- descriptive pericyte report ----
## No p-value. Group means, donor counts, and the effect this dataset COULD have
## detected at 80% power given its actual n, so "we saw nothing" is separable
## from "we could not have seen anything".
##
## THE SD IS MARGINAL, NOT RESIDUAL (label corrected 2026-09-07, P3-15). This
## header used to say "residual SD". `mde_80pct` below is built from
## `sd(AGTR1__expr)` across all donors in the compartment, before any model, so
## it includes between-group variance. Marginal >= residual, so the true MDE is
## SMALLER than reported and this table UNDER-claims the module's sensitivity --
## the approximation is conservative in the direction that matters. The column is
## named `sd_log1p_cp10k` and says so.
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
    paste0("Meets the >=", POWERED_MIN_DONORS,
           "-donor floor for COPD vs Control: "), paste(has_donor_floor_copd, collapse = ", "),
    paste0("Meets the >=", POWERED_MIN_DONORS,
           "-donor floor for IPF vs Control: "),  paste(has_donor_floor_ipf,  collapse = ", "),
    "",
    "NOTE: that list is a DONOR-COUNT GATE, not a power calculation. Four of the",
    "compartments it admits have AGTR1 identically zero in a whole arm (ATI COPD",
    "and IPF; ATII Control; Endothelial COPD; SMC Control). For actual",
    "sensitivity read agtr1_copd_mde.tsv.")
writeLines(readme, file.path(opt$outdir, "agtr1_copd_README.txt"))

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
