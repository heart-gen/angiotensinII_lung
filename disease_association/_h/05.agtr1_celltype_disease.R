## =============================================================================
## Is AGTR1 itself disease-regulated, and WHERE in the stroma?
##
## The disease figure's first three panels show that graded pericyte
## injury-program engagement tracks disease. This script asks the receptor-level
## question that follows: does AGTR1 expression itself move with disease, and is
## that movement pericyte-specific? The answer that comes out is the honest one --
## the disease-associated AGTR1 variation is carried by FIBROBLAST populations,
## not by pericytes -- so the script is built to make that comparison fair rather
## than to find a pericyte effect.
##
## Design, deliberately mirroring 03.disease_forest.R so the panels are readable
## on one scale:
##   UNIT      : donor x cell type (from 01.disease_association.R's
##               mean_expr/donor_metadata.tsv, >=10 cells/donor/cell type).
##   ENDPOINT  : donor mean AGTR1 (log-normalised), z-standardised WITHIN cell
##               type. Cell types differ ~5x in baseline AGTR1 (pericytes 0.80 vs
##               alveolar fibroblasts 0.17), so a raw-scale coefficient would
##               compare a big number in pericytes against a small one in
##               fibroblasts and call pericytes the stronger signal by
##               construction. Within-cell-type z puts every cell type on
##               "SD of its own donor distribution", which is the scale the
##               injury-program panels already use. Raw-scale estimates are
##               emitted alongside so nothing is hidden.
##   GROUPS    : Healthy / Fibrotic-ILD / Other, from `lung_condition` via the
##               same regex as 03.disease_forest.R. COPD excluded for the same
##               reason (single study).
##   MODEL     : per cell type, LMM  z_AGTR1 ~ disease_group [+ sex] + (1|dataset)
##               (lm when a cell type spans one dataset), contrasts against
##               Healthy, BH-corrected ACROSS CELL TYPES within each contrast.
##               `sex` is included only where it is recorded for every donor in
##               the cell type -- see sex_covar(). The per-row `covars` column
##               says which cell types are sex-adjusted.
##   SENSITIVITY: the same model on frac_AGTR1_pos (detection rate) -- AGTR1 is
##               dropout-prone, so a mean-expression effect that does not appear
##               in detection is worth seeing as such.
##
## LIMITATION carried into the outputs: `donor_metadata.tsv` is built after
## 01.disease_association.R's carcinoma filter, so the "Other" group here is
## COVID-19/chronic-rhinitis/pneumonia donors and does NOT include the carcinoma
## donors that sit in "Other" in 03.disease_forest.R. The Healthy and Fibrotic-ILD
## groups are unaffected.
##
## REBUILT 2026-07-30. Until this date 01.disease_association.R gated on
## `filter(age > 20)`, which drops missing age, and HLCA does not report age for
## 89% of Fibrotic/ILD stroma cells. Every fibrotic arm here was therefore ~6
## donors out of the 24-43 available, and Myofibroblasts (10 fibrotic donors ->
## 1) and Mesothelium (15 -> 1) fell below the >=3-donor gate and were never
## tested at all. That is why the pre-specified fibroblast-lineage family of the
## independent GSE136831 evaluation contained a compartment this script had no
## estimate for. The gate is now "exclude known minors" and both cell types are
## tested. The >=10-cell floor is UNCHANGED -- the floor was never what excluded
## them, and lowering it re-opens the low-cell dropout problem documented in
## AGTR1_DISEASE_DIRECTION.md.
## =============================================================================
suppressPackageStartupMessages({
    library(optparse); library(data.table); library(dplyr)
    library(lme4); library(lmerTest); library(emmeans)
})

TRI_LEVELS <- c("Healthy", "Fibrotic_ILD", "Other")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--donor-celltype", type = "character",
                default = "mean_expr/donor_metadata.tsv", dest = "donor_celltype"),
    make_option("--donor-map", type = "character",
                default = "mean_expr/donor_dataset_map.tsv", dest = "donor_map"),
    make_option("--outdir", type = "character", default = "mean_expr"),
    ## 3 is the smallest group that admits a within-group SD, so it is the point
    ## below which a contrast stops being estimable rather than an arbitrary
    ## power preference. Lowered from 5 on 2026-07-28; the tested set is
    ## UNCHANGED either way (the excluded cell types have 0-1 fibrotic donors,
    ## not 3-4), so this loosens the rule without altering any result.
    make_option("--min-donors", type = "integer", default = 3L, dest = "min_donors",
                help = "donors required in BOTH Healthy and Fibrotic-ILD for a cell type to be tested")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
wt <- function(x, f) fwrite(x, file.path(opt$outdir, f), sep = "\t")

## identical to 03.disease_forest.R / figures/_h/_fig_common.R
map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc)                                              ~ "Healthy",
        lc %in% c("COPD")                                                  ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis",
              lc, ignore.case = TRUE)                                      ~ "Fibrotic_ILD",
        TRUE                                                               ~ "Other")
}

## ---- assemble the donor x cell-type table ----------------------------------
dc  <- fread(opt$donor_celltype)
map <- fread(opt$donor_map)
map <- unique(map[, .(donor_id = patient, dataset, lung_condition)])

dc <- merge(dc, map, by = "donor_id", all.x = TRUE)
miss <- dc[is.na(dataset), uniqueN(donor_id)]
if (miss > 0) warning(miss, " donors had no dataset/lung_condition match and are dropped")
dc <- dc[!is.na(dataset) & !is.na(lung_condition)]

dc[, disease_group := map_disease_group(lung_condition)]
dc <- dc[disease_group %in% TRI_LEVELS & cell_type != "Unknown"]
## HLCA codes unreported sex as the STRING "unknown", which `!is.na()` does not
## catch, so it was silently entering the model as a third sex level. That level
## is not a sex -- it is a marker for the studies that withhold demographics, and
## after the 2026-07-30 age fix it is ~90% of the fibrotic arm against ~15% of the
## healthy arm. Left as a factor level it would absorb much of the very contrast
## this script estimates. Recoded to NA so the covariate logic in fit_ct() can
## make an explicit decision about it; the study structure it was standing in for
## is already carried by (1 | dataset).
dc[sex %in% c("unknown", "Unknown", ""), sex := NA_character_]
dc[, `:=`(disease_group = factor(disease_group, levels = TRI_LEVELS),
          dataset = factor(dataset), sex = factor(sex))]
cat("\n== donor-level sex availability by disease group ==\n")
print(dcast(dc[, .(n = uniqueN(donor_id)), by = .(disease_group, sex_known = !is.na(sex))],
            disease_group ~ sex_known, value.var = "n", fill = 0L))

cat("\n== donors per cell type x disease group ==\n")
inv <- dcast(dc, cell_type ~ disease_group, fun.aggregate = length, value.var = "donor_id")
print(inv); wt(inv, "agtr1_celltype_disease_inventory.tsv")

## Cell types are tested only where the primary contrast is actually estimable.
testable <- inv[Healthy >= opt$min_donors & Fibrotic_ILD >= opt$min_donors, cell_type]
cat("\ntestable cell types (>=", opt$min_donors, " donors in both Healthy and Fibrotic-ILD): ",
    paste(testable, collapse = ", "), "\n", sep = "")
dc <- dc[cell_type %in% testable]

## within-cell-type z (see header: makes cell types comparable)
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
dc[, `:=`(z_AGTR1 = z(AGTR1_mean), z_AGTR1_frac = z(frac_AGTR1_pos)), by = cell_type]

## Sex is adjusted for only where it is OBSERVED FOR EVERY DONOR in the cell type.
## Adjusting on a covariate recorded for a biased subset is worse than not
## adjusting: after the age fix the donors with known sex are overwhelmingly the
## healthy ones, so `sex` would carry part of the disease contrast rather than
## nuisance variation. Dropping donors to keep sex is the other way to satisfy the
## model, and it is the option that produced the 6-fibrotic-donor arms in the first
## place. The chosen covariate set is written to the `covars` column of every
## output row, so which cell types are sex-adjusted is visible, not implicit.
sex_covar <- function(d) {
    if (anyNA(d$sex)) return(character(0))
    if (uniqueN(d$sex) > 1) "sex" else character(0)
}

## ---- per-cell-type three-group model ---------------------------------------
fit_ct <- function(d, resp) {
    d <- d[is.finite(get(resp))]
    d[, `:=`(disease_group = droplevels(disease_group), dataset = droplevels(dataset))]
    if (uniqueN(d$disease_group) < 2) return(NULL)
    n_ds <- nlevels(d$dataset)
    covars <- sex_covar(d)
    rhs <- paste(c("disease_group", covars, if (n_ds >= 2) "(1 | dataset)"), collapse = " + ")
    form <- as.formula(sprintf("%s ~ %s", resp, rhs))
    fit <- try(if (n_ds >= 2) lmerTest::lmer(form, data = d, REML = TRUE)
               else lm(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    emm <- emmeans(fit, ~ disease_group)
    ct  <- as.data.table(as.data.frame(
        contrast(emm, "trt.vs.ctrl", ref = 1, adjust = "none")))
    setnames(ct, "SE", "se")
    ct[, `:=`(response = resp, n_donors = nrow(d),
              n_healthy  = sum(d$disease_group == "Healthy"),
              n_fibrotic = sum(d$disease_group == "Fibrotic_ILD"),
              n_other    = sum(d$disease_group == "Other"),
              ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se,
              covars = if (length(covars)) paste(covars, collapse = "+") else "none",
              model = if (n_ds >= 2) "LMM (1|dataset)" else "lm (single dataset)")]
    ct[]
}

## Omnibus disease effect (2 df) per cell type. The pairwise contrasts alone
## cannot support "AGTR1 varies more with disease in fibroblasts than in
## pericytes": a near-zero pericyte contrast with a wide CI is imprecision, not
## evidence of stability. The joint test of `disease_group` plus a partial
## eta-squared is the statistic that comparison actually needs, so it is computed
## and carried alongside rather than being read off the contrast table.
omnibus_ct <- function(d, resp) {
    d <- d[is.finite(get(resp))]
    d[, `:=`(disease_group = droplevels(disease_group), dataset = droplevels(dataset))]
    if (uniqueN(d$disease_group) < 2) return(NULL)
    n_ds <- nlevels(d$dataset)
    covars <- sex_covar(d)
    rhs <- paste(c("disease_group", covars, if (n_ds >= 2) "(1 | dataset)"), collapse = " + ")
    fit <- try(if (n_ds >= 2)
                   lmerTest::lmer(as.formula(sprintf("%s ~ %s", resp, rhs)), data = d, REML = TRUE)
               else lm(as.formula(sprintf("%s ~ %s", resp, rhs)), data = d), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    jt <- as.data.table(joint_tests(fit))
    jt <- jt[`model term` == "disease_group"]
    if (!nrow(jt)) return(NULL)
    setnames(jt, c("df1", "df2", "F.ratio", "p.value"), c("df1", "df2", "F", "p_omnibus"),
             skip_absent = TRUE)
    ## partial eta^2 from the F statistic and its dfs -- the share of the
    ## endpoint's variance attributable to disease group, on a 0-1 scale that is
    ## comparable across cell types with different n.
    jt[, partial_eta_sq := (F * df1) / (F * df1 + df2)]
    jt[, `:=`(response = resp, n_donors = nrow(d))]
    jt[, .(response, n_donors, df1, df2, F, p_omnibus, partial_eta_sq)]
}

run_omnibus <- function(resp) {
    res <- rbindlist(lapply(testable, function(ctp) {
        r <- omnibus_ct(dc[cell_type == ctp], resp)
        if (is.null(r)) return(NULL)
        r[, cell_type := ctp][]
    }), fill = TRUE)
    res[, p_omnibus_BH := p.adjust(p_omnibus, method = "BH")]
    setcolorder(res, "cell_type")
    res[order(p_omnibus)]
}

run_all <- function(resp) {
    res <- rbindlist(lapply(testable, function(ctp) {
        r <- fit_ct(dc[cell_type == ctp], resp)
        if (is.null(r)) return(NULL)
        r[, cell_type := ctp][]
    }), fill = TRUE)
    ## BH across cell types, separately within each contrast -- the family is
    ## "which cell type shows this contrast", not "all contrasts everywhere".
    res[, p_BH := p.adjust(p.value, method = "BH"), by = contrast]
    setcolorder(res, c("cell_type", "response", "contrast", "estimate", "se",
                       "ci_lo", "ci_hi", "df", "t.ratio", "p.value", "p_BH"))
    res[order(contrast, p.value)]
}

res_z <- run_all("z_AGTR1")
wt(res_z, "agtr1_celltype_disease_effects.tsv")
cat("\n== AGTR1 (within-cell-type SD units): contrasts vs Healthy ==\n"); print(res_z)

omni_z <- run_omnibus("z_AGTR1")
wt(omni_z, "agtr1_celltype_disease_omnibus.tsv")
cat("\n== omnibus disease effect on AGTR1 per cell type (2 df) ==\n"); print(omni_z)

res_raw <- run_all("AGTR1_mean")
wt(res_raw, "agtr1_celltype_disease_effects_rawscale.tsv")
cat("\n== AGTR1 (raw log-normalised scale) ==\n"); print(res_raw)

res_det <- run_all("z_AGTR1_frac")
wt(res_det, "agtr1_celltype_disease_effects_detection.tsv")
cat("\n== SENSITIVITY: AGTR1 detection rate (within-cell-type SD units) ==\n"); print(res_det)

## ---- headline comparison ----------------------------------------------------
## The claim the figure makes, written out as a table so it cannot drift: rank
## cell types by the magnitude of the Fibrotic-Healthy AGTR1 effect and show
## where pericytes land.
head_dt <- merge(res_z[grepl("Fibrotic", contrast)],
                 omni_z[, .(cell_type, F, df1, df2, p_omnibus, p_omnibus_BH,
                            partial_eta_sq)],
                 by = "cell_type", all.x = TRUE)
head_dt <- head_dt[order(-partial_eta_sq)]
## Lineage is now an EXPLICIT map rather than `grepl("fibroblast") ? F : Mural`.
## The old fallback was safe only while the tested set happened to contain nothing
## but fibroblasts, pericytes and vSMC. The age fix admits Mesothelium, which the
## regex would have filed under "Mural" and quietly folded into the mural arm of
## the fibroblast-vs-mural comparison this table exists to support. Myofibroblasts
## do belong with the fibroblasts -- that is also how the independent GSE136831
## evaluation defines its fibroblast-lineage family.
LINEAGE <- c("Adventitial fibroblasts"   = "Fibroblast",
             "Alveolar fibroblasts"      = "Fibroblast",
             "Peribronchial fibroblasts" = "Fibroblast",
             "Subpleural fibroblasts"    = "Fibroblast",
             "Myofibroblasts"            = "Fibroblast",
             "Pericytes"                 = "Mural",
             "Vascular smooth muscle"    = "Mural",
             "Mesothelium"               = "Mesothelial")
head_dt[, `:=`(rank_by_omnibus = seq_len(.N),
               lineage = LINEAGE[cell_type])]
if (head_dt[is.na(lineage), .N])
    stop("unmapped cell type(s) in LINEAGE: ",
         paste(head_dt[is.na(lineage), cell_type], collapse = ", "))
wt(head_dt, "agtr1_celltype_disease_ranking.tsv")
cat("\n== ranking by omnibus disease effect on AGTR1 (partial eta^2) ==\n")
print(head_dt[, .(rank_by_omnibus, cell_type, lineage, partial_eta_sq, p_omnibus,
                  p_omnibus_BH, estimate, ci_lo, ci_hi)])
cat("\nFibroblast vs mural mean partial eta^2:\n")
print(head_dt[, .(mean_partial_eta_sq = mean(partial_eta_sq)), by = lineage])

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
