## =============================================================================
## Disease association of pericyte injury-program engagement, done the
## reviewer-defensible way: a PRE-SPECIFIED Healthy-vs-Fibrotic/ILD contrast on
## ONE interpretable donor-level endpoint, with STUDY treated as a modeled factor
## instead of a footnote.
##
## Why this supersedes the pooled ANCOVA in niche_index/01:
##   * The pooled composite's only significant contrast was Healthy-vs-"Other"
##     (COVID/carcinoma grab-bag), not the mechanistically-motivated Fibrotic/ILD
##     contrast (which was borderline, P=0.055). The headline was off-target.
##   * The z-averaged composite (3 program scores + injury fraction + ...) is
##     arbitrary-weighted, and its most interpretable component (the injury
##     FRACTION) is degenerate: after the BM relabel `state_program` has only 3
##     dominant labels, so injury_frac collapses to "activated_migratory only".
##   * Disease and study are NOT fully confounded here: Healthy and Fibrotic/ILD
##     co-occur WITHIN several studies -- notably Banovich_Kropski_2020 (~21 H /
##     18 F) and Kaminski_2020 (~8 H / 23 F). That within-study signal is the
##     asset the pooled model wasted.
##
##     CAVEAT ON THAT SECOND EXAMPLE (P2-25, added 2026-09-08). The counts above
##     are cell-level and describe the cohort BEFORE the >= MIN_CELLS pericyte
##     filter. Kaminski_2020 does not survive it on the Healthy side: its Healthy
##     donors are wiped out entirely while its fibrotic donors largely survive,
##     so the study contributes no within-study contrast and does NOT enter the
##     forest. The forest is built from the studies that did qualify, and the
##     within-study argument stands on those -- but the header's own motivating
##     example is not among them, which a reader deserves to be told rather than
##     left to discover from the attrition table. `forest_eligibility_audit.tsv`
##     now states, per dataset, exactly which arm failed and by how much, and
##     `forest_floor_sensitivity.tsv` rebuilds the forest at lower floors to show
##     what a fourth study would do to the pooled estimate.
##
## Design:
##   PRIMARY endpoint  : donor-level mean pericyte injury-program score =
##                       mean of the z-standardised donor means of the
##                       inflammatory, activated/migratory, and fibrillar
##                       fibroblast-like per-cell scores (continuous, no fraction,
##                       no AGTR1, no basement-membrane -- BM is vascular-support).
##   PRIMARY contrast  : Healthy vs Fibrotic/ILD only. COPD/"Other" kept for
##                       description, excluded from the primary estimand.
##
##   WHAT "HEALTHY" MEANS HERE (P2-26, stated 2026-09-08). The Healthy arm is 42
##   donors, of which 12 (29%) are `Healthy (tumor adjacent)` -- histologically
##   normal lung resected alongside a tumour. This is deliberate and is standard
##   practice for the HLCA: tumour-adjacent normal lung is the largest source of
##   non-diseased human lung tissue in existence, and excluding it would cut the
##   reference arm by nearly a third for a distinction the `disease` field itself
##   does not draw (all 42 are coded `normal`). Note the interaction with the
##   carcinoma exclusion below: that filter operates on `disease`, so it removes
##   carcinoma CASES from "Other" while correctly leaving tumour-adjacent
##   CONTROLS in "Healthy". The two are not in conflict -- they are different
##   tissue.
##
##   The convention was previously stated in no script, legend or summary, which
##   is the actual defect. It is now stated here, flagged per donor in the
##   endpoint tables (`healthy_tumor_adjacent`), and backed by a labelled
##   sensitivity refit on the 30 unambiguous controls
##   (`primary_effect_tumor_adjacent_SENS.tsv`) so the reader can see the
##   convention is not load-bearing. The PRIMARY model is unchanged and keeps all
##   42.
##   PRIMARY model     : donor-level LMM  endpoint ~ disease + age + sex + (1|dataset)
##                       (study-adjusted pooled effect, Satterthwaite df).
##   HEADLINE FIGURE   : random-effects META-ANALYSIS forest across the studies
##                       that sampled BOTH groups -- per-study effect + DL pooled
##                       diamond + I^2. This is the direct answer to the
##                       "disease is confounded with batch" objection.
##   SENSITIVITY       : (a) each of the 3 injury components separately (was the
##                       composite driven by one program?), plus the 3 non-injury
##                       programs -- `vascular_stabilizing` and, since 2026-09-07,
##                       `synthetic_contractile` and `basement_membrane` (P2-33)
##                       -- as negative controls that are reported but never
##                       composited. All SIX panel scores are now tested on
##                       identical footing; (b) min-cells/donor
##                       threshold sweep; (c) SMOKING -- adaptive: adjust the
##                       disease contrast for smoking only if estimable, else fall
##                       back to a within-Healthy smoking effect and a
##                       never-smoker-restricted contrast, with the availability
##                       table reported either way.
##
## ---------------------------------------------------------------------------
## ADDED 2026-07-28 (continuous-injury revision) -- the THREE-GROUP block.
## The manuscript disease figure now leads with graded engagement across
## Healthy / Fibrotic-ILD / Other-disease rather than with the two-group forest,
## so this script also emits, on ONE internally consistent scale:
##   * `threegroup_*`      -- sex + study-adjusted LMM over the three groups,
##                            with both contrasts against Healthy.
##   * `component_effects_3group.tsv`
##                         -- the same three-group model per program score, so
##                            the figure can show WHICH continuous programs
##                            carry the composite.
##   * `leave_one_study_out_3group.tsv`
##                         -- the three-group model refit dropping each dataset
##                            in turn, reporting the Fibrotic-Healthy contrast.
##                            This replaces the per-study forest as the
##                            robustness panel: a forest asks "is the effect
##                            reproduced within studies", LOSO asks "does any
##                            single study create it". The forest is retained
##                            for the supplement.
## COPD is EXCLUDED from the three-group block. Its 12 donors come from a single
## study (Kaminski_2020), so a COPD estimate here is inseparable from that study;
## COPD is instead evaluated in the independent GSE136831 dataset. "Other" is
## therefore COVID/carcinoma/etc., NOT COPD.
## SCALE: the three-group block z-standardises the programs ONCE, over the
## three-group donor set, and every model in the block (composite, components,
## LOSO) reads those fixed columns. Estimates within the block are therefore
## directly comparable; they are NOT on the same scale as the two-group primary
## above (which standardises over Healthy+Fibrotic only) and the two must not be
## quoted interchangeably.
## =============================================================================
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(lme4); library(lmerTest); library(emmeans); library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
META      <- parse_arg("--meta", "../../pericyte_states/_m/pericytes_states_metadata.tsv.gz")
OUTDIR    <- parse_arg("--outdir", "mixed_model_forest")
MIN_CELLS <- as.integer(parse_arg("--min-cells", "10"))   # pericytes/donor for a stable donor mean
MIN_GRP   <- as.integer(parse_arg("--min-per-group", "2"))# donors/group for a study to enter the forest
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
save_gg <- function(fn, p, w, h) for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)
wt <- function(x, f) fwrite(x, file.path(OUTDIR, f), sep = "\t")

## ---- disease grouping (identical regex to niche_index/01) -------------------
map_disease_group <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc)                                              ~ "Healthy",
        lc %in% c("COPD")                                                  ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis",
              lc, ignore.case = TRUE)                                      ~ "Fibrotic_ILD",
        TRUE                                                               ~ "Other")
}

## ---- load per-cell metadata, aggregate to donor -----------------------------
meta <- fread(META); setnames(meta, 1, "barcode")
if ("age_or_mean_of_age_range" %in% names(meta) && !("age" %in% names(meta)))
    setnames(meta, "age_or_mean_of_age_range", "age")

## resolve per-cell injury score columns (defaults, with grep fallback)
pick_col <- function(prog) {
    cand <- paste0(prog, "_score")
    if (cand %in% names(meta)) return(cand)
    hit <- grep(prog, names(meta), value = TRUE, ignore.case = TRUE)
    hit <- grep("score", hit, value = TRUE, ignore.case = TRUE)
    hit <- hit[!grepl("^mean_|^z_", hit)]
    if (length(hit)) hit[1] else NA_character_
}
INJURY   <- c("inflammatory", "activated_migratory", "fibroblast_like")
INJ_COLS <- setNames(vapply(INJURY, pick_col, ""), INJURY)

## NON-INJURY PROGRAMS, TESTED BUT NOT COMPOSITED (P2-33, added 2026-09-07).
## These are reported alongside the injury components and are NOT part of
## `injury_program_score`. Until 2026-09-07 only `vascular_stabilizing` was here,
## which left `synthetic_contractile` -- a named program of the state model --
## with no donor-level disease test anywhere in the repository. That is the P1-4
## hazard: an untested program acquires a "flat" label by association with the
## programs that were tested. It is now tested, and may be described only from
## `component_effects*.tsv`.
##
## `basement_membrane` joined them the same day. It had been excluded by design
## -- a matrix-stabilizing vascular function rather than an injury program, with
## its own disease arm in `basement_membrane/` (GSE136831 COPD) -- but a
## design-level exclusion is indistinguishable, in the output tables, from the
## P1-4 hazard above: a reader cannot tell an untested program from a null one.
## It is therefore tested HERE too, on the same HLCA donors and the same LMM as
## the other five, so the exclusion from `injury_program_score` rests on the
## composite alone and not on a missing test. Its own module's GSE136831 arm is
## a DIFFERENT question (cohort, disease, and direction) and is unaffected.
##
## All six scores are pre-specified single contrasts reported with their own P;
## as elsewhere in this script there is no across-program multiplicity
## adjustment, so read `component_effects*.tsv` as six descriptive arms.
CONTROL   <- c("vascular_stabilizing", "synthetic_contractile", "basement_membrane")
CTRL_COLS <- setNames(vapply(CONTROL, pick_col, ""), CONTROL)
CTRL_COLS <- CTRL_COLS[!is.na(CTRL_COLS) & nzchar(CTRL_COLS)]
stopifnot(!anyNA(INJ_COLS))
cat("injury score columns:\n"); print(INJ_COLS)
cat("non-injury (tested, not composited) score columns:\n"); print(CTRL_COLS)
cat("all six panel scores tested; only the 3 injury scores are composited\n")

## `disease` is carried so the carcinoma exclusion below can be applied here too.
## `study` is carried for the leave-one-STUDY-out arm (P2-8, added 2026-09-07):
## 5 of the 25 studies span more than one dataset (Sun_2020 x4, Regev_2021 x3,
## Meyer_2021 / Thienpont_2018 / Lafyatis_Rojas_2019 x2), so a dataset drop
## cannot remove those cohorts. The random effect stays `(1 | dataset)`; only the
## leave-one-out grouping gains a study-level arm.
need <- c("donor_id", "dataset", "study", "lung_condition", "age", "sex",
          "smoking_status", "disease")
have <- intersect(need, names(meta))
first_ok <- function(x) { y <- x[!is.na(x)]; if (length(y)) y[[1]] else x[NA_integer_] }

donor <- meta[, c(
    .(n_cells = .N),
    lapply(.SD, function(v) mean(v, na.rm = TRUE))
), by = donor_id, .SDcols = c(unname(INJ_COLS), unname(CTRL_COLS))]
setnames(donor, unname(INJ_COLS), paste0("m_", names(INJ_COLS)))
if (length(CTRL_COLS)) setnames(donor, unname(CTRL_COLS), paste0("m_", names(CTRL_COLS)))

covar <- meta[, lapply(.SD, first_ok), by = donor_id,
              .SDcols = setdiff(have, "donor_id")]
donor <- merge(donor, covar, by = "donor_id")

donor[, disease_group := relevel(factor(map_disease_group(lung_condition)), "Healthy")]
## P2-26: mark the tumour-adjacent controls. They REMAIN Healthy (see header);
## the flag exists so the composition of the reference arm is visible in every
## donor-level table instead of being recoverable only from `lung_condition`.
donor[, healthy_tumor_adjacent :=
          grepl("tumor|tumour", as.character(lung_condition), ignore.case = TRUE) &
          disease_group == "Healthy"]
donor[, sex := factor(sex)]
donor[, age := suppressWarnings(as.numeric(age))]
donor[, dataset := factor(dataset)]
if ("study" %in% names(donor)) donor[, study := factor(study)]

## ---- attrition + confounding tables (report BEFORE filtering) ---------------
attr_all <- dcast(donor, dataset ~ disease_group, fun.aggregate = length, value.var = "donor_id")
wt(attr_all, "attrition_dataset_by_disease_ALL.tsv")
cat(sprintf("\nAll donors with pericytes: %d\n", nrow(donor))); print(table(donor$disease_group))

## ---- carcinoma exclusion (added 2026-07-30) ---------------------------------
## 01.disease_association.R and 05.agtr1_celltype_disease.R have always dropped
## carcinoma donors; this script did not, so "Other" here was a
## COVID/carcinoma/pneumonia grab-bag while "Other" in the cell-type panel was
## COVID/pneumonia only. Two panels of the same figure were labelling different
## cohorts "Other". The exclusion is applied with the same patterns as 01 so the
## groups are now defined identically everywhere.
##
## SCOPE OF THE CHANGE, checked rather than assumed: all 11 carcinoma donors fall
## in "Other" (0 in Healthy, 0 in Fibrotic/ILD), so the primary Fibrotic-vs-Healthy
## comparison keeps every donor it had. What does move is the SD of the three-group
## endpoint, which is standardised over Healthy + Fibrotic + Other -- dropping 12
## Other donors rescales it, so the three-group SD-unit estimates shift slightly.
## The two-group forest is standardised over Healthy + Fibrotic only and is
## therefore numerically untouched.
CANCER_RE <- "carcinoma|adenocarcinoma|cancer"
if ("disease" %in% names(donor)) {
    n_canc <- donor[grepl(CANCER_RE, disease, ignore.case = TRUE), .N]
    cat(sprintf("\nExcluding %d carcinoma donors (by `disease`):\n", n_canc))
    print(table(donor[grepl(CANCER_RE, disease, ignore.case = TRUE), disease_group]))
    donor <- donor[!grepl(CANCER_RE, disease, ignore.case = TRUE)]
    print(table(donor$disease_group))
} else {
    warning("no `disease` column -- carcinoma donors NOT excluded")
}

## P1-15. Keep the UNFILTERED donor table. The min-cells sweep at the foot of
## this script re-thresholds donors, and it used to read `donor` -- which by then
## had already been cut at MIN_CELLS, so the sweep could only ever move upward.
## Its `min_cells_5` row was byte-identical to `min_cells_10` (same 66 donors,
## same 0.727854809057206), and that identity was then explained in the S16 legend
## as "no donor has 5-9 pericytes", which is false: 16 donors do in the metadata,
## 13 of them inside this script's Healthy + Fibrotic/ILD analysis set (7 Healthy,
## 6 Fibrotic/ILD; the other 3 are 2 adenocarcinoma and 1 COPD).
donor_all <- copy(donor)
donor <- donor[n_cells >= MIN_CELLS]
cat(sprintf("\nAfter min-cells >= %d: %d donors (%d retained below it for the sweep)\n",
            MIN_CELLS, nrow(donor), nrow(donor_all) - nrow(donor)))
print(table(donor$disease_group))
attr_flt <- dcast(donor, dataset ~ disease_group, fun.aggregate = length, value.var = "donor_id")
wt(attr_flt, "attrition_dataset_by_disease_FILTERED.tsv")

## ---- PRIMARY endpoint: z-standardised mean injury-program score -------------
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
prim <- donor[disease_group %in% c("Healthy", "Fibrotic_ILD")]
prim[, disease_group := droplevels(disease_group)]
for (p in names(INJ_COLS)) prim[[paste0("z_", p)]] <- z(prim[[paste0("m_", p)]])
prim[, injury_program_score := rowMeans(as.matrix(prim[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
for (p in names(CTRL_COLS))
    if (paste0("m_", p) %in% names(prim)) prim[[paste0("z_", p)]] <- z(prim[[paste0("m_", p)]])
wt(prim, "donor_endpoint_table.tsv")

## ===========================================================================
## PRIMARY: study-adjusted LMM  (Healthy vs Fibrotic/ILD)
## ===========================================================================
## covars defaults to sex only: age is missing for ~75% of fibrotic donors, so
## forcing age collapses the case group to n=6 (a missingness artifact, not a
## null effect). Age enters as an explicit sensitivity on the age-complete subset.
fit_lmm <- function(resp, data, covars = c("sex")) {
    keep <- is.finite(data[[resp]])
    for (cv in covars) keep <- keep & !is.na(data[[cv]])
    d <- data[keep]
    d[, disease_group := droplevels(disease_group)]
    n_ds <- nlevels(droplevels(d$dataset))
    rhs  <- paste(c("disease_group", covars, if (n_ds >= 2) "(1 | dataset)"), collapse = " + ")
    form <- as.formula(sprintf("%s ~ %s", resp, rhs))
    fit  <- if (n_ds >= 2) lmerTest::lmer(form, data = d, REML = TRUE) else lm(form, data = d)
    list(fit = fit, d = d, mixed = n_ds >= 2, covars = paste(covars, collapse = "+"))
}
summ_disease <- function(m, label) {
    emm <- emmeans(m$fit, ~ disease_group)
    ph  <- as.data.frame(pairs(emm, adjust = "none"))       # single pre-specified contrast
    co  <- if (m$mixed) coef(summary(m$fit)) else summary(m$fit)$coefficients
    row <- grep("Fibrotic", rownames(co), value = TRUE)[1]
    data.table(response = label, n = nrow(m$d),
               n_healthy  = sum(m$d$disease_group == "Healthy"),
               n_fibrotic = sum(m$d$disease_group == "Fibrotic_ILD"),
               beta_fib_minus_healthy = co[row, "Estimate"],
               se = co[row, "Std. Error"],
               df = if (m$mixed) co[row, "df"] else m$fit$df.residual,
               t = co[row, if (m$mixed) "t value" else "t value"],
               p = co[row, ncol(co)],
               contrast_estimate = ph$estimate[1], contrast_p = ph$p.value[1],
               covars = m$covars,
               model = if (m$mixed) "LMM (1|dataset)" else "lm (single study)")
}

m_prim <- fit_lmm("injury_program_score", prim)
res_prim <- summ_disease(m_prim, "injury_program_score")
emm_prim <- as.data.frame(emmeans(m_prim$fit, ~ disease_group))
wt(as.data.table(emm_prim), "primary_emmeans.tsv")
wt(res_prim, "primary_effect.tsv")
cat("\n== PRIMARY: injury_program_score, Healthy vs Fibrotic/ILD (sex + study-adjusted LMM) ==\n")
print(res_prim)

## age-adjusted SENSITIVITY (age-complete subset only -- age missing for ~75% of
## fibrotic donors, so this is underpowered by design; direction is what matters)
res_prim_age <- summ_disease(fit_lmm("injury_program_score", prim, covars = c("age", "sex")),
                             "injury_program_score")
wt(res_prim_age, "primary_effect_age_adjusted_SENS.tsv")
cat("\n== SENSITIVITY: same endpoint, +age (age-complete subset) ==\n"); print(res_prim_age)

## ---- SENSITIVITY: drop the tumour-adjacent controls (P2-26) ----------------
## The convention -- tumour-adjacent normal lung counts as Healthy -- is a
## deliberate design decision and the primary model above keeps all 42 controls.
## This arm exists to show the decision is not carrying the result: if the
## contrast holds on the 30 unambiguous controls, the objection closes for good.
## Reported as a labelled sensitivity, never as a replacement estimand.
n_ta <- prim[healthy_tumor_adjacent == TRUE, .N]
if (n_ta > 0) {
    prim_nota <- prim[healthy_tumor_adjacent == FALSE]
    prim_nota[, disease_group := droplevels(disease_group)]
    grp_ok <- all(c("Healthy", "Fibrotic_ILD") %in% levels(prim_nota$disease_group)) &&
        min(table(prim_nota$disease_group)) >= 3L
    if (grp_ok) {
        ## NOTE THE SCALE. `injury_program_score` was z-standardised over the FULL
        ## primary set, and those columns are reused here rather than recomputed,
        ## so this estimate is on the SAME SD scale as the primary and the two are
        ## directly comparable. Re-standardising over the reduced set would change
        ## the unit and make the comparison meaningless -- that is the scale trap
        ## this file warns about elsewhere.
        res_nota <- summ_disease(fit_lmm("injury_program_score", prim_nota),
                                 "injury_program_score")
        res_nota[, `:=`(arm = "excl_tumor_adjacent_controls",
                        n_dropped = n_ta,
                        scale_note = "z from the FULL primary set; comparable to primary_effect.tsv")]
        wt(res_nota, "primary_effect_tumor_adjacent_SENS.tsv")
        cat(sprintf("\n== SENSITIVITY: primary contrast excluding %d tumour-adjacent controls ==\n", n_ta))
        print(res_nota)
    } else {
        wt(data.table(arm = "excl_tumor_adjacent_controls", estimable = FALSE,
                      reason = "a group falls below 3 donors after the exclusion"),
           "primary_effect_tumor_adjacent_SENS.tsv")
    }
}
## Composition of the reference arm, so the 29% is a number in a file.
wt(prim[, .(n_donors = .N), by = .(disease_group, healthy_tumor_adjacent)],
   "healthy_arm_composition.tsv")

## component decomposition (was the composite driven by ONE program?)
comp_rows <- rbindlist(lapply(names(INJ_COLS), function(p) {
    prim2 <- copy(prim); prim2[, tmp := get(paste0("z_", p))]
    summ_disease(fit_lmm("tmp", prim2), paste0("z_", p))
}))
for (p in names(CTRL_COLS)) {
    zc <- paste0("z_", p)
    if (!zc %in% names(prim)) next
    prim2 <- copy(prim); prim2[, tmp := get(zc)]
    comp_rows <- rbind(comp_rows, summ_disease(fit_lmm("tmp", prim2), zc))
}
wt(comp_rows, "component_effects.tsv")
cat("\n== component effects (each program separately) ==\n"); print(comp_rows)

## ===========================================================================
## HEADLINE: random-effects meta-analysis forest (within-study effects)
## ===========================================================================
## per-study Fibrotic - Healthy difference in composite SD units (unadjusted;
## within-study => age/sex/study confounding largely internal). Study enters if
## it has >= MIN_GRP donors in BOTH groups after the min-cells filter.
per_study <- prim[, .(nH = sum(disease_group == "Healthy"),
                      nF = sum(disease_group == "Fibrotic_ILD")), by = dataset][
                      nH >= MIN_GRP & nF >= MIN_GRP]
cat(sprintf("\nForest: %d studies sampled BOTH groups (>= %d each)\n", nrow(per_study), MIN_GRP))
print(per_study)

## ---- WHY EACH DATASET DID OR DID NOT ENTER THE FOREST (P2-25) --------------
## The forest is a minority of the primary donors, and until now nothing said so.
## This table names, per dataset, the donors available BEFORE the min-cells
## filter, the donors surviving it, and which arm failed -- so an attrition that
## wipes one arm of a cohort is visible instead of being inferable only by
## comparing two other tables. It is a diagnostic: nothing downstream reads it.
elig_pre <- donor_all[disease_group %in% c("Healthy", "Fibrotic_ILD"),
                      .(nH_pre = sum(disease_group == "Healthy"),
                        nF_pre = sum(disease_group == "Fibrotic_ILD")), by = dataset]
elig_post <- prim[, .(nH_post = sum(disease_group == "Healthy"),
                      nF_post = sum(disease_group == "Fibrotic_ILD")), by = dataset]
elig <- merge(elig_pre, elig_post, by = "dataset", all = TRUE)
for (cc in c("nH_pre", "nF_pre", "nH_post", "nF_post"))
    elig[!is.finite(get(cc)), (cc) := 0L]
elig[, in_forest := nH_post >= MIN_GRP & nF_post >= MIN_GRP]
elig[, reason := fifelse(
    in_forest, "enters the forest",
    fifelse(nH_pre < MIN_GRP & nF_pre < MIN_GRP,
            sprintf("never sampled both arms (%d H / %d F before the filter)", nH_pre, nF_pre),
    fifelse(nH_post < MIN_GRP & nH_pre >= MIN_GRP,
            sprintf("Healthy arm lost to the >=%d-pericyte filter (%d -> %d)", MIN_CELLS, nH_pre, nH_post),
    fifelse(nF_post < MIN_GRP & nF_pre >= MIN_GRP,
            sprintf("Fibrotic arm lost to the >=%d-pericyte filter (%d -> %d)", MIN_CELLS, nF_pre, nF_post),
            sprintf("below %d donors in an arm (%d H / %d F)", MIN_GRP, nH_post, nF_post)))))]
setorder(elig, -in_forest, -nF_post)
wt(elig, "forest_eligibility_audit.tsv")
cat("\n== forest eligibility, per dataset (P2-25) ==\n"); print(elig)
cat(sprintf("Forest covers %d of %d primary donors (%.0f%%) across %d of %d datasets.\n",
            prim[dataset %in% per_study$dataset, .N], nrow(prim),
            100 * prim[dataset %in% per_study$dataset, .N] / nrow(prim),
            nrow(per_study), uniqueN(prim$dataset)))

study_eff <- rbindlist(lapply(per_study$dataset, function(ds) {
    d <- prim[dataset == ds]; d[, disease_group := droplevels(disease_group)]
    f <- lm(injury_program_score ~ disease_group, data = d)
    co <- summary(f)$coefficients
    row <- grep("Fibrotic", rownames(co), value = TRUE)[1]
    data.table(dataset = as.character(ds),
               nH = sum(d$disease_group == "Healthy"),
               nF = sum(d$disease_group == "Fibrotic_ILD"),
               yi = co[row, "Estimate"], sei = co[row, "Std. Error"])
}))
study_eff <- study_eff[is.finite(yi) & is.finite(sei) & sei > 0]

## DerSimonian-Laird random-effects pooling (self-contained; no metafor dep)
dl_pool <- function(yi, sei) {
    vi <- sei^2; w <- 1 / vi
    fe <- sum(w * yi) / sum(w)
    Q  <- sum(w * (yi - fe)^2); k <- length(yi); dfree <- k - 1
    C  <- sum(w) - sum(w^2) / sum(w)
    tau2 <- max(0, (Q - dfree) / C)
    ws <- 1 / (vi + tau2)
    re <- sum(ws * yi) / sum(ws); se_re <- sqrt(1 / sum(ws))
    I2 <- if (Q > dfree && Q > 0) max(0, (Q - dfree) / Q) * 100 else 0
    list(estimate = re, se = se_re, ci_lo = re - 1.96 * se_re, ci_hi = re + 1.96 * se_re,
         tau2 = tau2, Q = Q, df = dfree, p_Q = pchisq(Q, dfree, lower.tail = FALSE),
         I2 = I2, weights = ws / sum(ws) * 100)
}
pool <- dl_pool(study_eff$yi, study_eff$sei)
study_eff[, ci_lo := yi - 1.96 * sei][, ci_hi := yi + 1.96 * sei]
study_eff[, weight_pct := pool$weights]
wt(study_eff, "forest_per_study.tsv")
pool_row <- data.table(estimate = pool$estimate, se = pool$se, ci_lo = pool$ci_lo,
                       ci_hi = pool$ci_hi, tau2 = pool$tau2, Q = pool$Q, df = pool$df,
                       p_Q = pool$p_Q, I2 = pool$I2,
                       p_pooled = 2 * pnorm(-abs(pool$estimate / pool$se)))
wt(pool_row, "forest_pooled_RE.tsv")
cat("\n== random-effects pooled Fibrotic - Healthy (SD units) ==\n"); print(pool_row)

## ---- FOREST-FLOOR SENSITIVITY (P2-25) --------------------------------------
## The forest's membership is set by MIN_CELLS, a threshold chosen for the
## stability of a donor MEAN -- not for meta-analytic coverage. A cohort can be
## excluded from the forest for having few pericytes per donor even when it holds
## the largest fibrotic sample in the study. So: rebuild the whole forest at
## lower floors and report what enters and what the pooled estimate does.
##
## This is a SENSITIVITY, not a proposal to lower the primary floor. A lower
## floor buys studies at the cost of noisier donor means, and the trade is only
## worth making if the pooled estimate is stable across it -- which is exactly
## what this table lets a reader check. The primary forest remains MIN_CELLS.
##
## SCALE: `z` is recomputed within each rung's own Healthy+Fibrotic set, the same
## convention `mincells_sensitivity.tsv` uses, so each row is in that rung's SD
## units. Rows are comparable in DIRECTION and in study membership; treat the
## magnitudes as within-rung.
FOREST_RUNGS <- sort(unique(c(MIN_CELLS, 5L, 3L, 2L)), decreasing = TRUE)
fs_ps <- data.table()
forest_sens <- rbindlist(lapply(FOREST_RUNGS, function(mc) {
    dd <- donor_all[n_cells >= mc & disease_group %in% c("Healthy", "Fibrotic_ILD")]
    dd[, disease_group := droplevels(disease_group)]
    if (uniqueN(dd$disease_group) < 2L) return(NULL)
    for (pp in names(INJ_COLS)) dd[[paste0("z_", pp)]] <- z(dd[[paste0("m_", pp)]])
    dd[, injury_program_score :=
           rowMeans(as.matrix(dd[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
    ps <- dd[, .(nH = sum(disease_group == "Healthy"),
                 nF = sum(disease_group == "Fibrotic_ILD")), by = dataset][
                 nH >= MIN_GRP & nF >= MIN_GRP]
    if (!nrow(ps)) return(data.table(min_cells = mc, n_studies = 0L, n_donors_forest = 0L,
                                     n_donors_primary = nrow(dd), estimate = NA_real_,
                                     se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                                     I2 = NA_real_, p_pooled = NA_real_,
                                     studies = "", new_vs_primary = ""))
    se_dt <- rbindlist(lapply(ps$dataset, function(ds) {
        d <- dd[dataset == ds]; d[, disease_group := droplevels(disease_group)]
        f <- lm(injury_program_score ~ disease_group, data = d)
        co <- summary(f)$coefficients
        r <- grep("Fibrotic", rownames(co), value = TRUE)[1]
        data.table(dataset = as.character(ds), yi = co[r, "Estimate"], sei = co[r, "Std. Error"])
    }))
    se_dt <- se_dt[is.finite(yi) & is.finite(sei) & sei > 0]
    if (nrow(se_dt) < 2L) return(NULL)
    ## Per-study rows at every rung, with the MEDIAN PERICYTE COUNT per arm.
    ## Without that column a rung looks like it simply added a study, when what
    ## it actually added may be donor means estimated from a handful of cells --
    ## which is the thing MIN_CELLS exists to prevent, so it has to be visible.
    fs_ps <<- rbind(fs_ps, merge(se_dt, dd[, .(
        nH = sum(disease_group == "Healthy"), nF = sum(disease_group == "Fibrotic_ILD"),
        med_cells_H = as.numeric(median(n_cells[disease_group == "Healthy"])),
        med_cells_F = as.numeric(median(n_cells[disease_group == "Fibrotic_ILD"]))),
        by = dataset][, dataset := as.character(dataset)],
        by = "dataset")[, min_cells := mc][], fill = TRUE)
    pl <- dl_pool(se_dt$yi, se_dt$sei)
    data.table(min_cells = mc, n_studies = nrow(se_dt),
               n_donors_forest = dd[dataset %in% se_dt$dataset, .N],
               n_donors_primary = nrow(dd),
               estimate = pl$estimate, se = pl$se, ci_lo = pl$ci_lo, ci_hi = pl$ci_hi,
               I2 = pl$I2, p_pooled = 2 * pnorm(-abs(pl$estimate / pl$se)),
               studies = paste(sort(se_dt$dataset), collapse = "|"),
               new_vs_primary = paste(sort(setdiff(se_dt$dataset,
                                                   study_eff$dataset)), collapse = "|"))
}), fill = TRUE)
if (nrow(fs_ps)) {
    setcolorder(fs_ps, c("min_cells", "dataset", "nH", "nF",
                         "med_cells_H", "med_cells_F", "yi", "sei"))
    setorder(fs_ps, -min_cells, dataset)
    wt(fs_ps, "forest_floor_per_study.tsv")
    cat("\n== per-study effects at each forest floor (P2-25) ==\n"); print(fs_ps)
}
if (nrow(forest_sens)) {
    wt(forest_sens, "forest_floor_sensitivity.tsv")
    cat("\n== forest rebuilt at lower pericyte floors (P2-25) ==\n")
    print(forest_sens[, .(min_cells, n_studies, n_donors_forest, estimate,
                          ci_lo, ci_hi, I2, p_pooled, new_vs_primary)])
}

## ---- forest plot (manuscript style: no title, direct labels) ----------------
study_lab <- sprintf("%s  (%d H / %d F)", study_eff$dataset, study_eff$nH, study_eff$nF)
pool_lab  <- sprintf("RE pooled  (I2=%.0f%%)", pool$I2)
fp <- rbind(
    data.table(label = study_lab, y = study_eff$yi, lo = study_eff$ci_lo,
               hi = study_eff$ci_hi, w = study_eff$weight_pct, kind = "study"),
    data.table(label = pool_lab, y = pool$estimate, lo = pool$ci_lo, hi = pool$ci_hi,
               w = max(study_eff$weight_pct), kind = "pooled"))
## studies ordered by effect size, pooled diamond pinned to the bottom row
ord <- c(study_lab[order(study_eff$yi)], pool_lab)
fp[, label := factor(label, levels = rev(ord))]
pal <- c(study = "#0072B2", pooled = "#D55E00")
pf <- ggplot(fp, aes(y = label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = kind), height = 0.22, linewidth = 0.6) +
    geom_point(aes(x = y, colour = kind, size = w, shape = kind)) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c(study = 16, pooled = 18), guide = "none") +
    scale_size_continuous(range = c(2.5, 6.5), guide = "none") +
    labs(x = "Fibrotic/ILD minus Healthy: injury-program score (SD units)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
save_gg(file.path(OUTDIR, "forest_injury_program"), pf, 7.2, 0.6 + 0.5 * nrow(fp))

## ===========================================================================
## THREE-GROUP BLOCK: graded injury engagement across Healthy / Fibrotic-ILD /
## Other-disease  (main disease figure, panels A-C)
## ===========================================================================
## COPD dropped: single-study (Kaminski_2020), so not separable from study here.
TRI_LEVELS <- c("Healthy", "Fibrotic_ILD", "Other")
tri <- donor[disease_group %in% TRI_LEVELS]
tri[, disease_group := factor(as.character(disease_group), levels = TRI_LEVELS)]
## z ONCE over the three-group set; every model below reads these fixed columns.
for (p in names(INJ_COLS)) tri[[paste0("z_", p)]] <- z(tri[[paste0("m_", p)]])
for (p in names(CTRL_COLS))
    if (paste0("m_", p) %in% names(tri)) tri[[paste0("z_", p)]] <- z(tri[[paste0("m_", p)]])
tri[, injury_program_score := rowMeans(as.matrix(tri[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
wt(tri, "donor_endpoint_table_3group.tsv")
cat("\n== THREE-GROUP set (COPD excluded) ==\n"); print(table(tri$disease_group))

## Fit once, return BOTH contrasts vs Healthy plus the adjusted marginal means.
## Contrasts come from emmeans (not the coefficient table) so the reference level
## is explicit and the same call serves the composite and each component.
fit_tri <- function(resp, data, covars = c("sex")) {
    keep <- is.finite(data[[resp]])
    for (cv in covars) keep <- keep & !is.na(data[[cv]])
    d <- data[keep]
    d[, disease_group := droplevels(disease_group)]
    n_ds <- nlevels(droplevels(d$dataset))
    rhs  <- paste(c("disease_group", covars, if (n_ds >= 2) "(1 | dataset)"), collapse = " + ")
    form <- as.formula(sprintf("%s ~ %s", resp, rhs))
    fit  <- if (n_ds >= 2) lmerTest::lmer(form, data = d, REML = TRUE) else lm(form, data = d)
    list(fit = fit, d = d, mixed = n_ds >= 2, covars = paste(covars, collapse = "+"))
}
summ_tri <- function(m, label) {
    emm <- emmeans(m$fit, ~ disease_group)
    ## "trt.vs.ctrl" against Healthy (level 1); unadjusted -- two pre-specified
    ## comparisons, each reported with its own P.
    ct <- as.data.frame(contrast(emm, "trt.vs.ctrl", ref = 1, adjust = "none"))
    ct <- as.data.table(ct)
    setnames(ct, c("estimate", "SE"), c("estimate", "se"), skip_absent = TRUE)
    ct[, `:=`(response = label, n = nrow(m$d),
              ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se,
              covars = m$covars,
              model = if (m$mixed) "LMM (1|dataset)" else "lm (single study)")]
    ns <- m$d[, .N, by = disease_group]
    list(contrasts = ct[], emmeans = as.data.table(as.data.frame(emm))[
            , `:=`(response = label)][], n_by_group = ns)
}

m_tri  <- fit_tri("injury_program_score", tri)
r_tri  <- summ_tri(m_tri, "injury_program_score")
wt(r_tri$contrasts, "threegroup_effects.tsv")
wt(r_tri$emmeans,   "threegroup_emmeans.tsv")
wt(r_tri$n_by_group, "threegroup_n_by_group.tsv")
cat("\n== THREE-GROUP composite: contrasts vs Healthy (sex + study-adjusted LMM) ==\n")
print(r_tri$contrasts)
cat("\n-- adjusted marginal means --\n"); print(r_tri$emmeans)

## component decomposition on the SAME three-group model
comp_tri <- rbindlist(lapply(c(names(INJ_COLS),
                               names(CTRL_COLS)[paste0("z_", names(CTRL_COLS)) %in% names(tri)]),
                             function(p) {
    t2 <- copy(tri); t2[, tmp := get(paste0("z_", p))]
    summ_tri(fit_tri("tmp", t2), paste0("z_", p))$contrasts
}))
wt(comp_tri, "component_effects_3group.tsv")
cat("\n== THREE-GROUP component effects (each program separately) ==\n"); print(comp_tri)

## ---- leave-one-study-out on the three-group model --------------------------
## Refit dropping each dataset in turn and re-report the Fibrotic-Healthy
## contrast. Datasets are dropped one at a time from the FULL three-group set
## (including datasets that contribute only one group -- removing those still
## changes the study random effect and the Other arm, so they are legitimate
## refits). The z columns are held FIXED at their full-set values so every refit
## is on one scale.
## GROUPING FIXED 2026-09-07 (P2-8). This looped over `dataset` while writing
## `leave_one_study_out_3group.tsv`. Five studies span more than one dataset, so
## dropping a dataset leaves the rest of its study in the fit and cannot answer
## "is this carried by one cohort?". Both arms are now emitted under names that
## say what they drop; the study-level arm is the primary.
## factor() rather than droplevels(): `dataset` is coerced to a factor upstream
## but `study` is not, and droplevels() has no character method.
loso_by <- function(by) rbindlist(lapply(levels(droplevels(factor(tri[[by]]))), function(g) {
    d <- tri[get(by) != g]
    d[, disease_group := droplevels(disease_group)]
    if (!all(c("Healthy", "Fibrotic_ILD") %in% levels(d$disease_group))) return(NULL)
    if (d[disease_group == "Healthy", .N] < 3 || d[disease_group == "Fibrotic_ILD", .N] < 3)
        return(NULL)
    r <- summ_tri(fit_tri("injury_program_score", d), "injury_program_score")$contrasts
    r <- r[grepl("Fibrotic", contrast)]
    ## `dropped_level` holds the VALUE that was dropped, and `loso_by` names the
    ## grouping it came from. The first version of this had `dropped_level = by`,
    ## i.e. the constant string "study" -- a column whose name promised the level
    ## and delivered the level's TYPE. It broke the S16A figure panel, which
    ## reasonably read `dropped_level` for a label and got 17 identical values.
    ## Keep both: `dropped_level` is the arm-independent column a consumer should
    ## read, and `dropped_<by>` stays for back-compatibility.
    r[, `:=`(dropped = g, dropped_level = g, loso_by = by,
             n_dropped = tri[get(by) == g, .N],
             n_studies_left = uniqueN(d$dataset), n_left = nrow(d))]
    setnames(r, "dropped", paste0("dropped_", by))
    r[]
}), fill = TRUE)

loso        <- loso_by("study")
loso_datset <- loso_by("dataset")
setcolorder(loso, c("dropped_study", "n_dropped", "contrast", "estimate", "se",
                    "ci_lo", "ci_hi"))
setcolorder(loso_datset, c("dropped_dataset", "n_dropped", "contrast", "estimate", "se",
                           "ci_lo", "ci_hi"))
wt(loso, "leave_one_study_out_3group.tsv")
wt(loso_datset, "leave_one_dataset_out_3group.tsv")
cat(sprintf("\n== LOSO (three-group, Fibrotic - Healthy): %d STUDY-level refits, %d with P < 0.05 ==\n",
            nrow(loso), sum(loso$p.value < 0.05)))
cat(sprintf("   (dataset-level arm: %d refits, %d with P < 0.05 -> leave_one_dataset_out_3group.tsv)\n",
            nrow(loso_datset), sum(loso_datset$p.value < 0.05)))
cat(sprintf("   estimate range %.3f to %.3f (full-data %.3f)\n",
            min(loso$estimate), max(loso$estimate),
            r_tri$contrasts[grepl("Fibrotic", contrast), estimate]))
print(loso[order(estimate)])

## ===========================================================================
## SMOKING sensitivity (adaptive) -- major reviewer comment
## ===========================================================================
norm_smoke <- function(x) {
    s <- tolower(trimws(as.character(x)))
    dplyr::case_when(
        is.na(s) | s %in% c("", "nan", "na", "unknown", "not available", "none") ~ NA_character_,
        grepl("never|non.?smoker|no", s)                                          ~ "never",
        grepl("former|ex.?smoker|past|previous", s)                               ~ "ever",
        grepl("current|active|yes|smoker|cigar", s)                               ~ "ever",
        TRUE                                                                       ~ NA_character_)
}
smoke_note <- c()
if ("smoking_status" %in% names(prim)) {
    prim[, smoke := norm_smoke(smoking_status)]
    avail <- prim[, .(n = .N, n_known = sum(!is.na(smoke)),
                      n_never = sum(smoke == "never", na.rm = TRUE),
                      n_ever  = sum(smoke == "ever",  na.rm = TRUE)), by = disease_group]
    wt(avail, "smoking_availability_by_disease.tsv")
    cat("\n== smoking availability by disease group ==\n"); print(avail)

    kn <- prim[!is.na(smoke)]
    grp_ok <- kn[, .(nlev = uniqueN(smoke), n = .N), by = disease_group]
    estimable <- nrow(grp_ok) == 2 && all(grp_ok$n >= 3) &&
                 kn[, uniqueN(smoke)] >= 2 &&
                 all(grp_ok$disease_group %in% c("Healthy", "Fibrotic_ILD"))

    ## (a) disease contrast ADJUSTED for smoking -- only if estimable
    if (estimable) {
        kn2 <- copy(kn); kn2[, smoke := factor(smoke)]
        nds <- nlevels(droplevels(kn2$dataset))
        form <- if (nds >= 2) injury_program_score ~ disease_group + age + sex + smoke + (1 | dataset)
                else            injury_program_score ~ disease_group + age + sex + smoke
        fit <- if (nds >= 2) lmerTest::lmer(form, data = kn2) else lm(form, data = kn2)
        co <- if (nds >= 2) coef(summary(fit)) else summary(fit)$coefficients
        row <- grep("Fibrotic", rownames(co), value = TRUE)[1]
        smk_adj <- data.table(analysis = "disease adjusted for smoking",
                              n = nrow(kn2),
                              beta_fib_minus_healthy = co[row, "Estimate"],
                              se = co[row, "Std. Error"], p = co[row, ncol(co)])
        wt(smk_adj, "smoking_adjusted_disease.tsv")
        cat("\n== disease contrast adjusted for smoking ==\n"); print(smk_adj)
        smoke_note <- c(smoke_note, "smoking-adjusted disease contrast WAS estimable (see smoking_adjusted_disease.tsv)")
    } else {
        smoke_note <- c(smoke_note,
            "smoking-adjusted disease contrast NOT estimable (smoking missing/constant within a disease group)")
        cat("\nsmoking-adjusted disease contrast NOT estimable given availability above.\n")
    }

    ## (b) within-Healthy smoking effect (is the endpoint sensitive to smoking at all?)
    h <- kn[disease_group == "Healthy"]
    if (h[, uniqueN(smoke)] >= 2 && nrow(h) >= 6) {
        nds <- nlevels(droplevels(h$dataset))
        form <- if (nds >= 2) injury_program_score ~ smoke + age + sex + (1 | dataset)
                else            injury_program_score ~ smoke + age + sex
        fit <- if (nds >= 2) lmerTest::lmer(form, data = h) else lm(form, data = h)
        co <- if (nds >= 2) coef(summary(fit)) else summary(fit)$coefficients
        row <- grep("smoke", rownames(co), value = TRUE)[1]
        smk_h <- data.table(analysis = "within-Healthy smoking effect", n = nrow(h),
                            beta_ever_minus_never = co[row, "Estimate"],
                            se = co[row, "Std. Error"], p = co[row, ncol(co)])
        wt(smk_h, "smoking_within_healthy.tsv")
        cat("\n== within-Healthy smoking effect on the endpoint ==\n"); print(smk_h)
    }

    ## (c) never-smoker-restricted disease contrast
    ns <- kn[smoke == "never"]; ns[, disease_group := droplevels(disease_group)]
    if (all(c("Healthy", "Fibrotic_ILD") %in% levels(ns$disease_group)) &&
        ns[disease_group == "Healthy", .N] >= 3 && ns[disease_group == "Fibrotic_ILD", .N] >= 3) {
        m_ns <- fit_lmm("injury_program_score", ns)
        smk_ns <- summ_disease(m_ns, "never-smokers only")
        wt(smk_ns, "smoking_neversmoker_restricted.tsv")
        cat("\n== never-smoker-restricted disease contrast ==\n"); print(smk_ns)
        smoke_note <- c(smoke_note, "never-smoker-restricted contrast reported")
    } else {
        smoke_note <- c(smoke_note, "never-smoker-restricted contrast not estimable (too few never-smokers in a group)")
    }
} else {
    smoke_note <- "no smoking_status column found"
}
writeLines(smoke_note, file.path(OUTDIR, "smoking_sensitivity_NOTE.txt"))

## ===========================================================================
## min-cells threshold sensitivity (attrition robustness)
## ===========================================================================
## Reads `donor_all`, NOT `donor`: see the P1-15 note where donor_all is taken.
## The rungs below MIN_CELLS are the whole point of the sweep -- they are the ones
## that ADD donors, and therefore the ones that stress the result rather than
## re-confirming it on a subset.
SWEEP_RUNGS <- c(5, 10, 15, 20, 30)
dsw <- donor_all[disease_group %in% c("Healthy", "Fibrotic_ILD")]
cat("\n== donors available at each sweep rung (Healthy + Fibrotic/ILD) ==\n")
print(data.table(
    min_cells  = SWEEP_RUNGS,
    n_donors   = sapply(SWEEP_RUNGS, function(mc) dsw[n_cells >= mc, .N]),
    n_healthy  = sapply(SWEEP_RUNGS, function(mc) dsw[n_cells >= mc & disease_group == "Healthy", .N]),
    n_fibrotic = sapply(SWEEP_RUNGS, function(mc) dsw[n_cells >= mc & disease_group == "Fibrotic_ILD", .N])))
sweep <- rbindlist(lapply(SWEEP_RUNGS, function(mc) {
    d <- donor_all[n_cells >= mc & disease_group %in% c("Healthy", "Fibrotic_ILD")]
    d[, disease_group := droplevels(disease_group)]
    if (d[, uniqueN(disease_group)] < 2) return(NULL)
    for (p in names(INJ_COLS)) d[[paste0("z_", p)]] <- z(d[[paste0("m_", p)]])
    d[, injury_program_score := rowMeans(as.matrix(d[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
    r <- summ_disease(fit_lmm("injury_program_score", d), sprintf("min_cells_%d", mc))
    r[, min_cells := mc][]
}), fill = TRUE)
## The guard that makes P1-15 impossible to repeat silently. A rung below
## MIN_CELLS must admit strictly more donors than MIN_CELLS does, unless the
## cohort genuinely has no donor in that band -- in which case say so, rather than
## letting a duplicated row be read as robustness.
low <- SWEEP_RUNGS[SWEEP_RUNGS < MIN_CELLS]
if (length(low) > 0) {
    band <- donor_all[disease_group %in% c("Healthy", "Fibrotic_ILD") &
                      n_cells >= min(low) & n_cells < MIN_CELLS, .N]
    n_lo <- sweep[min_cells == min(low), n]
    n_at <- sweep[min_cells == MIN_CELLS, n]
    if (band > 0 && length(n_lo) && length(n_at) && n_lo <= n_at)
        stop("min-cells sweep: rung ", min(low), " admits no more donors than rung ",
             MIN_CELLS, " (", n_lo, " vs ", n_at, ") although ", band, " donors sit ",
             "in the band. The sweep is reading an already-filtered table -- this ",
             "is defect P1-15.")
    cat(sprintf("\n[min-cells sweep] %d donors sit in the %d-%d band and are ADMITTED by the low rung (n %d -> %d).\n",
                band, min(low), MIN_CELLS - 1L, n_at, n_lo))
}
wt(sweep, "mincells_sensitivity.tsv")
cat("\n== min-cells threshold sweep (primary endpoint) ==\n"); print(sweep)

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
