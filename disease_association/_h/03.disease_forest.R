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
## Design:
##   PRIMARY endpoint  : donor-level mean pericyte injury-program score =
##                       mean of the z-standardised donor means of the
##                       inflammatory, activated/migratory, and fibrillar
##                       fibroblast-like per-cell scores (continuous, no fraction,
##                       no AGTR1, no basement-membrane -- BM is vascular-support).
##   PRIMARY contrast  : Healthy vs Fibrotic/ILD only. COPD/"Other" kept for
##                       description, excluded from the primary estimand.
##   PRIMARY model     : donor-level LMM  endpoint ~ disease + age + sex + (1|dataset)
##                       (study-adjusted pooled effect, Satterthwaite df).
##   HEADLINE FIGURE   : random-effects META-ANALYSIS forest across the studies
##                       that sampled BOTH groups -- per-study effect + DL pooled
##                       diamond + I^2. This is the direct answer to the
##                       "disease is confounded with batch" objection.
##   SENSITIVITY       : (a) each of the 3 injury components separately (was the
##                       composite driven by one program?); (b) min-cells/donor
##                       threshold sweep; (c) SMOKING -- adaptive: adjust the
##                       disease contrast for smoking only if estimable, else fall
##                       back to a within-Healthy smoking effect and a
##                       never-smoker-restricted contrast, with the availability
##                       table reported either way.
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
STAB_COL <- pick_col("vascular_stabilizing")
stopifnot(!anyNA(INJ_COLS))
cat("injury score columns:\n"); print(INJ_COLS)
cat("stability score column:", STAB_COL, "\n")

need <- c("donor_id", "dataset", "lung_condition", "age", "sex", "smoking_status")
have <- intersect(need, names(meta))
first_ok <- function(x) { y <- x[!is.na(x)]; if (length(y)) y[[1]] else x[NA_integer_] }

donor <- meta[, c(
    .(n_cells = .N),
    lapply(.SD, function(v) mean(v, na.rm = TRUE))
), by = donor_id, .SDcols = c(unname(INJ_COLS), if (!is.na(STAB_COL)) STAB_COL)]
setnames(donor, unname(INJ_COLS), paste0("m_", names(INJ_COLS)))
if (!is.na(STAB_COL)) setnames(donor, STAB_COL, "m_vascular_stabilizing")

covar <- meta[, lapply(.SD, first_ok), by = donor_id,
              .SDcols = setdiff(have, "donor_id")]
donor <- merge(donor, covar, by = "donor_id")

donor[, disease_group := relevel(factor(map_disease_group(lung_condition)), "Healthy")]
donor[, sex := factor(sex)]
donor[, age := suppressWarnings(as.numeric(age))]
donor[, dataset := factor(dataset)]

## ---- attrition + confounding tables (report BEFORE filtering) ---------------
attr_all <- dcast(donor, dataset ~ disease_group, fun.aggregate = length, value.var = "donor_id")
wt(attr_all, "attrition_dataset_by_disease_ALL.tsv")
cat(sprintf("\nAll donors with pericytes: %d\n", nrow(donor))); print(table(donor$disease_group))

donor <- donor[n_cells >= MIN_CELLS]
cat(sprintf("\nAfter min-cells >= %d: %d donors\n", MIN_CELLS, nrow(donor)))
print(table(donor$disease_group))
attr_flt <- dcast(donor, dataset ~ disease_group, fun.aggregate = length, value.var = "donor_id")
wt(attr_flt, "attrition_dataset_by_disease_FILTERED.tsv")

## ---- PRIMARY endpoint: z-standardised mean injury-program score -------------
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
prim <- donor[disease_group %in% c("Healthy", "Fibrotic_ILD")]
prim[, disease_group := droplevels(disease_group)]
for (p in names(INJ_COLS)) prim[[paste0("z_", p)]] <- z(prim[[paste0("m_", p)]])
prim[, injury_program_score := rowMeans(as.matrix(prim[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
if ("m_vascular_stabilizing" %in% names(prim)) prim[, z_vascular_stabilizing := z(m_vascular_stabilizing)]
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

## component decomposition (was the composite driven by ONE program?)
comp_rows <- rbindlist(lapply(names(INJ_COLS), function(p) {
    prim2 <- copy(prim); prim2[, tmp := get(paste0("z_", p))]
    summ_disease(fit_lmm("tmp", prim2), paste0("z_", p))
}))
if ("z_vascular_stabilizing" %in% names(prim)) {
    prim2 <- copy(prim); prim2[, tmp := z_vascular_stabilizing]
    comp_rows <- rbind(comp_rows, summ_disease(fit_lmm("tmp", prim2), "z_vascular_stabilizing"))
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
sweep <- rbindlist(lapply(c(5, 10, 15, 20, 30), function(mc) {
    d <- donor[n_cells >= mc & disease_group %in% c("Healthy", "Fibrotic_ILD")]
    d[, disease_group := droplevels(disease_group)]
    if (d[, uniqueN(disease_group)] < 2) return(NULL)
    for (p in names(INJ_COLS)) d[[paste0("z_", p)]] <- z(d[[paste0("m_", p)]])
    d[, injury_program_score := rowMeans(as.matrix(d[, paste0("z_", names(INJ_COLS)), with = FALSE]))]
    r <- summ_disease(fit_lmm("injury_program_score", d), sprintf("min_cells_%d", mc))
    r[, min_cells := mc][]
}), fill = TRUE)
wt(sweep, "mincells_sensitivity.tsv")
cat("\n== min-cells threshold sweep (primary endpoint) ==\n"); print(sweep)

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
