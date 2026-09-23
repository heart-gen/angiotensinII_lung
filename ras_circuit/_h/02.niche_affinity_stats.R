## Figure 5B -- which component of the alveolar-capillary niche tracks pericyte AGTR1?
##
## Unit: the DONOR (pericytes averaged within donor; >= --min-cells pericytes).
## Exposure: donor pericyte AGTR1 pseudobulk (expm1 -> mean -> log1p over the same
## pericytes whose affinities are averaged), z within dataset. `_detect` arm uses
## the detection fraction instead.
## Outcome: the four niche-affinity axes (01.niche_affinity.py), z within dataset.
##
## PRIMARY MODEL (long format, one row per donor x compartment):
##   affinity_z ~ compartment * AGTR1_z + sex + mean_log10_total_counts
##               + (1 | study) + (1 | donor_id)
## Global test: likelihood-ratio test of the compartment x AGTR1 interaction (3 df,
## ML fits). Per-compartment slopes: emtrends, BH over the four. Two pre-specified
## contrasts on the slopes: aerocyte - general capillary, and epithelial
## (AT1, AT2) - endothelial (aerocyte, gCap).
##
## NO `+ age` IN THE PRIMARY. In this atlas age missingness is a STUDY property,
## so `+ age` deletes whole studies (memory: age-covariate-study-confounding). It
## is fitted as the labelled `_ageadj` arm, which records the donors and studies
## it drops. `(1 | study)` is load-bearing and is always present.
##
## DETECTION-MATCHED NULL. A single-gene exposure regressed on a continuous
## outcome has a non-zero null (memory: gene-set-score-null-not-zero). Every
## primary statistic is refitted with each detection-matched gene in place of
## AGTR1, and p_emp is referred to the null's own centre.
##
## Also written: the legacy aggregate result beside its study-guarded refit, so the
## change in model is visible rather than silent.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(lmerTest)
    library(emmeans)
})
source("../_h/_stats_common.R")
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--affinity", default = "./pericyte_niche_affinity.tsv.gz"),
    make_option("--gene-pb", dest = "gene_pb",
                default = "./niche_affinity_gene_pseudobulk.tsv.gz"),
    make_option("--legacy", default = paste0("../../localization/airspace_analysis/",
                                             "_m/airspace/airspace_effect_AGTR1.csv")),
    make_option("--outdir", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 10L, dest = "min_cells"),
    make_option("--max-null", type = "integer", default = 0L, dest = "max_null")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

FOUR <- c("AT1", "AT2", "EC_aerocyte", "EC_gcap")
REF  <- c("SMC", "EC_arterial")

cells <- read_req(opt$affinity)
setnames(cells, 1, "index")
gpb   <- read_req(opt$gene_pb)
message(sprintf("pericytes: %d | donors: %d", nrow(cells), uniqueN(cells$donor_id)))

## ---- donor table -------------------------------------------------------------
aff_cols <- paste0("affinity_", c(FOUR, REF))
stopifnot(all(aff_cols %in% names(cells)))
donor <- cells[, c(lapply(.SD, mean),
                   .(n_cells = .N,
                     mean_log10_total_counts = mean(log10_total_counts),
                     frac_AGTR1_pos = mean(AGTR1_detect > 0),
                     mean_airspace_score = mean(airspace_score),
                     study = study[1], dataset = dataset[1],
                     sex = as.character(sex[1]),
                     lung_condition = as.character(lung_condition[1]),
                     age = suppressWarnings(as.numeric(as.character(
                         age_or_mean_of_age_range[1]))))),
               by = donor_id, .SDcols = aff_cols]
## Unknown sex is a level, not a deletion: dropping it could filter by study the
## same way `+ age` did.
donor[is.na(sex) | !nzchar(sex) | sex %in% c("unknown", "nan", "NA"), sex := "unknown"]
## Disease group is a NUISANCE covariate in the `_disease_adj` arm and a filter in
## `_healthy`; this repository makes no disease claims. Same regex as
## figures/_h/_fig_common.R::map_disease so the groups cannot disagree.
donor[, disease_group := fcase(
    grepl("^Healthy", lung_condition), "Healthy",
    lung_condition == "COPD", "COPD",
    grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lung_condition,
          ignore.case = TRUE), "Fibrotic_ILD",
    default = "Other")]

agtr1 <- gpb[gene == "AGTR1", .(donor_id, AGTR1_pb = expr, AGTR1_pb_detect = detect)]
donor <- merge(donor, agtr1, by = "donor_id")
fwrite(donor, file.path(opt$outdir, "niche_affinity_donor.tsv"), sep = "\t")

prep <- function(d) {
    d <- copy(d)
    for (v in c(aff_cols, "mean_airspace_score"))
        d[[paste0(v, "_z")]] <- z_within_dataset(d[[v]], d$dataset)
    d[, AGTR1_z := z_within_dataset(AGTR1_pb, dataset)]
    d[, AGTR1_detect_z := z_within_dataset(AGTR1_pb_detect, dataset)]
    d[, sex := factor(sex)]
    d
}

to_long <- function(d, xcol) {
    L <- melt(d, id.vars = c("donor_id", "study", "dataset", "sex", "disease_group",
                             "mean_log10_total_counts", xcol),
              measure.vars = paste0("affinity_", FOUR, "_z"),
              variable.name = "compartment", value.name = "affinity_z")
    L[, compartment := factor(sub("^affinity_(.*)_z$", "\\1", compartment),
                              levels = FOUR)]
    setnames(L, xcol, "X")
    L
}

## ---- the primary fit, reused for every arm and every null gene -------------------
fit_long <- function(d, xcol, extra = NULL, quiet = TRUE) {
    L <- to_long(d, xcol)
    extra_new <- setdiff(extra, names(L))
    if (length(extra_new)) L <- merge(L, d[, c("donor_id", extra_new), with = FALSE],
                                      by = "donor_id")
    cov <- paste(c("sex", "mean_log10_total_counts", extra), collapse = " + ")
    ## A sex factor with a single level cannot enter the model.
    if (nlevels(droplevels(L$sex)) < 2) cov <- sub("sex \\+ ", "", cov)
    f1 <- as.formula(paste("affinity_z ~ compartment * X +", cov,
                           "+ (1 | study) + (1 | donor_id)"))
    f0 <- as.formula(paste("affinity_z ~ compartment + X +", cov,
                           "+ (1 | study) + (1 | donor_id)"))
    m1 <- try(suppressMessages(lmer(f1, data = L)), silent = quiet)
    if (inherits(m1, "try-error")) return(NULL)
    m1ml <- suppressMessages(update(m1, REML = FALSE))
    m0ml <- suppressMessages(lmer(f0, data = L, REML = FALSE))
    lrt <- anova(m0ml, m1ml)
    tr <- emtrends(m1, ~ compartment, var = "X")
    sl <- as.data.table(summary(tr, infer = c(TRUE, TRUE)))
    lv <- levels(L$compartment); stopifnot(identical(lv, FOUR))
    ct <- as.data.table(summary(contrast(tr, method = list(
        "EC_aerocyte - EC_gcap" = c(0, 0, 1, -1),
        "epithelial - endothelial" = c(0.5, 0.5, -0.5, -0.5))), infer = c(TRUE, TRUE)))
    list(lrt_chisq = lrt$Chisq[2], lrt_df = lrt$Df[2], lrt_p = lrt$`Pr(>Chisq)`[2],
         slopes = sl, contrasts = ct, singular = isSingular(m1),
         n_donors = uniqueN(L$donor_id), n_studies = uniqueN(L$study))
}

## ---- arms ------------------------------------------------------------------------
base <- prep(donor[n_cells >= opt$min_cells])
d20  <- prep(donor[n_cells >= 20])
age_ok <- base[!is.na(age)]
arms <- list(
    list(arm = "primary", d = base, x = "AGTR1_z", extra = NULL),
    list(arm = "_detect", d = base, x = "AGTR1_detect_z", extra = NULL),
    list(arm = "_mincells20", d = d20, x = "AGTR1_z", extra = NULL),
    list(arm = "_ageadj", d = prep(donor[n_cells >= opt$min_cells & !is.na(age)]),
         x = "AGTR1_z", extra = "age"),
    list(arm = "_disease_adj", d = base, x = "AGTR1_z", extra = "disease_group"),
    list(arm = "_healthy", d = prep(donor[n_cells >= opt$min_cells &
                                          disease_group == "Healthy"]),
         x = "AGTR1_z", extra = NULL))
studies_all <- sort(unique(base$study))

glob <- list(); slopes <- list(); contr <- list()
for (A in arms) {
    r <- fit_long(A$d, A$x, A$extra, quiet = FALSE)
    if (is.null(r)) { message("arm ", A$arm, ": fit failed"); next }
    dropped <- setdiff(studies_all, unique(A$d$study))
    meta <- data.table(arm = A$arm, exposure = A$x, n_donors = r$n_donors,
                       n_studies = r$n_studies, singular = r$singular,
                       studies_dropped = paste(dropped, collapse = ";"),
                       model = paste0("lmer(affinity_z ~ compartment * ", A$x,
                                      " + sex + mean_log10_total_counts",
                                      if (!is.null(A$extra)) paste0(" + ", A$extra),
                                      " + (1|study) + (1|donor_id))"))
    glob[[A$arm]] <- cbind(meta, lrt_chisq = r$lrt_chisq, lrt_df = r$lrt_df,
                           lrt_p = r$lrt_p)
    s <- r$slopes[, .(compartment, slope = X.trend, SE, df, lower.CL, upper.CL,
                      t_ratio = t.ratio, p_value = p.value)]
    s[, p_BH := p.adjust(p_value, "BH")]
    slopes[[A$arm]] <- cbind(meta[rep(1, nrow(s))], s)
    c2 <- r$contrasts[, .(contrast, estimate, SE, df, lower.CL, upper.CL,
                          t_ratio = t.ratio, p_value = p.value)]
    contr[[A$arm]] <- cbind(meta[rep(1, nrow(c2))], c2)
    message(sprintf("[%s] n=%d donors / %d studies | interaction LRT chisq=%.2f df=%d P=%.3g",
                    A$arm, r$n_donors, r$n_studies, r$lrt_chisq, r$lrt_df, r$lrt_p))
}
glob <- rbindlist(glob, fill = TRUE)
slopes <- rbindlist(slopes, fill = TRUE)
contr <- rbindlist(contr, fill = TRUE)

## ---- detection-matched null (primary arm only) --------------------------------
null_genes <- unique(gpb[role == "null", gene])
if (opt$max_null > 0) null_genes <- head(null_genes, opt$max_null)
message(sprintf("null: refitting the primary model for %d detection-matched genes",
                length(null_genes)))
nullres <- rbindlist(lapply(seq_along(null_genes), function(i) {
    g <- null_genes[i]
    if (i %% 25 == 0) message(sprintf("  null %d/%d", i, length(null_genes)))
    gd <- gpb[gene == g, .(donor_id, AGTR1_pb = expr, AGTR1_pb_detect = detect)]
    d <- merge(donor[n_cells >= opt$min_cells, !c("AGTR1_pb", "AGTR1_pb_detect")],
               gd, by = "donor_id")
    r <- fit_long(prep(d), "AGTR1_z")
    if (is.null(r)) return(NULL)
    rbind(data.table(gene = g, stat = paste0("slope_", r$slopes$compartment),
                     value = r$slopes$X.trend),
          data.table(gene = g, stat = paste0("contrast_", r$contrasts$contrast),
                     value = r$contrasts$estimate),
          data.table(gene = g, stat = "lrt_chisq", value = r$lrt_chisq))
}), fill = TRUE)
fwrite(nullres, file.path(opt$outdir, "niche_affinity_null.tsv"), sep = "\t")

attach_null <- function(dt, key, prefix) {
    dt[, `:=`(null_mean = NA_real_, null_sd = NA_real_, n_null = NA_integer_,
              p_emp = NA_real_)]
    for (i in which(dt$arm == "primary")) {
        nv <- nullres[stat == paste0(prefix, dt[[key]][i]), value]
        dt[i, `:=`(null_mean = mean(nv), null_sd = sd(nv), n_null = length(nv),
                   p_emp = emp_p(dt[[if (prefix == "slope_") "slope" else "estimate"]][i], nv))]
    }
    dt
}
slopes <- attach_null(slopes, "compartment", "slope_")
contr  <- attach_null(contr, "contrast", "contrast_")
lrt_null <- nullres[stat == "lrt_chisq", value]
glob[, `:=`(lrt_null_median = median(lrt_null), n_null = length(lrt_null),
            lrt_p_emp = fifelse(arm == "primary",
                                (1 + sum(lrt_null >= lrt_chisq[1])) / (length(lrt_null) + 1),
                                NA_real_))]

write_tsv_safe(glob,   file.path(opt$outdir, "niche_affinity_global.tsv"))
write_tsv_safe(slopes, file.path(opt$outdir, "niche_affinity_agtr1_models.tsv"))
write_tsv_safe(contr,  file.path(opt$outdir, "niche_affinity_contrasts.tsv"))
message("\n---- per-compartment AGTR1 slopes (primary) ----")
print(slopes[arm == "primary", .(compartment, slope, SE, p_value, p_BH, null_mean, p_emp)])
print(contr[arm == "primary", .(contrast, estimate, SE, p_value, null_mean, p_emp)])

## ---- reference axes (separate per-axis donor models; not in the four-axis family)
ref <- rbindlist(lapply(c(FOUR, REF), function(ax) {
    y <- paste0("affinity_", ax, "_z")
    f <- as.formula(paste(y, "~ AGTR1_z + sex + mean_log10_total_counts + (1 | study)"))
    fit <- try(suppressMessages(lmer(f, data = base)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    r <- tidy_row(fit, "AGTR1_z", "donor", "lmer(+1|study)", nrow(base),
                  uniqueN(base$donor_id), spec = "per_axis")
    r[, `:=`(axis = ax, family = fifelse(ax %in% FOUR, "four_axis", "reference"))]
}))
ref[, p_BH := p.adjust(p_value, "BH"), by = family]
write_tsv_safe(ref, file.path(opt$outdir, "niche_affinity_per_axis_models.tsv"))

## ---- cell-level support -----------------------------------------------------------
cl <- rbindlist(lapply(c(FOUR, REF), function(ax) {
    y <- paste0("affinity_", ax)
    f <- as.formula(paste(y, "~ AGTR1_expr + log10_total_counts + (1 | donor_id) + (1 | study)"))
    fit <- try(suppressMessages(lmer(f, data = cells)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    r <- tidy_row(fit, "AGTR1_expr", "cell", "lmer(+1|donor)+(1|study)", nrow(cells),
                  uniqueN(cells$donor_id), spec = "cell_support")
    r[, `:=`(axis = ax, family = fifelse(ax %in% FOUR, "four_axis", "reference"))]
}))
cl[, p_BH := p.adjust(p_value, "BH"), by = family]
write_tsv_safe(cl, file.path(opt$outdir, "niche_affinity_cell_models.tsv"))

## ---- legacy aggregate, published vs study-guarded ---------------------------------
leg <- list()
if (file.exists(opt$legacy)) {
    pub <- fread(opt$legacy)
    leg[[1]] <- data.table(source = "published (localization/01.airspace_analysis.py)",
                           model = "OLS mean_airspace_score ~ frac_AGTR1_pos + age + C(sex)",
                           estimate = pub$estimate[1], SE = pub$se[1],
                           p_value = pub$pval[1], n_donors = NA_integer_)
}
fl <- lm(mean_airspace_score ~ frac_AGTR1_pos + age + sex, data = age_ok)
co <- summary(fl)$coefficients["frac_AGTR1_pos", ]
leg[[2]] <- data.table(source = "refit here, same OLS on donors with age",
                       model = "lm(mean_airspace_score ~ frac_AGTR1_pos + age + sex)",
                       estimate = co[1], SE = co[2], p_value = co[4],
                       n_donors = nrow(age_ok))
for (x in c("frac_AGTR1_pos", "AGTR1_z")) {
    f <- as.formula(paste("mean_airspace_score ~", x,
                          "+ sex + mean_log10_total_counts + (1 | study)"))
    fit <- suppressMessages(lmer(f, data = base))
    co <- summary(fit)$coefficients[x, ]
    leg[[length(leg) + 1]] <- data.table(
        source = "study-guarded refit, no age (primary convention)",
        model = paste0("lmer(mean_airspace_score ~ ", x,
                       " + sex + mean_log10_total_counts + (1|study))"),
        estimate = co[["Estimate"]], SE = co[["Std. Error"]],
        p_value = co[["Pr(>|t|)"]], n_donors = nrow(base))
}
write_tsv_safe(rbindlist(leg, fill = TRUE),
               file.path(opt$outdir, "niche_affinity_legacy_comparison.tsv"))

cat("\nReproducibility information:\n"); print(sessionInfo())
