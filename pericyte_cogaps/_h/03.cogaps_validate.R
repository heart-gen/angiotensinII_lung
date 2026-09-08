## Validate / refine the curated pericyte-state model with CoGAPS patterns.
##
## For a chosen nPatterns, asks whether the unsupervised CoGAPS programs
## corroborate the curated 6-program NVU-pattern model and where the data diverge:
##   (A) gene overlap: top pattern markers vs each curated panel (Jaccard + phyper)
##   (B) score concordance: cell-level Spearman of pattern weights vs curated scores
##   (C) state concordance: dominant CoGAPS pattern per cell vs stable Leiden state
##   (D) clinical: per-pattern AGTR1 correlation + donor-level disease association
##
## A pattern that maps cleanly to a curated program (high Jaccard + score corr)
## corroborates it; an extra data-driven pattern with no curated match flags a
## refinement the marker panels miss. Note the curated argmax program annotation
## collapsed to a single program (synthetic_contractile dominates by raw
## magnitude); this data-driven check does not depend on that annotation.
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr); library(ggplot2)
    library(lme4); library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INDIR     <- parse_arg("--indir", "../_m")
META      <- parse_arg("--meta", "../../pericyte_states/_m/pericytes_states_metadata.tsv.gz")
HVGINFO   <- parse_arg("--hvg-info", "../_m/cogaps_hvg_info.tsv.gz")
NP        <- as.integer(parse_arg("--npatterns", "7"))
TOPK      <- as.integer(parse_arg("--top-markers", "50"))
OUTDIR    <- parse_arg("--outdir", "../_m/validation")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

## All six curated NVU-pattern panels. basement_membrane was added after the
## original 5-program run; its genes ARE present in cogaps_hvg_info (panel column,
## incl. the shared fibroblast_like;basement_membrane gene), so the top-marker
## overlap test below now scores it too. grepl() substring-matches the ";"-joined
## multi-panel labels, so a gene in a shared panel counts for each of its panels.
PROGRAMS <- c("vascular_stabilizing", "inflammatory", "synthetic_contractile",
              "activated_migratory", "fibroblast_like", "basement_membrane")
INJURY   <- c("inflammatory", "fibroblast_like", "activated_migratory")

save_gg <- function(fn, p, w, h)
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    case_when(
        grepl("^Healthy", lc) ~ "Healthy",
        lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|Systemic sclerosis",
              lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
        TRUE ~ "Other")
}

## ---- Load --------------------------------------------------------------------
pat <- fread(file.path(INDIR, sprintf("patterns_np%d.tsv.gz", NP)))
pcols <- grep("^Pattern_", names(pat), value = TRUE)
cat(sprintf("loaded %d cells x %d patterns\n", nrow(pat), length(pcols)))

meta <- fread(META)
setnames(meta, 1, "barcode")
meta[, disease_group := factor(map_disease_group(lung_condition))]
meta[, age := suppressWarnings(as.numeric(age_or_mean_of_age_range))]

df <- merge(pat, meta, by = "barcode")
cat(sprintf("merged: %d cells (%.1f%% of patterns matched metadata)\n",
            nrow(df), 100 * nrow(df) / nrow(pat)))
score_cols <- grep("_score$", names(df), value = TRUE)

## ---- (A) Pattern markers vs curated panels (Jaccard + hypergeometric) --------
hvg <- fread(HVGINFO)
bg <- unique(hvg$gene); nbg <- length(bg)
panel_sets <- lapply(PROGRAMS, function(p)
    intersect(hvg$gene[grepl(p, hvg$panel)], bg))
names(panel_sets) <- PROGRAMS

pm_file <- file.path(INDIR, sprintf("pattern_markers_np%d.tsv.gz", NP))
if (file.exists(pm_file)) {
    pm <- fread(pm_file)
    overlap <- rbindlist(lapply(unique(pm$pattern), function(pp) {
        markers <- head(pm[pattern == pp]$gene, TOPK)
        rbindlist(lapply(PROGRAMS, function(prog) {
            ps <- panel_sets[[prog]]; ov <- length(intersect(markers, ps))
            data.table(pattern = pp, program = prog,
                       n_markers = length(markers), panel_size = length(ps),
                       overlap = ov,
                       jaccard = ov / length(union(markers, ps)),
                       phyper = phyper(ov - 1, length(ps), nbg - length(ps),
                                       length(markers), lower.tail = FALSE))
        }))
    }))
    overlap[, padj := p.adjust(phyper, "BH")]
    fwrite(overlap, file.path(OUTDIR, "pattern_panel_overlap.tsv.gz"), sep = "\t")
    pg <- ggplot(overlap, aes(program, pattern, fill = jaccard)) +
        geom_tile() +
        geom_text(aes(label = sprintf("%d\n%.0e", overlap, phyper)), size = 2.6) +
        scale_fill_viridis_c() +
        labs(title = sprintf("Top-%d pattern markers vs curated panels", TOPK),
             x = "", y = "") +
        theme_bw(base_size = 11) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_gg(file.path(OUTDIR, "pattern_panel_overlap"), pg, 7, 5)
}

## ---- (B) Cell-level Spearman: pattern weights vs curated scores --------------
cormat <- outer(pcols, score_cols, Vectorize(function(a, b)
    suppressWarnings(cor(df[[a]], df[[b]], method = "spearman"))))
dimnames(cormat) <- list(pcols, score_cols)
cor_long <- as.data.table(as.table(cormat))
setnames(cor_long, c("pattern", "score", "rho"))
fwrite(cor_long, file.path(OUTDIR, "pattern_score_spearman.tsv.gz"), sep = "\t")
cg <- ggplot(cor_long, aes(gsub("_score$", "", score), pattern, fill = rho)) +
    geom_tile() + geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Pattern weight vs curated state score (Spearman)",
         x = "", y = "") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_gg(file.path(OUTDIR, "pattern_score_spearman"), cg, 7, 5)

## ---- (C) Dominant pattern per cell vs stable Leiden state --------------------
df[, dominant_pattern := pcols[max.col(as.matrix(.SD), ties.method = "first")],
   .SDcols = pcols]
ct <- df[, .N, by = .(pericyte_state, dominant_pattern)]
ct[, frac := N / sum(N), by = pericyte_state]
fwrite(ct, file.path(OUTDIR, "state_vs_dominant_pattern.tsv.gz"), sep = "\t")
sg <- ggplot(ct, aes(factor(pericyte_state), dominant_pattern, fill = frac)) +
    geom_tile() + geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
    scale_fill_viridis_c() +
    labs(title = "Stable Leiden state -> dominant CoGAPS pattern",
         x = "pericyte_state (stable cluster)", y = "") +
    theme_bw(base_size = 11)
save_gg(file.path(OUTDIR, "state_vs_dominant_pattern"), sg, 7, 5)

## ---- (D) Per-pattern AGTR1 correlation + donor-level disease association -----
if ("AGTR1_expr" %in% names(df)) {
    agtr1 <- rbindlist(lapply(pcols, function(p)
        data.table(pattern = p,
                   rho_AGTR1 = suppressWarnings(
                       cor(df[[p]], df$AGTR1_expr, method = "spearman")))))
    fwrite(agtr1, file.path(OUTDIR, "pattern_AGTR1_spearman.tsv.gz"), sep = "\t")
}

## donor-level mean pattern weight ~ disease
##
## WHY THIS MODEL NO LONGER CARRIES `+ age`, AND WHY IT GAINED `(1 | study)`
## (changed 2026-09-07, defect P2-16 in writings/TODO.md)
##
## This fit used to be `lm(pattern ~ disease_group + age + sex)`. Both halves of
## that were wrong, and fixing either alone would have been worse than fixing
## neither -- the lesson recorded in pericyte_states/_h/01.state_stats.R (P1-2).
##
##   1. `+ age` was never a covariate adjustment. Age missingness in the HLCA is
##      a STUDY property (17 of 18 studies are all-or-nothing), so requiring it
##      deleted whole cohorts and did so unevenly across disease groups:
##      69 donors -> 36, Fibrotic/ILD 17 -> 4, Other 20 -> 3. The survivors were
##      concentrated in single studies (3 of 4 fibrotic from Lafyatis_2019; all 3
##      "Other" from Budinger_2020, which contributes no Healthy donor here), so
##      the disease contrast was confounded with study by construction.
##
##   2. There was no study term. Dropping `age` restores study-clustered donors,
##      and without a guard those donors can MANUFACTURE an effect --
##      pericyte_states saw exactly that at P = 0.0015.
##
## niche_index is the precedent for what to expect: its Healthy-vs-Other contrast
## was P = 0.018 on this design (3 donors, one study) and P = 0.25 once both
## halves were fixed.
##
## Three arms are written, labelled in an `arm` column. Quote `primary`.
donor <- df[, c(lapply(.SD, mean, na.rm = TRUE),
                .(disease_group = first(disease_group), sex = first(sex),
                  study = first(study), age = mean(age, na.rm = TRUE), n = .N)),
            by = donor_id, .SDcols = pcols]
donor <- donor[n >= 20]
donor[, disease_group := relevel(droplevels(factor(disease_group)), "Healthy")]
donor[, study := factor(study)]

cat(sprintf("\n(D) disease donors (>=20 cells): %d across %d studies\n",
            nrow(donor), nlevels(droplevels(donor$study))))
print(table(donor$disease_group))
cat(sprintf("(D) age-complete subset would be %d of %d donors:\n",
            sum(!is.na(donor$age)), nrow(donor)))
print(table(donor$disease_group[!is.na(donor$age)]))

## Fit one arm and return tidy disease-term rows. `guard` adds (1 | study);
## `ageadj` restricts to age-reporting donors and adds age as a covariate.
fit_arm <- function(p, arm, guard, ageadj) {
    sub <- donor[is.finite(get(p)) & !is.na(sex)]
    if (ageadj) sub <- sub[!is.na(age)]
    sub <- sub[, disease_group := droplevels(disease_group)]
    if (nlevels(sub$disease_group) < 2) return(NULL)
    terms <- c("disease_group", "sex", if (ageadj) "age")
    n_stud <- dplyr::n_distinct(sub$study)
    ## A random intercept needs >= 2 levels; fall back and say so if it does not.
    use_guard <- guard && n_stud > 1
    if (use_guard) {
        fit <- suppressMessages(lmerTest::lmer(
            reformulate(c(terms, "(1 | study)"), p), data = sub))
        co <- as.data.frame(summary(fit)$coefficients)
        est <- co[, "Estimate"]; pv <- co[, "Pr(>|t|)"]
        sd_study <- as.data.frame(lme4::VarCorr(fit))$sdcor[1]
        sing <- lme4::isSingular(fit)
    } else {
        fit <- lm(reformulate(terms, p), data = sub)
        co <- as.data.frame(summary(fit)$coefficients)
        est <- co[, 1]; pv <- co[, 4]
        sd_study <- NA_real_; sing <- NA
    }
    keep <- grepl("disease_group", rownames(co))
    if (!any(keep)) return(NULL)
    data.table(pattern = p, term = rownames(co)[keep],
               estimate = est[keep], p_value = pv[keep],
               arm = arm, n_donors = nrow(sub), n_studies = n_stud,
               study_sd = sd_study, singular = sing,
               n_donors_group = paste(names(table(sub$disease_group)),
                                      as.integer(table(sub$disease_group)),
                                      sep = "=", collapse = ";"),
               model = if (use_guard)
                   paste0("lmer(", paste(terms, collapse = " + "), " + (1 | study))")
               else paste0("lm(", paste(terms, collapse = " + "), ")"))
}

ARMS <- list(
    list(arm = "primary", guard = TRUE,  ageadj = FALSE),
    list(arm = "unguarded -- comparison only", guard = FALSE, ageadj = FALSE),
    list(arm = "_ageadj -- RESTRICTED to age-reporting cohorts, not an age adjustment",
         guard = TRUE, ageadj = TRUE))
dis_res <- rbindlist(lapply(ARMS, function(a)
    rbindlist(lapply(pcols, fit_arm, arm = a$arm, guard = a$guard,
                     ageadj = a$ageadj), fill = TRUE)), fill = TRUE)
if (nrow(dis_res)) {
    ## BH within arm -- the three arms are not one family.
    dis_res[, padj := p.adjust(p_value, "BH"), by = arm]
    fwrite(dis_res, file.path(OUTDIR, "pattern_disease_ols.tsv.gz"), sep = "\t")
    cat("\n(D) disease terms by arm (quote `primary`):\n")
    print(dis_res[padj < 0.05, .(pattern, term, estimate, p_value, padj, arm)])
    if (!nrow(dis_res[padj < 0.05])) cat("  none reach BH < 0.05 in any arm\n")
}

cat("\n---- sessionInfo ----\n")
sessioninfo::session_info()
