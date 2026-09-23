## Figure 5D -- does the pericyte AT1R-response programme track AGTR1, the
## continuum and the subclusters?
##
## Unit: donor x pericyte_state (Leiden subcluster, panel-independent), >= 5 cells,
## as in basement_membrane/_h/10.agtr1_count_models.R. Scores are z within dataset.
##
## (A) ARBITER -- AGTR1 integer counts as the RESPONSE of an NB GLMM with a
##     library-size offset (no imputation):
##       agtr1_sum ~ score_z + mean_log10_counts + (1|study) + (1|donor_id)
##                   + offset(log(total_sum))
##     Poisson + OLRE fallback, labelled. Binomial (detected / not) robustness.
##     For the signature, the SAME fit on 1,000 detection-matched null panels gives
##     p_emp referred to the null's own centre. (Standing rule: the count model
##     arbitrates; the denoised lens is a concordant sensitivity.)
## (B) DENOISED LENS -- lmer(AGTR1_scvi_z ~ score_z + depth + (1|donor)+(1|study)).
## (C) CONTINUUM -- per-donor Spearman of each score against DPT pseudotime (root:
##     vascular-stabilizing, fixed 2026-09-02), >= 20 cells, raw and partialled on
##     depth + ambient tracer; one-sample Wilcoxon on donor rhos, BH within family.
##     The signature's median donor rho is also referred to 200 null panels,
##     because every panel score except BM falls along this axis and that may be a
##     score-magnitude gradient rather than biology.
## (D) SUBCLUSTERS -- lmer(score_z ~ pericyte_state + depth + (1|donor)+(1|study)),
##     emmeans + BH pairwise; units with < 15 donors are flagged.
## (E) Donor-level Spearman with the six state-programme scores (descriptive).

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(lmerTest)
    library(emmeans); library(parallel)
})
source("../_h/_stats_common.R")
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--meta", default = "./at1r_response_metadata.tsv.gz"),
    make_option("--counts", default = "../../basement_membrane/_m/agtr1_count_input.tsv.gz"),
    make_option("--null-pb", dest = "null_pb", default = "./at1r_null_pseudobulk.tsv.gz"),
    make_option("--null-cells", dest = "null_cells", default = "./at1r_null_cells.tsv.gz"),
    make_option("--outdir", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--max-null", type = "integer", default = 0L, dest = "max_null")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
CORES <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

m <- read_req(opt$meta); setnames(m, 1, "index")
cnt <- read_req(opt$counts); setnames(cnt, 1, "index")
cnt <- cnt[, .(index, AGTR1_count, raw_total_counts)]
d <- merge(m, cnt, by = "index")
message(sprintf("cells: %d of %d matched to AGTR1 counts", nrow(d), nrow(m)))
if (nrow(d) < 0.99 * nrow(m)) stop("count input does not cover the scored cells")

SCORES <- c("at1r_response_score", "at1r_legacy_curated_score",
            grep("^(progeny_|tf_)", names(d), value = TRUE))
family_of <- function(s) fcase(
    s == "at1r_response_score", "signature",
    s == "at1r_legacy_curated_score", "legacy_curated",
    s %in% paste0("progeny_", c("MAPK", "NFkB", "TGFb", "EGFR", "JAK-STAT")), "progeny_prespecified",
    grepl("^progeny_", s), "progeny_other",
    default = "tf")

## ---- unit pseudobulk --------------------------------------------------------------
pb <- d[, c(.(n_cells = .N, agtr1_sum = sum(AGTR1_count), n_pos = sum(AGTR1_count > 0),
              total_sum = sum(raw_total_counts),
              mean_log10_counts = mean(log10_total_counts),
              AGTR1_scvi = mean(AGTR1_scvi, na.rm = TRUE),
              study = study[1], dataset = dataset[1]),
            lapply(.SD, mean, na.rm = TRUE)),
        by = .(donor_id, pericyte_state),
        .SDcols = c(SCORES, "basement_membrane_score", "fibrillar_collagen_score",
                    "fibrillar_ecm_score", "ambient_tracer_score")]
pb <- pb[n_cells >= opt$min_cells]
pb[, pericyte_state := factor(as.character(pericyte_state))]
for (s in c(SCORES, "AGTR1_scvi", "ambient_tracer_score"))
    pb[[paste0(s, "_z")]] <- z_within_dataset(pb[[s]], pb$dataset)
message(sprintf("units (>= %d cells): %d from %d donors, %d studies",
                opt$min_cells, nrow(pb), uniqueN(pb$donor_id), uniqueN(pb$study)))
fwrite(pb, "./at1r_unit_pseudobulk.tsv.gz", sep = "\t")

## ---- (A) arbiter --------------------------------------------------------------------
## Errors fall back; WARNINGS DO NOT. Gating on warnings silently discarded 18
## null genes in basement_membrane/16 -- almost all the highest-expressed ones --
## because glmer.nb warns there (memory: gene-set-score-null-not-zero). Warnings
## are muffled, counted and written to `note`; convergence is graded on max|grad|.
fit_nb <- function(pred, data, level = "pseudobulk") {
    f_nb <- as.formula(sprintf("agtr1_sum ~ %s + mean_log10_counts + (1|study) + (1|donor_id) + offset(log(total_sum))", pred))
    warn <- character(0)
    fit <- tryCatch(withCallingHandlers(suppressMessages(glmer.nb(f_nb, data = data)),
                                        warning = function(w) {
                                            warn <<- c(warn, conditionMessage(w))
                                            invokeRestart("muffleWarning")
                                        }),
                    error = function(e) NULL)
    if (!is.null(fit)) {
        g <- tryCatch(max(abs(fit@optinfo$derivs$gradient)), error = function(e) NA_real_)
        r <- tidy_row(fit, pred, level, "NB GLMM", nrow(data), uniqueN(data$donor_id),
                      note = if (length(warn)) paste(unique(warn), collapse = " | ") else "")
        if (is.null(r)) return(NULL)
        return(r[, `:=`(max_abs_grad = g, n_warnings = length(warn))])
    }
    dd <- copy(data)[, olre := factor(seq_len(.N))]
    f_po <- update(f_nb, . ~ . + (1 | olre))
    fit <- tryCatch(suppressMessages(glmer(f_po, data = dd, family = poisson)),
                    error = function(e) NULL)
    tidy_row(fit, pred, level, "Poisson+OLRE", nrow(data), uniqueN(data$donor_id),
             note = "glmer.nb did not converge")
}
arb <- rbindlist(lapply(SCORES, function(s) {
    r <- fit_nb(paste0(s, "_z"), pb)
    fb <- tryCatch(suppressMessages(glmer(as.formula(sprintf(
        "cbind(n_pos, n_cells - n_pos) ~ %s_z + mean_log10_counts + (1|study) + (1|donor_id)", s)),
        data = pb, family = binomial)), error = function(e) NULL)
    rb <- tidy_row(fb, paste0(s, "_z"), "pseudobulk", "binomial GLMM", nrow(pb),
                   uniqueN(pb$donor_id), spec = "detection")
    out <- rbindlist(list(r, rb), fill = TRUE)
    if (nrow(out)) out[, `:=`(score = s, family = family_of(s))]
    out
}), fill = TRUE)

## null for the signature
nl <- read_req(opt$null_pb)
nl[, pericyte_state := as.character(pericyte_state)]
null_cols <- grep("^null_", names(nl), value = TRUE)
if (opt$max_null > 0) null_cols <- head(null_cols, opt$max_null)
pbn <- merge(pb[, .(donor_id = as.character(donor_id),
                    pericyte_state = as.character(pericyte_state),
                    agtr1_sum, total_sum, mean_log10_counts, study, dataset,
                    AGTR1_scvi_z, basement_membrane_score, fibrillar_collagen_score)],
             nl[, c("donor_id", "pericyte_state", null_cols), with = FALSE],
             by = c("donor_id", "pericyte_state"))
if (nrow(pbn) != nrow(pb)) stop("null pseudobulk covers ", nrow(pbn), " of ", nrow(pb), " units")
for (cl in null_cols) pbn[[paste0(cl, "_z")]] <- z_within_dataset(pbn[[cl]], pbn$dataset)
message(sprintf("arbiter null: %d panels on %d cores", length(null_cols), CORES))
null_nb <- rbindlist(mclapply(null_cols, function(cl) {
    r <- fit_nb(paste0(cl, "_z"), pbn)
    if (is.null(r)) return(NULL)
    r[, panel := cl]
}, mc.cores = CORES), fill = TRUE)
null_scvi <- rbindlist(mclapply(null_cols, function(cl) {
    fit <- try(suppressMessages(lmer(as.formula(sprintf(
        "AGTR1_scvi_z ~ %s_z + mean_log10_counts + (1|donor_id) + (1|study)", cl)),
        data = pbn)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    data.table(panel = cl, estimate = fixef(fit)[[paste0(cl, "_z")]])
}, mc.cores = CORES))
fwrite(null_nb, file.path(opt$outdir, "at1r_vs_agtr1_count_null.tsv"), sep = "\t")
fwrite(null_scvi, file.path(opt$outdir, "at1r_vs_agtr1_scvi_null.tsv"), sep = "\t")

sig_row <- arb[score == "at1r_response_score" & model %in% c("NB GLMM", "Poisson+OLRE")]
arb[, `:=`(null_mean = NA_real_, null_sd = NA_real_, n_null = NA_integer_, p_emp = NA_real_)]
if (nrow(sig_row)) {
    nv <- null_nb[model == sig_row$model[1], estimate]
    arb[score == "at1r_response_score" & model == sig_row$model[1],
        `:=`(null_mean = mean(nv), null_sd = sd(nv), n_null = length(nv),
             p_emp = emp_p(estimate, nv))]
    message(sprintf("ARBITER signature: beta %.3f (model %s, P %.3g); null centre %.3f, p_emp %.3g",
                    sig_row$estimate, sig_row$model, sig_row$p_value, mean(nv),
                    emp_p(sig_row$estimate, nv)))
}
arb[, p_BH := p.adjust(p_value, "BH"), by = .(family, model)]
write_tsv_safe(arb, file.path(opt$outdir, "at1r_vs_agtr1_count.tsv"))

## ---- (B) denoised lens ----------------------------------------------------------------
scvi <- rbindlist(lapply(SCORES, function(s) {
    fit <- try(suppressMessages(lmer(as.formula(sprintf(
        "AGTR1_scvi_z ~ %s_z + mean_log10_counts + (1|donor_id) + (1|study)", s)),
        data = pb)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    r <- tidy_row(fit, paste0(s, "_z"), "pseudobulk", "lmer denoised lens", nrow(pb),
                  uniqueN(pb$donor_id), spec = "sensitivity")
    r[, `:=`(score = s, family = family_of(s))]
}), fill = TRUE)
nv <- null_scvi$estimate
## emp_p() is scalar; calling it inside fifelse() passes the whole estimate column
## and errors. Assign on the matching row instead.
scvi[, `:=`(null_mean = NA_real_, null_sd = NA_real_, p_emp = NA_real_)]
scvi[score == "at1r_response_score",
     `:=`(null_mean = mean(nv), null_sd = sd(nv), p_emp = emp_p(estimate, nv))]
scvi[, p_BH := p.adjust(p_value, "BH"), by = family]
write_tsv_safe(scvi, file.path(opt$outdir, "at1r_vs_agtr1_scvi.tsv"))

## ---- (C) continuum ---------------------------------------------------------------------
spear <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
part <- function(v, dd) resid(lm(v ~ log10_total_counts + ambient_tracer_score, data = dd,
                                 na.action = na.exclude))
cont_rows <- d[is.finite(dpt_pseudotime)]
donor_rho <- cont_rows[, {
    if (.N < 20) NULL else {
        pt_r <- part(dpt_pseudotime, .SD)
        c(list(n = .N),
          setNames(lapply(SCORES, function(s) spear(get(s), dpt_pseudotime)),
                   paste0("rho__", SCORES)),
          setNames(lapply(SCORES, function(s) spear(part(get(s), .SD), pt_r)),
                   paste0("prho__", SCORES)))
    }
}, by = donor_id]
fwrite(donor_rho, file.path(opt$outdir, "at1r_continuum_donor_rho.tsv"), sep = "\t")
cs <- rbindlist(lapply(SCORES, function(s) rbindlist(lapply(c("raw", "partial"), function(k) {
    v <- donor_rho[[paste0(if (k == "raw") "rho__" else "prho__", s)]]
    v <- v[is.finite(v)]
    if (length(v) < 5) return(NULL)
    data.table(score = s, family = family_of(s), kind = k, n_donors = length(v),
               median_rho = median(v), q25 = quantile(v, 0.25), q75 = quantile(v, 0.75),
               p_wilcox = suppressWarnings(wilcox.test(v))$p.value)
}))))
cs[, p_BH := p.adjust(p_wilcox, "BH"), by = .(family, kind)]
## signature median rho vs the null panels' median rho
if (file.exists(opt$null_cells)) {
    nc <- fread(opt$null_cells); setnames(nc, 1, "index")
    ncols <- grep("^null_", names(nc), value = TRUE)
    nd <- merge(d[, .(index, donor_id, dpt_pseudotime)], nc, by = "index")[is.finite(dpt_pseudotime)]
    keep_d <- donor_rho$donor_id
    null_med <- vapply(ncols, function(cl) {
        r <- nd[donor_id %in% keep_d, .(r = spear(get(cl), dpt_pseudotime)), by = donor_id]$r
        median(r, na.rm = TRUE)
    }, numeric(1))
    fwrite(data.table(panel = ncols, median_rho = null_med),
           file.path(opt$outdir, "at1r_continuum_null.tsv"), sep = "\t")
    i <- cs$score == "at1r_response_score" & cs$kind == "raw"
    cs[i, `:=`(null_mean = mean(null_med), null_sd = sd(null_med),
               p_emp = emp_p(median_rho, null_med))]
}
write_tsv_safe(cs, file.path(opt$outdir, "at1r_continuum_summary.tsv"))
print(cs[family %in% c("signature", "legacy_curated") ])

## ---- (D) subclusters ---------------------------------------------------------------------
nd_state <- pb[, .(n_donors = uniqueN(donor_id), n_units = .N), by = pericyte_state]
emm <-rbindlist(lapply(c("at1r_response_score", "at1r_legacy_curated_score"), function(s) {
    fit <- suppressMessages(lmer(as.formula(sprintf(
        "%s_z ~ pericyte_state + mean_log10_counts + (1|donor_id) + (1|study)", s)), data = pb))
    e <- merge(as.data.table(emmeans(fit, "pericyte_state")), nd_state, by = "pericyte_state")
    e[, `:=`(score = s, small_group = n_donors < 15)]
}), fill = TRUE)
prs <- rbindlist(lapply(c("at1r_response_score", "at1r_legacy_curated_score"), function(s) {
    fit <- suppressMessages(lmer(as.formula(sprintf(
        "%s_z ~ pericyte_state + mean_log10_counts + (1|donor_id) + (1|study)", s)), data = pb))
    pr <- as.data.table(pairs(emmeans(fit, "pericyte_state"), adjust = "BH"))
    small <- nd_state[n_donors < 15, as.character(pericyte_state)]
    pr[, `:=`(score = s,
              touches_small_group = vapply(strsplit(as.character(contrast), " - "),
                  function(x) any(sub("^pericyte_state", "", trimws(x)) %in% small), logical(1)))]
}), fill = TRUE)
write_tsv_safe(emm, file.path(opt$outdir, "at1r_by_subcluster_emmeans.tsv"))
write_tsv_safe(prs, file.path(opt$outdir, "at1r_by_subcluster_posthoc.tsv"))

## ---- (E) donor-level correlation with the state-programme scores -------------------------------
prog <- c("vascular_stabilizing_score", "inflammatory_score", "synthetic_contractile_score",
          "activated_migratory_score", "fibroblast_like_score", "basement_membrane_score")
prog <- intersect(prog, names(d))
dm <- d[, lapply(.SD, mean, na.rm = TRUE), by = donor_id,
        .SDcols = c("at1r_response_score", prog)][d[, .N, by = donor_id][N >= 10], on = "donor_id"]
vp <- rbindlist(lapply(prog, function(p) {
    ct <- suppressWarnings(cor.test(dm$at1r_response_score, dm[[p]], method = "spearman"))
    data.table(program = p, rho = unname(ct$estimate), p_value = ct$p.value,
               n_donors = nrow(dm),
               note = "descriptive; panel scores share a positive general-expression null")
}))
write_tsv_safe(vp, file.path(opt$outdir, "at1r_vs_programs.tsv"))

cat("\nReproducibility information:\n"); print(sessionInfo())
