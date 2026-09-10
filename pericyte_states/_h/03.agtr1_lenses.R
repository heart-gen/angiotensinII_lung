## AGTR1 reported separately from the injury composite, via three lenses.
##
## Design point (revision): AGTR1 is the focal receptor and is dropout-prone in
## pericytes, so it is kept OUT of the injury-stromal composite (see niche_index/)
## and instead reported on its own. We ask one question -- does AGTR1 mark a
## particular pericyte program? -- through three increasingly dropout-robust lenses:
##   (1) AGTR1_expr     mean log-normalized expression (raw; dropout-affected)
##   (2) AGTR1_scvi     scVI-denoised AGTR1 (from localization/airspace_analysis;
##                      ambient/dropout-corrected -- the robust readout)
##   (3) AGTR1_detect   binary detectability fraction (SENSITIVITY only)
##
## If the SAME program ranks top across all three lenses (esp. the denoised one),
## the "AGTR1 marks the vascular-stabilizing / vulnerable mural pole" claim is not
## a dropout artifact. Each lens: donor-aware emmeans by state_program (lmer with
## donor random intercept) + pairwise BH; plus a donor-level disease association.
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr); library(ggplot2)
    library(lme4); library(lmerTest); library(emmeans)
})

## The by-program fits are CELL-level (11,680 pericytes), which is above emmeans'
## default 3000-observation ceiling, so d.f. calculation was being disabled and the
## emmeans/pairs output fell back to asymptotic (z) inference with df = Inf. Raise
## the ceiling so Satterthwaite d.f. are actually computed.
##
## pbkrtest.limit is deliberately NOT raised. emmeans prefers Kenward-Roger for
## lmerMod and only falls back to Satterthwaite above that limit; KR on 11,680 rows
## builds the full adjusted covariance and is the expensive path for no gain here
## (one 3-level fixed factor, one donor intercept). Leaving pbkrtest.limit at its
## default keeps the Satterthwaite route while lmerTest.limit makes it usable.
emm_options(lmerTest.limit = 50000)

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
META    <- parse_arg("--meta", "pericytes_states_metadata.tsv.gz")
DENOISE <- parse_arg("--denoise",
                     "../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv")
## the denoising TSV carries two scVI models per cell; use the pericyte-specific
## one (trained on the same population) as the denoised lens for pericyte states.
DEN_MODEL <- parse_arg("--den-model", "Pericyte-only-trained")
OUTDIR  <- parse_arg("--outdir", "stats_data")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
save_gg <- function(fn, p, w, h) for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)

map_disease_group <- function(lc) {
    lc <- as.character(lc)
    case_when(grepl("^Healthy", lc) ~ "Healthy",
              lc %in% c("COPD") ~ "COPD",
              grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc, ignore.case = TRUE) ~ "Fibrotic_ILD",
              TRUE ~ "Other")
}

## ---- load + merge (all TSV; no h5ad) ---------------------------------------
meta <- fread(META); setnames(meta, 1, "barcode")
den  <- fread(DENOISE, select = c("index", "Model", "AGTR1_scvi"))
setnames(den, "index", "barcode")
den  <- den[Model == DEN_MODEL]                       # one model -> unique barcodes
stopifnot(!anyDuplicated(den$barcode))
cat(sprintf("denoised model: %s (%d cells)\n", DEN_MODEL, nrow(den)))
df <- merge(meta, den[, .(barcode, AGTR1_scvi)], by = "barcode")
cat(sprintf("merged %d cells (%.1f%% of pericyte states have a denoised AGTR1)\n",
            nrow(df), 100 * nrow(df) / nrow(meta)))
df[, disease_group := relevel(factor(map_disease_group(lung_condition)), "Healthy")]
df[, age := suppressWarnings(as.numeric(age_or_mean_of_age_range))]
df[, state_program := factor(state_program)]
df[, donor_id := factor(donor_id)]

LENSES <- c(AGTR1_expr = "AGTR1_expr", AGTR1_scvi = "AGTR1_scvi", AGTR1_detect = "AGTR1_detect")

## ---- (A) by program: donor-aware emmeans + pairwise BH ---------------------
emm_all <- list(); ph_all <- list()
for (nm in names(LENSES)) {
    col <- LENSES[[nm]]; if (!col %in% names(df)) next
    d <- df[is.finite(get(col))]; d[, w := get(col)]
    fit <- tryCatch(suppressMessages(lmer(w ~ state_program + (1 | donor_id), data = d)),
                    error = function(e) lm(w ~ state_program, data = d))
    e <- as.data.frame(emmeans(fit, ~ state_program)); e$lens <- nm
    ph <- as.data.frame(pairs(emmeans(fit, ~ state_program), adjust = "BH")); ph$lens <- nm
    emm_all[[nm]] <- e; ph_all[[nm]] <- ph
    cat("\n== ", nm, " by program (emmeans) ==\n", sep = ""); print(e[, c("state_program","emmean","SE")])
}
emm <- rbindlist(emm_all, fill = TRUE); ph <- rbindlist(ph_all, fill = TRUE)
fwrite(emm, file.path(OUTDIR, "agtr1_lenses_by_program_emmeans.tsv"), sep = "\t")
fwrite(ph,  file.path(OUTDIR, "agtr1_lenses_by_program_posthoc.tsv"), sep = "\t")

## rank table: which program is top for each lens
rank_tab <- emm[, .(top_program = state_program[which.max(emmean)],
                    top_emmean = max(emmean)), by = lens]
fwrite(rank_tab, file.path(OUTDIR, "agtr1_lenses_top_program.tsv"), sep = "\t")
cat("\n== Top program per lens ==\n"); print(rank_tab)

## ---- (B) donor-level disease association per lens --------------------------
##
## FIXED 2026-09-08 (P1-23). This section used to fit
##
##     sub <- dd[!is.na(age) & !is.na(sex)]
##     lm(val ~ disease_group + age + sex, data = sub)
##
## which is the `+ age` cohort filter (P1-2), with no study term, on the endpoint
## the filter distorts most -- a disease contrast. Age missingness in the HLCA is
## a STUDY property, not a donor property, so `!is.na(age)` deletes whole cohorts
## rather than thinning uniformly: the shipped table read df = 31, i.e. 36 donors,
## against 89 in the age-free, study-guarded fits elsewhere in this project.
##
## This was the NINTH instance repo-wide and the last one unfixed. It survived
## the P1-2 sweep, the P1-8 sweep and the 2026-09-07 reconciliation because
## nothing quotes its output -- no number from `agtr1_lenses_disease_emmeans.tsv`
## appears in any .md/.tex/.txt and it is in no supplementary table part. Being
## uncited is not being correct; it is why it went unnoticed.
##
## The fix is the same shape as P1-22 and P2-16: `age` leaves the primary and
## comes back as a LABELLED `_ageadj` arm, and `(1 | study)` guards the primary.
## Both arms run through one function so they cannot drift apart.
COVARS_PRIMARY  <- c("disease_group", "sex", "(1 | study)")
COVARS_AGE_SENS <- c("disease_group", "sex", "age", "(1 | study)")
MIN_CELLS_DONOR <- 20L   # unchanged; a cell-count floor, unrelated to the above

fit_disease <- function(covars, data) {
    f <- reformulate(covars, "val")
    if (any(grepl("\\|", covars)) && uniqueN(data$study) >= 2)
        list(fit = suppressMessages(lmerTest::lmer(f, data = data)),
             model = paste0("lmer(", deparse(f[[3]]), ")"))
    else {
        ## Fall back rather than fail, but SAY SO in the model column: a
        ## single-study cohort cannot estimate a study term, and silently
        ## dropping the guard is how P1-9 shipped an unguarded fit as primary.
        f2 <- reformulate(setdiff(covars, "(1 | study)"), "val")
        list(fit = lm(f2, data = data),
             model = paste0("lm(", deparse(f2[[3]]), ") -- NO STUDY GUARD (<2 studies)"))
    }
}

run_disease_arm <- function(covars, arm_sfx, drop_age_rows) {
    emm_rows <- list(); ph_rows <- list()
    for (nm in names(LENSES)) {
        col <- LENSES[[nm]]; if (!col %in% names(df)) next
        dd <- df[is.finite(get(col)), .(val = mean(get(col), na.rm = TRUE),
                                        disease_group = first(disease_group),
                                        study = first(study),
                                        age = mean(age, na.rm = TRUE), sex = first(sex),
                                        n = .N), by = donor_id]
        dd <- dd[n >= MIN_CELLS_DONOR]
        sub <- if (drop_age_rows) dd[!is.na(age) & !is.na(sex)] else dd[!is.na(sex)]
        sub[, disease_group := relevel(droplevels(factor(disease_group)), "Healthy")]
        if (nlevels(sub$disease_group) < 2) next

        r   <- fit_disease(covars, sub)
        fit <- r$fit
        vc  <- if (inherits(fit, "merMod")) as.data.frame(lme4::VarCorr(fit)) else NULL
        ann <- function(x) {
            x$lens <- nm; x$arm <- if (nzchar(arm_sfx)) sub("^_", "", arm_sfx) else "primary"
            x$model <- r$model; x$n_donors <- nrow(sub); x$n_studies <- uniqueN(sub$study)
            x$study_sd <- if (is.null(vc)) NA_real_ else vc$sdcor[vc$grp == "study"][1]
            x$singular <- if (inherits(fit, "merMod")) lme4::isSingular(fit) else NA
            x$n_donors_group <- paste(names(table(sub$disease_group)),
                                      as.integer(table(sub$disease_group)),
                                      sep = "=", collapse = ";")
            x
        }
        e  <- ann(as.data.frame(emmeans(fit, ~ disease_group)))
        ## The old output was emmeans only -- three marginal means and no test,
        ## which invites reading a difference off overlapping CIs (the P2-10
        ## defect). Emit the BH-adjusted pairwise contrasts alongside them.
        ph <- ann(as.data.frame(pairs(emmeans(fit, ~ disease_group), adjust = "BH")))
        emm_rows[[nm]] <- e; ph_rows[[nm]] <- ph
        cat(sprintf("[%s%s] %d donors, %d studies, %s\n", nm,
                    if (nzchar(arm_sfx)) arm_sfx else "", nrow(sub),
                    uniqueN(sub$study), r$model))
    }
    list(emm = if (length(emm_rows)) rbindlist(emm_rows, fill = TRUE) else NULL,
         ph  = if (length(ph_rows))  rbindlist(ph_rows,  fill = TRUE) else NULL)
}

## WHICH LENS MAY CARRY A CONTRAST -- checked by the run, not remembered.
##
## Standing rule from [scvi-denoised-agtr1-invalid] and the `citable_as` column
## of `basement_membrane/_h/14.agtr1_null_models.R`: the scVI-denoised lens is
## CONCORDANT SENSITIVITY ONLY. It may never be the sole support for a contrast
## and is never the arbiter -- the count model is. Removing the `+ age` filter
## here (P1-23) made that rule load-bearing rather than decorative, so it is
## enforced in the output instead of being left to whoever reads the table.
mark_citable <- function(ph) {
    if (is.null(ph)) return(ph)
    sig <- ph[, .(any_sig = any(p.value < 0.05, na.rm = TRUE)), by = lens]
    raw_ok <- any(sig[lens %in% c("AGTR1_expr", "AGTR1_detect"), any_sig])
    ph[, citable_as := fifelse(
        lens != "AGTR1_scvi", "citable",
        fifelse(raw_ok,
                "CONCORDANT SENSITIVITY -- a raw-scale lens agrees; still not the arbiter",
                paste("NOT CITABLE as evidence: the denoised lens is the ONLY lens that",
                      "moves. Raw expression and detection are both null on this contrast,",
                      "and the count model (disease_association mean_expr -- that module is",
                      "now in heart-gen/lung-pericyte-analysis -- z_AGTR1 in",
                      "Pericytes, 77 donors) is a well-powered null: -0.152, BH = 0.85.",
                      "A denoiser-only result is the signature the standing rule exists",
                      "to catch.")))]
    ph[]
}

for (arm in list(list(COVARS_PRIMARY,  "",         FALSE),
                 list(COVARS_AGE_SENS, "_ageadj",  TRUE))) {
    out <- run_disease_arm(arm[[1]], arm[[2]], arm[[3]])
    out$ph <- mark_citable(out$ph)
    if (!is.null(out$emm))
        fwrite(out$emm, file.path(OUTDIR, paste0("agtr1_lenses_disease_emmeans", arm[[2]], ".tsv")), sep = "\t")
    if (!is.null(out$ph))
        fwrite(out$ph,  file.path(OUTDIR, paste0("agtr1_lenses_disease_posthoc",  arm[[2]], ".tsv")), sep = "\t")

    ## Print the rule when it fires, so a run that changes this cannot be silent.
    if (!is.null(out$ph)) {
        den <- out$ph[lens == "AGTR1_scvi" & p.value < 0.05]
        oth <- out$ph[lens != "AGTR1_scvi" & p.value < 0.05]
        if (nrow(den) && !nrow(oth))
            cat(sprintf(paste0("\n!! STANDING RULE [%s arm]: %d denoised contrast(s) reach BH < 0.05 ",
                              "while NO raw or detection contrast does.\n",
                              "   Do not report this as a disease effect on AGTR1. See the ",
                              "citable_as column.\n"),
                        if (nzchar(arm[[2]])) sub("^_", "", arm[[2]]) else "primary", nrow(den)))
    }
}

## ---- figure: by-program emmeans across the three lenses --------------------
emm[, lens := factor(lens, levels = names(LENSES))]
pg <- ggplot(emm, aes(reorder(state_program, emmean), emmean)) +
    geom_col(fill = "#3B6FB6") +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.3) +
    coord_flip() + facet_wrap(~ lens, scales = "free_x") +
    labs(x = "", y = "donor-aware emmean",
         title = "AGTR1 across pericyte programs, three lenses") +
    theme_bw(base_size = 11)
save_gg(file.path(OUTDIR, "agtr1_lenses_by_program"), pg, 9, 4)

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
