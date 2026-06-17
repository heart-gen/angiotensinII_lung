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
dis_rows <- list()
for (nm in names(LENSES)) {
    col <- LENSES[[nm]]; if (!col %in% names(df)) next
    dd <- df[is.finite(get(col)), .(val = mean(get(col), na.rm = TRUE),
                                    disease_group = first(disease_group),
                                    age = mean(age, na.rm = TRUE), sex = first(sex),
                                    n = .N), by = donor_id]
    dd <- dd[n >= 20]; dd[, disease_group := relevel(droplevels(disease_group), "Healthy")]
    sub <- dd[!is.na(age) & !is.na(sex)]
    if (nlevels(droplevels(sub$disease_group)) < 2) next
    fit <- lm(val ~ disease_group + age + sex, data = sub)
    e <- as.data.frame(emmeans(fit, ~ disease_group)); e$lens <- nm
    dis_rows[[nm]] <- e
}
if (length(dis_rows)) {
    fwrite(rbindlist(dis_rows, fill = TRUE),
           file.path(OUTDIR, "agtr1_lenses_disease_emmeans.tsv"), sep = "\t")
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
