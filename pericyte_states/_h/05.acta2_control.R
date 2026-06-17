## CONTROL / BENCHMARK (not a mechanistic pillar): does AGTR1 simply track the
## canonical ACTA2+ contractile mural identity? We put ACTA2 through the SAME
## donor-aware framework used for AGTR1 (03.agtr1_lenses.R) and ask three things:
##
##   (A) ACTA2 across pericyte programs (donor-aware emmeans by state_program;
##       raw expr + detectability) -- the contractile-marker benchmark. Compared
##       side-by-side with the AGTR1 lenses already computed in 03.
##   (B) AGTR1 vs ACTA2 at the donor x program pseudobulk level -- if AGTR1 is just
##       contractile identity, the denoised AGTR1 should co-vary tightly with ACTA2
##       across donor-program units. Pearson + Spearman, raw-AGTR1 and denoised-AGTR1
##       vs raw ACTA2.
##   (C) Leave-ACTA2-out synthetic/contractile score (04.acta2_control.py): a broader
##       contractile program with the single focal gene removed. We (i) confirm it
##       still ranks the synthetic/contractile program top across programs, then
##       (ii) ask whether AGTR1 (raw + denoised) tracks THIS broader contractile axis
##       at donor x program pseudobulk -- testing the "contractile identity" claim
##       beyond ACTA2 alone.
##
## Interpretation guard: the robust lens for AGTR1 is the scVI-denoised one
## (AGTR1_scvi); raw AGTR1/ACTA2 are dropout-affected. A weak denoised-AGTR1 vs
## contractile correlation means AGTR1 is NOT reducible to contractile mural identity.
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr); library(ggplot2)
    library(lme4); library(lmerTest); library(emmeans)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) { i <- which(args == flag); if (length(i)) args[i + 1] else default }
META    <- parse_arg("--meta", "pericytes_states_metadata.tsv.gz")
DENOISE <- parse_arg("--denoise",
                     "../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv")
DEN_MODEL <- parse_arg("--den-model", "Pericyte-only-trained")
NOACTA2 <- parse_arg("--noacta2", "synth_contr_noACTA2.tsv.gz")
OUTDIR  <- parse_arg("--outdir", "stats_data")
MIN_CELLS <- as.integer(parse_arg("--min-cells", "5"))   # per donor x program unit
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
save_gg <- function(fn, p, w, h) for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)

## ---- load + merge (all TSV; no h5ad) ---------------------------------------
meta <- fread(META); setnames(meta, 1, "barcode")
den  <- fread(DENOISE, select = c("index", "Model", "AGTR1_scvi"))
setnames(den, "index", "barcode"); den <- den[Model == DEN_MODEL]
stopifnot(!anyDuplicated(den$barcode))
noa  <- fread(NOACTA2)                                  # barcode + synth_contr_noACTA2_score
df <- meta |>
    merge(den[, .(barcode, AGTR1_scvi)], by = "barcode") |>
    merge(noa, by = "barcode")
cat(sprintf("merged %d cells\n", nrow(df)))
df[, state_program := factor(state_program)]
df[, donor_id := factor(donor_id)]

## ===========================================================================
## (A) ACTA2 across pericyte programs -- donor-aware emmeans (benchmark)
## ===========================================================================
ACTA2_LENSES <- c(ACTA2_expr = "ACTA2_expr", ACTA2_detect = "ACTA2_detect")
emm_all <- list(); ph_all <- list()
for (nm in names(ACTA2_LENSES)) {
    col <- ACTA2_LENSES[[nm]]; if (!col %in% names(df)) next
    d <- df[is.finite(get(col))]; d[, w := get(col)]
    fit <- tryCatch(suppressMessages(lmer(w ~ state_program + (1 | donor_id), data = d)),
                    error = function(e) lm(w ~ state_program, data = d))
    e <- as.data.frame(emmeans(fit, ~ state_program)); e$lens <- nm
    ph <- as.data.frame(pairs(emmeans(fit, ~ state_program), adjust = "BH")); ph$lens <- nm
    emm_all[[nm]] <- e; ph_all[[nm]] <- ph
    cat("\n== ", nm, " by program (emmeans) ==\n", sep = ""); print(e[, c("state_program","emmean","SE")])
}
emmA <- rbindlist(emm_all, fill = TRUE); phA <- rbindlist(ph_all, fill = TRUE)
fwrite(emmA, file.path(OUTDIR, "acta2_by_program_emmeans.tsv"), sep = "\t")
fwrite(phA,  file.path(OUTDIR, "acta2_by_program_posthoc.tsv"), sep = "\t")
rankA <- emmA[, .(top_program = state_program[which.max(emmean)],
                  top_emmean = max(emmean)), by = lens]
fwrite(rankA, file.path(OUTDIR, "acta2_by_program_top.tsv"), sep = "\t")
cat("\n== ACTA2 top program per lens ==\n"); print(rankA)

## ===========================================================================
## (B,C) donor x program pseudobulk + correlations
## ===========================================================================
pb <- df[, .(AGTR1_expr = mean(AGTR1_expr, na.rm = TRUE),
             AGTR1_scvi = mean(AGTR1_scvi, na.rm = TRUE),
             ACTA2_expr = mean(ACTA2_expr, na.rm = TRUE),
             synth_contr = mean(synthetic_contractile_score, na.rm = TRUE),
             synth_contr_noACTA2 = mean(synth_contr_noACTA2_score, na.rm = TRUE),
             n = .N),
         by = .(donor_id, state_program)]
pb <- pb[n >= MIN_CELLS]
cat(sprintf("\npseudobulk: %d donor x program units (>= %d cells)\n", nrow(pb), MIN_CELLS))
fwrite(pb, file.path(OUTDIR, "acta2_control_pseudobulk.tsv"), sep = "\t")

cor_row <- function(x, y, xn, yn) {
    ok <- is.finite(x) & is.finite(y)
    pe <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    sp <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    data.table(x = xn, y = yn, n = sum(ok),
               pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
               spearman_rho = unname(sp$estimate), spearman_p = sp$p.value)
}
cors <- rbindlist(list(
    ## (B) AGTR1 vs ACTA2 directly
    cor_row(pb$AGTR1_expr, pb$ACTA2_expr, "AGTR1_expr", "ACTA2_expr"),
    cor_row(pb$AGTR1_scvi, pb$ACTA2_expr, "AGTR1_scvi", "ACTA2_expr"),
    ## (C) AGTR1 vs broader leave-ACTA2-out contractile program
    cor_row(pb$AGTR1_expr, pb$synth_contr_noACTA2, "AGTR1_expr", "synth_contr_noACTA2"),
    cor_row(pb$AGTR1_scvi, pb$synth_contr_noACTA2, "AGTR1_scvi", "synth_contr_noACTA2"),
    ## reference: AGTR1 vs full contractile score (with ACTA2), and ACTA2 vs leave-out
    cor_row(pb$AGTR1_scvi, pb$synth_contr, "AGTR1_scvi", "synth_contr_full"),
    cor_row(pb$ACTA2_expr, pb$synth_contr_noACTA2, "ACTA2_expr", "synth_contr_noACTA2")
))
cors[, pearson_p_BH := p.adjust(pearson_p, "BH")]
cors[, spearman_p_BH := p.adjust(spearman_p, "BH")]
fwrite(cors, file.path(OUTDIR, "acta2_control_correlations.tsv"), sep = "\t")
cat("\n== donor x program pseudobulk correlations ==\n"); print(cors)

## ---- (C) leave-ACTA2-out contractile score across programs (sanity) --------
d <- df[is.finite(synth_contr_noACTA2_score)]
fit <- tryCatch(suppressMessages(lmer(synth_contr_noACTA2_score ~ state_program + (1 | donor_id), data = d)),
                error = function(e) lm(synth_contr_noACTA2_score ~ state_program, data = d))
emmC <- as.data.frame(emmeans(fit, ~ state_program))
fwrite(emmC, file.path(OUTDIR, "synth_contr_noACTA2_by_program_emmeans.tsv"), sep = "\t")
cat("\n== leave-ACTA2-out contractile score by program (emmeans) ==\n")
print(emmC[, c("state_program","emmean","SE")])

## ---- figure: AGTR1 (denoised) vs ACTA2 and vs leave-out contractile --------
plt <- melt(pb, id.vars = c("donor_id","state_program","AGTR1_scvi"),
            measure.vars = c("ACTA2_expr","synth_contr_noACTA2"),
            variable.name = "contractile", value.name = "value")
pg <- ggplot(plt, aes(value, AGTR1_scvi)) +
    geom_point(aes(color = state_program), alpha = 0.6, size = 1.4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    facet_wrap(~ contractile, scales = "free_x") +
    labs(x = "contractile readout (donor x program pseudobulk)",
         y = "denoised AGTR1 (AGTR1_scvi)",
         title = "Control: denoised AGTR1 vs contractile identity",
         color = "program") +
    theme_bw(base_size = 11)
save_gg(file.path(OUTDIR, "acta2_control_agtr1_vs_contractile"), pg, 9, 4)

cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
