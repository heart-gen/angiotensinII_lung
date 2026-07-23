## Where does the local lung renin-angiotensin system actually live?
##
## Two questions the collaborator asked, in one analysis:
##   (1) which cell types produce AGT, and which carry AGTR1 -- i.e. the source
##       and receiver ends of the angiotensin axis;
##   (2) how AGT relates to the other RAS components across cell types.
##
## The motivating problem: LIANA and NicheNet both score AGT -> AGTR1 as if it
## were a direct ligand-receptor pair, but angiotensinogen is a SUBSTRATE. It has
## to be cleaved by renin to Ang I and by ACE to Ang II before AGTR1 can be
## engaged. So an AGT -> AGTR1 "interaction" is only mechanistically meaningful if
## the processing machinery is present somewhere in the same tissue. Mapping which
## cell type carries which step is therefore not a side analysis -- it is what
## licenses or qualifies the interaction call.
##
## Unit of analysis is the donor x cell-type pseudobulk, dataset-standardized,
## with (1|study) and depth as a covariate, matching the BM selectivity module.

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(data.table)
    library(lme4)
    library(lmerTest)
    library(emmeans)
})
emm_options(lmerTest.limit = 50000, pbkrtest.limit = 50000)

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", type = "character"),
    make_option("--panels", type = "character"),
    make_option("--outdir", type = "character"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--min-donors", type = "integer", default = 5L, dest = "min_donors"),
    make_option("--detect-thr", type = "double", default = 0.05, dest = "detect_thr")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(x, file) {
    if (inherits(x, "emmGrid")) x <- as.data.frame(x)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)
}

pb <- fread(opt$pseudobulk)
panel <- fread(opt$panels)
pb <- pb[n_cells >= opt$min_cells]
keep <- pb[, .(nd = uniqueN(donor_id)), by = ccc_group][nd >= opt$min_donors, ccc_group]
pb <- pb[ccc_group %in% keep]
message(sprintf("Units: %d; cell types: %d; donors: %d",
                nrow(pb), uniqueN(pb$ccc_group), uniqueN(pb$donor_id)))

z_within_dataset <- function(x, g) {
    out <- numeric(length(x)); g <- as.character(g)
    for (lev in unique(g)) {
        ix <- which(g == lev); v <- x[ix]
        s <- stats::sd(v, na.rm = TRUE)
        out[ix] <- if (is.na(s) || s == 0) 0 else (v - mean(v, na.rm = TRUE)) / s
    }
    out
}

genes <- intersect(panel$gene,
                   sub("__expr$", "", grep("__expr$", names(pb), value = TRUE)))
message("Modelling ", length(genes), " RAS/comparator genes")

## ------------------------------------------------- per-gene cell-type map ----
profile <- rbindlist(lapply(genes, function(g) {
    ecol <- paste0(g, "__expr"); dcol <- paste0(g, "__detect")
    d <- copy(pb)
    d[, y := get(ecol)]
    d[, y_z := z_within_dataset(y, dataset)]
    fit <- try(suppressMessages(lmer(
        y_z ~ ccc_group + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- as.data.frame(emmeans(fit, specs = "ccc_group"))
    det <- d[, .(detect = mean(get(dcol), na.rm = TRUE),
                 expr_raw = mean(get(ecol), na.rm = TRUE)), by = ccc_group]
    m <- merge(as.data.table(e), det, by = "ccc_group")
    m[, gene := g][]
}), fill = TRUE)
write_tsv_safe(profile, file.path(opt$outdir, "ras_celltype_profile.tsv"))

## Rank each cell type per gene -- the readable "who makes what" table.
profile[, rank_in_gene := frank(-emmean, ties.method = "min"), by = gene]
top <- profile[rank_in_gene <= 3][order(gene, rank_in_gene),
                                  .(gene, ccc_group, emmean, detect, rank_in_gene)]
write_tsv_safe(top, file.path(opt$outdir, "ras_top_celltypes.tsv"))

## --------------------------------- AGT source and AGTR1 receiver contrasts ----
contrast_vs <- function(gene, ref) {
    ecol <- paste0(gene, "__expr")
    if (!ecol %in% names(pb) || !ref %in% pb$ccc_group) return(NULL)
    d <- copy(pb)
    d[, y := get(ecol)]
    d[, y_z := z_within_dataset(y, dataset)]
    d[, ccc_group := relevel(factor(ccc_group), ref = ref)]
    fit <- try(suppressMessages(lmer(
        y_z ~ ccc_group + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- emmeans(fit, specs = "ccc_group")
    ct <- as.data.frame(contrast(e, "trt.vs.ctrl", ref = ref, adjust = "BH"))
    ct$estimate <- -ct$estimate
    data.frame(gene = gene, reference = ref, ct, row.names = NULL)
}
## AGT against its top source; AGTR1 against pericytes.
agt_src <- contrast_vs("AGT", "Vascular smooth muscle")
agtr1_rec <- contrast_vs("AGTR1", "Pericytes")
write_tsv_safe(rbindlist(list(agt_src, agtr1_rec), fill = TRUE),
               file.path(opt$outdir, "agt_source_posthoc.tsv"))

## ------------------------------------------ RAS circuit completeness ----------
## Is any single cell type capable of running the axis alone? Score presence of
## each step by donor-level detection above a threshold.
steps <- list(
    substrate = "AGT",
    renin_step = "REN",
    ace_step = "ACE",
    chymase_step = c("CMA1", "CTSG"),
    receptor_AT1 = "AGTR1")
det_by_group <- profile[, .(detect = mean(detect, na.rm = TRUE)),
                        by = .(ccc_group, gene)]
comp <- rbindlist(lapply(names(steps), function(s) {
    gg <- intersect(steps[[s]], det_by_group$gene)
    if (!length(gg)) return(NULL)
    det_by_group[gene %in% gg, .(step = s, detect = max(detect, na.rm = TRUE)),
                 by = ccc_group]
}), fill = TRUE)
comp_wide <- dcast(comp, ccc_group ~ step, value.var = "detect", fill = 0)
step_names <- setdiff(names(comp_wide), "ccc_group")
comp_wide[, n_steps_present := rowSums(.SD >= opt$detect_thr), .SDcols = step_names]
## "Autonomous" requires a substrate, SOME Ang II-generating protease, and AT1R.
comp_wide[, has_substrate := substrate >= opt$detect_thr]
comp_wide[, has_protease := (renin_step >= opt$detect_thr) |
              (ace_step >= opt$detect_thr) | (chymase_step >= opt$detect_thr)]
comp_wide[, has_receptor := receptor_AT1 >= opt$detect_thr]
comp_wide[, autonomous_circuit := has_substrate & has_protease & has_receptor]
setorder(comp_wide, -n_steps_present)
write_tsv_safe(comp_wide, file.path(opt$outdir, "ras_circuit_completeness.tsv"))

n_auto <- sum(comp_wide$autonomous_circuit, na.rm = TRUE)
renin_max <- det_by_group[gene == "REN", max(detect, na.rm = TRUE)]

readme <- c(
    "Local RAS landscape -- generated summary",
    sprintf("Units (>=%d cells): %d; cell types: %d; donors: %d",
            opt$min_cells, nrow(pb), uniqueN(pb$ccc_group), uniqueN(pb$donor_id)),
    sprintf("Detection threshold for 'step present': %.2f", opt$detect_thr),
    "",
    sprintf("Cell types with an autonomous AGT->AngII->AT1R circuit: %d", n_auto),
    sprintf("Maximum REN (renin) detection across all cell types: %.4f", renin_max),
    "",
    "Top 3 cell types per gene:",
    paste(utils::capture.output(print(top)), collapse = "\n"),
    "",
    "Circuit completeness by cell type:",
    paste(utils::capture.output(print(comp_wide)), collapse = "\n"))
writeLines(readme, file.path(opt$outdir, "ras_landscape_README.txt"))
message(paste(readme, collapse = "\n"))

sessioninfo::session_info()
