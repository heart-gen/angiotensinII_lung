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

## ---------------------------------------------------------------------------
## emmeans' trt.vs.ctrl returns (other - REF); every table here reports
## (REF - other), so the WHOLE contrast has to be reversed -- estimate,
## t.ratio and the label together.
##
## The idiom this replaces negated `estimate` alone. p-values are two-sided so
## they stayed correct, and SE/df were never affected, but `t.ratio` kept
## pointing the other way in 882 of 882 shipped rows across four supplementary
## tables: a reader inferring direction from `t.ratio` got every comparison
## backwards. Flipping in one place is what stops the three from drifting apart
## again.
flip_contrast <- function(ct, ref) {
    ct$estimate <- -ct$estimate
    ## emmeans names this z.ratio when df are infinite; handle both.
    if ("t.ratio" %in% names(ct)) ct$t.ratio <- -ct$t.ratio
    if ("z.ratio" %in% names(ct)) ct$z.ratio <- -ct$z.ratio
    ct$contrast <- paste0(ref, " - ",
                          sub(" - .*$", "", as.character(ct$contrast)))
    ct
}
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
## A stratum defined by a gene cannot be scored for that gene. The AT2 groups
## are split on AGTR2 detectability
## (cell_communication/_h/00.prepare_ccc_input.py:164,
## `np.where(agtr2 > 0, "AT2_AGTR2det", "AT2_AGTR2undet")`), so their AGTR2
## detection is 1.000 and 0.000 by construction, not by measurement.
##
## Suppressing the two ROWS afterwards is not enough, which is why the filter
## sits here rather than at the write step. Expression is z-scored within
## dataset across every unit BEFORE the model is fit, and the det stratum's raw
## AGTR2 (1.746) is 117x the next-highest cell type (peribronchial fibroblasts,
## 0.015). It therefore dominated the standardizing SD and squeezed all 21
## genuine populations into a 0.19-z band -- the scale the within-gene ranking
## is read off was set by the artifact. The circular units have to leave before
## z_within_dataset(), not after emmeans().
##
## The rows are kept, not deleted: `detect` and `expr_raw` are still the honest
## description of those strata, and a reader who does not know how the groups
## were built needs to see why the model declined to score them. `emmean` and
## its interval are NA, and `circular_by_construction` says so.
CIRCULAR_UNITS <- list(AGTR2 = c("AT2_AGTR2det", "AT2_AGTR2undet"))

profile <- rbindlist(lapply(genes, function(g) {
    ecol <- paste0(g, "__expr"); dcol <- paste0(g, "__detect")
    circ <- intersect(CIRCULAR_UNITS[[g]], unique(pb$ccc_group))
    d <- copy(pb)
    det <- d[, .(detect = mean(get(dcol), na.rm = TRUE),
                 expr_raw = mean(get(ecol), na.rm = TRUE)), by = ccc_group]
    if (length(circ)) {
        message(sprintf("  %s: excluding %d unit(s) defined by this gene (%s)",
                        g, length(circ), paste(circ, collapse = ", ")))
        d <- d[!ccc_group %in% circ]
        d[, ccc_group := droplevels(factor(ccc_group))]
    }
    d[, y := get(ecol)]
    d[, y_z := z_within_dataset(y, dataset)]
    fit <- try(suppressMessages(lmer(
        y_z ~ ccc_group + mean_log10_total_counts + (1 | donor_id) + (1 | study),
        data = d)), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    e <- as.data.frame(emmeans(fit, specs = "ccc_group"))
    ## merge = all.y keeps the circular strata with NA model columns.
    m <- merge(as.data.table(e), det, by = "ccc_group", all.y = TRUE)
    m[, gene := g]
    m[, circular_by_construction := ccc_group %in% circ][]
}), fill = TRUE)
stopifnot(!any(profile$circular_by_construction & !is.na(profile$emmean)))
write_tsv_safe(profile, file.path(opt$outdir, "ras_celltype_profile.tsv"))

## Rank each cell type per gene -- the readable "who makes what" table.
## Circular strata carry no emmean, so they rank NA and cannot take a slot.
profile[, rank_in_gene := NA_integer_]
profile[!is.na(emmean),
        rank_in_gene := frank(-emmean, ties.method = "min"), by = gene]
top <- profile[!is.na(rank_in_gene) & rank_in_gene <= 3][
    order(gene, rank_in_gene), .(gene, ccc_group, emmean, detect, rank_in_gene)]
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
    ## This site previously negated `estimate` and left BOTH t.ratio and the
    ## `contrast` label describing the opposite comparison -- all 42 rows of
    ## agt_source_posthoc.tsv read backwards from their own estimate.
    ct <- flip_contrast(ct, ref)
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
    "Circular units (stratum defined by the gene being scored; excluded before",
    "within-dataset standardization, retained with emmean = NA and the flag",
    "circular_by_construction in ras_celltype_profile.tsv):",
    paste0("  ", names(CIRCULAR_UNITS), ": ",
           vapply(CIRCULAR_UNITS, paste, character(1), collapse = ", ")),
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
