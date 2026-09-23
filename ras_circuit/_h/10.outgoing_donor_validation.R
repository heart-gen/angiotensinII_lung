## Figure 5F -- is the predicted pericyte -> neighbour signalling visible across donors?
##
## For each target with a NicheNet run (09.nichenet_outgoing.R): donors whose
## pericytes express more of the top-K prioritized ligands should show more of the
## target's response programme IN THE TARGET CELLS. Unit = donor.
##
##   primary      lmer(program_z ~ ligand_z + receiver_depth + sender_depth + (1|study))
##   _fixed       lm(program_z ~ ligand_z + receiver_depth + sender_depth + dataset)
##   _disease_adj primary + disease_group as a NUISANCE covariate (not interpreted;
##                this repository makes no disease claims)
##   _healthy     primary on Healthy donors only
##
## Null: 1,000 composites of the same size, each ligand replaced by a random gene
## from its detection-matched pool of pericyte-expressed NON-ligand genes
## (09 wrote the pools). p_emp is referred to the null's own centre. A ligand
## composite and a programme are scored in different cells, but ambient RNA can
## still couple them across donors -- which is exactly what the matched null
## prices in.
##
## Pseudobulk convention: log1p(mean(expm1(x))) per donor x cell type
## (basement_membrane/_h/02.niche_pseudobulk.py). The two AT2 strata are
## recombined cell-weighted on the linear scale.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4); library(lmerTest)
})
source("../_h/_stats_common.R")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", default = "./outgoing_pseudobulk.tsv.gz"),
    make_option("--design", default = "./nichenet_outgoing/outgoing_validation_design.tsv"),
    make_option("--outdir", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--n-null", type = "integer", default = 1000L, dest = "n_null"),
    make_option("--seed", type = "integer", default = 13L)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

pb <- read_req(opt$pseudobulk)
des <- read_req(opt$design)
pb <- pb[n_cells >= opt$min_cells]

## recombine AT2 strata (linear scale, cell-weighted)
ecols <- grep("__expr$", names(pb), value = TRUE)
at2 <- pb[ccc_group %in% c("AT2_AGTR2det", "AT2_AGTR2undet")]
if (nrow(at2)) {
    lin <- at2[, lapply(.SD, function(v) log1p(sum(expm1(v) * n_cells) / sum(n_cells))),
               by = .(donor_id), .SDcols = ecols]
    meta <- at2[, .(n_cells = sum(n_cells), study = study[1], dataset = dataset[1],
                    disease_group = disease_group[1],
                    mean_log10_total_counts = sum(mean_log10_total_counts * n_cells) / sum(n_cells)),
                by = donor_id]
    pb <- rbind(pb, merge(meta, lin, by = "donor_id")[, ccc_group := "AT2"], fill = TRUE)
}

composite <- function(units, genes) {
    genes <- genes[paste0(genes, "__expr") %in% names(units)]
    if (!length(genes)) return(rep(NA_real_, nrow(units)))
    m <- sapply(genes, function(g) z_within_dataset(units[[paste0(g, "__expr")]], units$dataset))
    rowMeans(as.matrix(m), na.rm = TRUE)
}

fit_arms <- function(d) {
    out <- list()
    f <- program_z ~ ligand_z + receiver_depth + sender_depth + (1 | study)
    m <- try(suppressMessages(lmer(f, data = d)), silent = TRUE)
    if (!inherits(m, "try-error"))
        out$primary <- tidy_row(m, "ligand_z", "donor", "lmer(+1|study)", nrow(d),
                                uniqueN(d$donor_id), spec = "primary")
    if (uniqueN(d$dataset) > 1) {
        m <- lm(program_z ~ ligand_z + receiver_depth + sender_depth + dataset, data = d)
        co <- summary(m)$coefficients["ligand_z", ]
        out$fixed <- data.table(level = "donor", model = "lm(+dataset fixed)", spec = "_fixed",
                                predictor = "ligand_z", estimate = co[1], SE = co[2],
                                z_value = co[3], p_value = co[4], n_units = nrow(d),
                                n_donors = uniqueN(d$donor_id), singular = NA, note = "")
    }
    if (uniqueN(d$disease_group) > 1) {
        m <- try(suppressMessages(lmer(update(f, . ~ . + disease_group), data = d)), silent = TRUE)
        if (!inherits(m, "try-error"))
            out$dz <- tidy_row(m, "ligand_z", "donor", "lmer(+1|study)+disease_group",
                               nrow(d), uniqueN(d$donor_id), spec = "_disease_adj")
    }
    h <- d[disease_group == "Healthy"]
    if (nrow(h) >= 15) {
        m <- try(suppressMessages(lmer(f, data = h)), silent = TRUE)
        if (!inherits(m, "try-error"))
            out$h <- tidy_row(m, "ligand_z", "donor", "lmer(+1|study)", nrow(h),
                              uniqueN(h$donor_id), spec = "_healthy")
    }
    rbindlist(out, fill = TRUE)
}

res <- list(); nulls <- list()
for (tgt in unique(des$target)) {
    dd <- des[target == tgt]
    ligs <- dd[role == "ligand", gene]
    prog_genes <- dd[role == "program", gene]
    sender <- pb[ccc_group == "Pericytes"]
    recv <- pb[ccc_group == tgt]
    if (!nrow(recv)) { message(tgt, ": no receiver units"); next }
    recv[, program_z := composite(recv, prog_genes)]
    base <- merge(sender[, .(donor_id, study, dataset, disease_group,
                             sender_depth = mean_log10_total_counts)],
                  recv[, .(donor_id, program_z, receiver_depth = mean_log10_total_counts)],
                  by = "donor_id")
    sender[, ligand_z := composite(sender, ligs)]
    d <- merge(base, sender[, .(donor_id, ligand_z)], by = "donor_id")
    d <- d[is.finite(program_z) & is.finite(ligand_z)]
    message(sprintf("%s: %d donors, %d ligands (%s), %d programme genes", tgt,
                    nrow(d), length(ligs), paste(ligs, collapse = ","), length(prog_genes)))
    if (nrow(d) < 15) { message("  < 15 donors; skipped"); next }
    r <- fit_arms(d)
    r[, `:=`(target = tgt, program = dd$program[1], ligands = paste(ligs, collapse = ";"),
             n_studies = uniqueN(d$study))]

    ## null composites
    pools <- lapply(ligs, function(l) dd[role == "null_pool" & matched_to == l, gene])
    ok <- lengths(pools) > 0
    nv <- rep(NA_real_, opt$n_null)
    if (all(ok)) {
        for (b in seq_len(opt$n_null)) {
            pick <- vapply(pools, function(p) p[sample.int(length(p), 1)], "")
            s2 <- copy(sender); s2[, ligand_z := composite(s2, pick)]
            dn <- merge(base, s2[, .(donor_id, ligand_z)], by = "donor_id")
            dn <- dn[is.finite(program_z) & is.finite(ligand_z)]
            m <- try(suppressMessages(lmer(program_z ~ ligand_z + receiver_depth +
                                               sender_depth + (1 | study), data = dn)),
                     silent = TRUE)
            if (!inherits(m, "try-error")) nv[b] <- fixef(m)[["ligand_z"]]
        }
        nulls[[tgt]] <- data.table(target = tgt, draw = seq_len(opt$n_null), estimate = nv)
    } else message("  a ligand has no matched pool; null skipped")
    r[, `:=`(null_mean = mean(nv, na.rm = TRUE), null_sd = sd(nv, na.rm = TRUE),
             n_null = sum(is.finite(nv)))]
    r[spec == "primary", p_emp := emp_p(estimate, nv)]
    res[[tgt]] <- r
}
res <- rbindlist(res, fill = TRUE)
res[spec == "primary", p_BH := p.adjust(p_value, "BH")]
write_tsv_safe(res, file.path(opt$outdir, "outgoing_donor_validation.tsv"))
write_tsv_safe(rbindlist(nulls), file.path(opt$outdir, "outgoing_donor_validation_null.tsv"))
print(res[, .(target, spec, estimate, SE, p_value, n_donors, null_mean, p_emp)])
cat("\nReproducibility information:\n"); print(sessionInfo())
