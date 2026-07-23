## How does AGT relate to the other ligands reaching pericytes?
##
## The frozen NicheNet run places AGT at rank 11 of the ligands predicted to
## regulate the pericyte target program (aupr_corrected 0.103), behind TGFB2
## (0.208) and TGFB1 (0.203). That is a point estimate off a single gene set, and
## a bare ordering invites over-reading, so this script asks three questions with
## uncertainty attached:
##
##   (A) RANK STABILITY. Bootstrap the target gene set and re-score, giving a
##       confidence interval on AGT's AUPR and rank. Answers "is AGT reliably
##       mid-tier, or could it be top-5 / bottom-30 depending on which targets
##       happen to be in the set?"
##   (B) CO-EXPRESSION. Across donors, within each sender cell type, does AGT
##       track TGF-beta -- i.e. are these the same axis or independent inputs?
##       Uses donor pseudobulk, partialling out depth.
##   (C) TARGET CONVERGENCE. Do AGT and TGF-beta act on overlapping target genes?
##       Jaccard plus a hypergeometric test on the NicheNet ligand-target links.
##
## Prior-network caveat, stated in the output: NicheNet regulatory potential is
## derived from a curated prior, not estimated per donor. (A) and (C) inherit that
## prior; only (B) is a measurement from this dataset.

suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(optparse)
    library(dplyr)
    library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--activities", type = "character"),
    make_option("--links", type = "character"),
    make_option("--pseudobulk", type = "character"),
    make_option("--priors", type = "character", default = NA_character_),
    make_option("--frac-file", type = "character", default = NA_character_,
                dest = "frac_file"),
    make_option("--outdir", type = "character"),
    make_option("--receiver", type = "character", default = "Pericytes"),
    make_option("--nboot", type = "integer", default = 500L),
    make_option("--seed", type = "integer", default = 13L),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

write_tsv_safe <- function(x, file)
    write.table(as.data.frame(x, check.names = FALSE), file = file, sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)

## ---------------------------------------------- (A) rank with uncertainty ----
act <- fread(opt$activities)
setorder(act, -aupr_corrected)
act[, rank := .I]
agt_rank <- act[test_ligand == "AGT", rank]
message("AGT point-estimate rank: ", ifelse(length(agt_rank), agt_rank, NA))

boot_done <- FALSE
if (!is.na(opt$priors) && dir.exists(opt$priors) && !is.na(opt$frac_file) &&
    file.exists(opt$frac_file) && requireNamespace("nichenetr", quietly = TRUE)) {
    suppressPackageStartupMessages(library(nichenetr))
    ltm <- readRDS(file.path(opt$priors, "ligand_target_matrix_nsga2r_final.rds"))
    lrn <- readRDS(file.path(opt$priors, "lr_network_human_21122021.rds")) |>
        distinct(from, to)
    frac <- fread(opt$frac_file, data.table = FALSE)
    rownames(frac) <- frac[[1]]; frac[[1]] <- NULL
    thr <- 0.10
    expressed_in <- function(g) if (!g %in% colnames(frac)) character(0) else
        rownames(frac)[frac[[g]] >= thr]
    recv <- expressed_in(opt$receiver)
    background <- intersect(recv, rownames(ltm))
    senders <- setdiff(colnames(frac), opt$receiver)
    expr_send <- unique(unlist(lapply(senders, expressed_in)))
    pot <- lrn |> filter(from %in% intersect(unique(lrn$from), expr_send),
                         to %in% intersect(unique(lrn$to), recv)) |>
        pull(from) |> unique()
    pot <- intersect(pot, colnames(ltm))

    ## Same target program as the frozen run (union of its categories).
    TARGET_PROGRAM <- c(
        "COL1A1","COL1A2","COL3A1","COL4A1","COL5A1","FN1","SPARC","LOX","POSTN",
        "MMP2","MMP3","MMP9","MMP14","TIMP1","TIMP2","LUM","DCN","BGN","FBN1",
        "ACTA2","TAGLN","MYH11","CNN1","MYL9","TPM1","TPM2",
        "ADAMTS1","THBS1","SERPINE1","CCN2","CCN1","PDGFRB","ITGA5","ITGB1",
        "IL6","CXCL1","CXCL2","CXCL8","CXCL10","CCL2","ICAM1","VCAM1","NFKBIA",
        "SOD2","CCL20","MKI67","PCNA","TOP2A","CCND1","CCNB1","BIRC5")
    geneset <- intersect(TARGET_PROGRAM, background)
    message("Bootstrapping ", opt$nboot, " target-set resamples over ",
            length(geneset), " genes")

    boot <- vapply(seq_len(opt$nboot), function(i) {
        gs <- unique(sample(geneset, length(geneset), replace = TRUE))
        if (length(gs) < 5) return(c(NA_real_, NA_real_))
        a <- suppressMessages(predict_ligand_activities(
            geneset = gs, background_expressed_genes = background,
            ligand_target_matrix = ltm, potential_ligands = pot))
        a <- a[order(-a$aupr_corrected), ]
        r <- which(a$test_ligand == "AGT")[1]
        v <- a$aupr_corrected[a$test_ligand == "AGT"][1]
        c(if (length(r)) r else NA_real_, if (length(v)) v else NA_real_)
    }, numeric(2))

    rk <- boot[1, ]; au <- boot[2, ]
    res <- data.table(
        ligand = "AGT", n_boot = sum(is.finite(rk)),
        rank_point = agt_rank,
        rank_median = median(rk, na.rm = TRUE),
        rank_lo = quantile(rk, 0.025, na.rm = TRUE),
        rank_hi = quantile(rk, 0.975, na.rm = TRUE),
        aupr_point = act[test_ligand == "AGT", aupr_corrected],
        aupr_lo = quantile(au, 0.025, na.rm = TRUE),
        aupr_hi = quantile(au, 0.975, na.rm = TRUE))
    write_tsv_safe(res, file.path(opt$outdir, "agt_ligand_rank_bootstrap.tsv"))
    message("AGT bootstrap rank 95% CI: ", res$rank_lo, " - ", res$rank_hi)
    boot_done <- TRUE
} else {
    message("priors/fraction table unavailable; skipping rank bootstrap")
}

## ------------------------------------------------------- (B) co-expression ----
## Donor-level, within sender cell type. If AGT and TGF-beta were the same axis,
## donors high in one would be high in the other.
pb <- fread(opt$pseudobulk)
pb <- pb[n_cells >= opt$min_cells]
partial_cor <- function(x, y, z) {
    ok <- is.finite(x) & is.finite(y) & is.finite(z)
    if (sum(ok) < 6) return(c(NA_real_, NA_real_, sum(ok)))
    rx <- residuals(lm(x[ok] ~ z[ok])); ry <- residuals(lm(y[ok] ~ z[ok]))
    ct <- suppressWarnings(cor.test(rx, ry, method = "spearman"))
    c(unname(ct$estimate), ct$p.value, sum(ok))
}
partners <- intersect(c("TGFB1__expr", "TGFB2__expr", "TGFB3__expr",
                        "CCN2__expr", "PDGFB__expr"), names(pb))
coex <- rbindlist(lapply(unique(pb$ccc_group), function(g) {
    d <- pb[ccc_group == g]
    if (nrow(d) < 6 || !"AGT__expr" %in% names(d)) return(NULL)
    rbindlist(lapply(partners, function(p) {
        v <- partial_cor(d$AGT__expr, d[[p]], d$mean_log10_total_counts)
        data.table(ccc_group = g, partner = sub("__expr$", "", p),
                   partial_rho = v[1], p_value = v[2], n_donors = v[3])
    }), fill = TRUE)
}), fill = TRUE)
if (nrow(coex)) {
    coex[, p_BH := p.adjust(p_value, method = "BH")]
    setorder(coex, p_BH)
    write_tsv_safe(coex, file.path(opt$outdir, "agt_ligand_coexpression.tsv"))
    message("Top AGT/partner co-expression:")
    print(head(coex, 10))
}

## -------------------------------------------------- (C) target convergence ----
if (file.exists(opt$links)) {
    lk <- fread(opt$links)
    lig_col <- intersect(c("ligand", "from"), names(lk))[1]
    tgt_col <- intersect(c("target", "to"), names(lk))[1]
    tg <- function(l) unique(lk[[tgt_col]][lk[[lig_col]] == l])
    universe <- length(unique(lk[[tgt_col]]))
    agt_t <- tg("AGT")
    conv <- rbindlist(lapply(c("TGFB1", "TGFB2", "TGFB3", "CCN2"), function(l) {
        o <- tg(l)
        if (!length(agt_t) || !length(o)) return(NULL)
        inter <- length(intersect(agt_t, o))
        jac <- inter / length(union(agt_t, o))
        ph <- phyper(inter - 1, length(o), universe - length(o),
                     length(agt_t), lower.tail = FALSE)
        data.table(ligand = l, n_agt_targets = length(agt_t),
                   n_other_targets = length(o), n_shared = inter,
                   jaccard = jac, hyper_p = ph, universe = universe)
    }), fill = TRUE)
    if (nrow(conv)) {
        conv[, hyper_p_BH := p.adjust(hyper_p, method = "BH")]
        write_tsv_safe(conv, file.path(opt$outdir, "agt_target_overlap.tsv"))
        message("AGT vs TGF-beta target convergence:")
        print(conv)
    }
}

readme <- c(
    "AGT versus other pericyte ligands -- generated summary",
    sprintf("AGT point-estimate rank in the frozen NicheNet run: %s",
            ifelse(length(agt_rank), agt_rank, "NA")),
    sprintf("Rank bootstrap run: %s", boot_done),
    "",
    "NOTE: (A) rank and (C) target convergence inherit the NicheNet prior network",
    "and are hypothesis-generating. (B) co-expression is measured in this dataset.",
    "",
    if (exists("conv") && is.data.frame(conv) && nrow(conv))
        paste(utils::capture.output(print(conv)), collapse = "\n") else "",
    "",
    if (nrow(coex))
        paste(utils::capture.output(print(head(coex, 10))), collapse = "\n") else "")
writeLines(readme, file.path(opt$outdir, "agt_vs_ligands_README.txt"))

sessioninfo::session_info()
