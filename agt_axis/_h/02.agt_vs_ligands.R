## How does AGT relate to the other ligands reaching pericytes?
##
## The frozen NicheNet run places AGT at rank 11 of the ligands predicted to
## regulate the pericyte target program (aupr_corrected 0.103), behind TGFB2
## (0.208) and TGFB1 (0.203). That is a point estimate off a single gene set, and
## a bare ordering invites over-reading, so this script asks three questions with
## uncertainty attached:
##
##   (A) RANK STABILITY. Subsample the target gene set at fixed size and
##       re-score, giving a confidence interval on AGT's AUPR and rank. Answers
##       "is AGT reliably mid-tier, or could it be top-5 / bottom-30 depending
##       on which targets happen to be in the set?" The INTERVAL is the result;
##       the median is not a bias-corrected rank (see the scheme note below).
##   (B) CO-EXPRESSION. Across donors, within each sender cell type, does AGT
##       track TGF-beta -- i.e. are these the same axis or independent inputs?
##       Uses donor pseudobulk, partialling out depth.
##   (C) TARGET OVERLAP. Do AGT and TGF-beta act on overlapping target genes?
##       Jaccard on the NicheNet ligand-target links, read against the overlap
##       of every other ligand pair in the same table and against a
##       degree-preserving null. The hypergeometric test this block used to
##       report is retained but flagged: see the note above section (C).
##
## Prior-network caveat, stated in the output: NicheNet regulatory potential is
## derived from a curated prior, not estimated per donor. (A) and (C) inherit
## that prior; only (B) is a measurement from this dataset. (B) and (C) are
## therefore NOT two independent lines of evidence, and must not be written up
## as though they were -- (C) is a consistency check on the prior that produced
## the ranking in (A), not a second observation of the same biology.

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
    make_option("--subsample-frac", type = "double", default = 0.8,
                dest = "subsample_frac"),
    make_option("--seed", type = "integer", default = 13L),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    ## Sender AGT detection floor for the co-expression arm (P2-18). See the long
    ## note at (B). 0.01 keeps every group with genuine AGT expression and drops
    ## the near-zero-detection groups whose correlations are dropout structure.
    make_option("--min-agt-detect", type = "double", default = 0.01,
                dest = "min_agt_detect")
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

    ## RESAMPLING SCHEME. This previously drew
    ##   unique(sample(geneset, length(geneset), replace = TRUE))
    ## which is a multiset bootstrap collapsed back to a set. A target program
    ## is a SET -- a gene counted twice contributes nothing to an AUPR -- so the
    ## duplicates are discarded and each replicate scores a smaller program
    ## than the real one: with 30 genes, 19.1 +/- 1.7 (95% range 16-22). Two
    ## consequences, and they pull in different directions:
    ##
    ##   (i)  the ~36% shrinkage is systematic, so `rank_median` is not a
    ##        bias-corrected version of `rank_point`; the 11 -> 18 gap is part
    ##        instability and part arithmetic, and nothing in the old output
    ##        let a reader tell which;
    ##   (ii) the replicate SIZE also varies, adding a nuisance variance that
    ##        widened the interval on top of shifting it.
    ##
    ## Fixed-size subsampling without replacement removes both: every replicate
    ## scores exactly `m` genes, so the spread is set-membership uncertainty
    ## alone. `size_ref` then scores the SAME m on the real program's top-m
    ## by prior weight -- a size-matched reference that says what rank an
    ## m-gene set costs on its own, so the shrinkage is measurable rather than
    ## silently mixed into the interval.
    m <- max(5L, round(opt$subsample_frac * length(geneset)))
    message(sprintf(paste("Subsampling %d replicates of m = %d of %d target",
                          "genes (fixed size, without replacement)"),
                    opt$nboot, m, length(geneset)))

    score_set <- function(gs) {
        if (length(gs) < 5) return(c(NA_real_, NA_real_))
        a <- suppressMessages(predict_ligand_activities(
            geneset = gs, background_expressed_genes = background,
            ligand_target_matrix = ltm, potential_ligands = pot))
        a <- a[order(-a$aupr_corrected), ]
        r <- which(a$test_ligand == "AGT")[1]
        v <- a$aupr_corrected[a$test_ligand == "AGT"][1]
        c(if (length(r)) r else NA_real_, if (length(v)) v else NA_real_)
    }

    boot <- vapply(seq_len(opt$nboot), function(i)
        score_set(sample(geneset, m, replace = FALSE)), numeric(2))
    rk <- boot[1, ]; au <- boot[2, ]

    ## Per-replicate draws. Without these the bias in a resampling scheme is
    ## not diagnosable after the fact -- the previous run wrote quantiles only.
    write_tsv_safe(data.table(replicate = seq_len(opt$nboot), m = m,
                              rank = rk, aupr_corrected = au),
                   file.path(opt$outdir, "agt_ligand_rank_bootstrap_draws.tsv"))

    ## Size-matched reference: the m highest-weight genes of the real program.
    ## Deterministic, so it isolates the cost of scoring m genes instead of n.
    agt_w <- if ("AGT" %in% colnames(ltm))
        ltm[intersect(geneset, rownames(ltm)), "AGT"] else NULL
    size_ref <- if (!is.null(agt_w) && length(agt_w) >= m)
        score_set(names(sort(agt_w, decreasing = TRUE))[seq_len(m)]) else
            c(NA_real_, NA_real_)

    res <- data.table(
        ligand = "AGT", n_boot = sum(is.finite(rk)),
        n_geneset = length(geneset), m_subsample = m,
        scheme = "fixed-size subsample without replacement",
        rank_point = agt_rank,
        ## Reported for completeness. NOT a corrected point estimate: it is the
        ## centre of an m-gene distribution, and `rank_size_ref` is the
        ## size-matched benchmark it must be read against.
        rank_median = median(rk, na.rm = TRUE),
        rank_size_ref = size_ref[1],
        rank_lo = quantile(rk, 0.025, na.rm = TRUE),
        rank_hi = quantile(rk, 0.975, na.rm = TRUE),
        aupr_point = act[test_ligand == "AGT", aupr_corrected],
        aupr_size_ref = size_ref[2],
        aupr_lo = quantile(au, 0.025, na.rm = TRUE),
        aupr_hi = quantile(au, 0.975, na.rm = TRUE))
    write_tsv_safe(res, file.path(opt$outdir, "agt_ligand_rank_bootstrap.tsv"))
    message("AGT subsample rank 95% CI: ", res$rank_lo, " - ", res$rank_hi,
            "; median ", res$rank_median, "; size-matched reference ",
            res$rank_size_ref)
    boot_done <- TRUE
} else {
    message("priors/fraction table unavailable; skipping rank bootstrap")
}

## ------------------------------------------------------- (B) co-expression ----
## Donor-level, within sender cell type. If AGT and TGF-beta were the same axis,
## donors high in one would be high in the other.
##
## THE EXPRESSION FLOOR (P2-18, added 2026-09-08).
##
## `partial_cor` requires 6 FINITE observations, not 6 donors in which AGT is
## actually detected. A donor pseudobulk of a cell type that essentially never
## expresses AGT is still a finite number -- it is a number near zero whose
## variation is dropout and depth, not transcription. Correlating that against a
## partner gene measures how the two genes' zeros co-occur across donors, which
## is a library-composition statistic wearing a co-expression label.
##
## The consequence was not subtle. Ranked by BH, the top of the old output was:
##
##   Non-classical monocytes  AGT detection 0.00034   |rho| 0.809   BH ~ 0
##   DC2                                    0.00059         0.774        ~ 0
##   Mast cells                             0.00038         0.769        ~ 0
##   Classical monocytes                    0.00009         0.766        ~ 0
##   Vascular smooth muscle                 0.101           0.144         0.144
##
## -- 35 of 42 BH-significant rows came from cell types detecting AGT in under
## 1% of cells, and the file is `setorder`ed by `p_BH`, so those rows were the
## first thing a reader saw. The ONE sender with real AGT expression ranked last.
## `AGT_SUMMARY.md` caught this and refused to report it, but the file itself
## carried no marker, so the refusal lived only in prose.
##
## The floor is applied to the sender's AGT detection, computed cell-weighted
## over exactly the donors that enter the test (not over the whole cohort), so
## the number printed beside a row describes that row's own fit. At the default
## 0.01 seven groups qualify:
##
##   Vascular smooth muscle 0.101, Adventitial fib. 0.045, Peribronchial fib.
##   0.033, Myofibroblasts 0.022, Pericytes 0.018, Alveolar fib. 0.017,
##   AT2_AGTR2det 0.011
##
## and BH tightens from 110 tests to 35, which STRENGTHENS the rows that matter
## rather than weakening them -- the artifact rows were consuming the correction.
##
## Rows below the floor are still WRITTEN, with `tested = FALSE` and a reason, so
## the exclusion is visible in the output instead of being a silent filter. They
## are excluded from the BH family. Two smaller defects close with the same
## change: the lone NA row (Subpleural fibroblasts x PDGFB, 12 donors, no
## variance) no longer inflates the BH denominator -- `p.adjust` takes `n` at
## call time, so an NA had been counted as a test -- and that group sits at
## detection 0.0099, below the floor, so it leaves on its own terms too.
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

## Cell-weighted sender AGT detection over the donors that survive --min-cells.
agt_det <- if ("AGT__detect" %in% names(pb))
    pb[, .(agt_detect_group = sum(AGT__detect * n_cells) / sum(n_cells),
           n_cells_group = sum(n_cells)), by = ccc_group] else NULL

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
    if (!is.null(agt_det)) coex <- merge(coex, agt_det, by = "ccc_group", all.x = TRUE)
    else coex[, `:=`(agt_detect_group = NA_real_, n_cells_group = NA_integer_)]

    ## Spearman P underflows to a literal 0 at these sample sizes, which is not a
    ## representable P and reads as certainty. Flag it and clamp to the double
    ## epsilon so BH sees a number rather than a zero. Every underflow row in the
    ## previous run was an artifact row, so this is bookkeeping, not a rescue.
    coex[, p_underflow := is.finite(p_value) & p_value <= 0]
    coex[p_underflow == TRUE, p_value := .Machine$double.eps]

    coex[, tested := is.finite(p_value) &
             is.finite(agt_detect_group) & agt_detect_group >= opt$min_agt_detect]
    coex[, excluded_reason := fifelse(
        tested, "",
        fifelse(!is.finite(p_value), "no variance in one arm (P not estimable)",
                sprintf("sender AGT detection %.4g < floor %.4g",
                        agt_detect_group, opt$min_agt_detect)))]

    ## BH over the tested rows ONLY. `p.adjust` evaluates `n` at call time, so
    ## passing the full column would keep the excluded rows in the denominator.
    coex[, p_BH := NA_real_]
    coex[tested == TRUE, p_BH := p.adjust(p_value, method = "BH")]
    setorder(coex, -tested, p_BH, na.last = TRUE)
    write_tsv_safe(coex, file.path(opt$outdir, "agt_ligand_coexpression.tsv"))

    ## The floor is a judgement call, so record what other floors would have
    ## done. A reader can then see whether a claim depends on where it was set.
    if (!is.null(agt_det)) {
        sweep <- rbindlist(lapply(c(0.005, 0.01, 0.02, 0.05), function(f) {
            k <- coex[is.finite(p_value) & agt_detect_group >= f]
            if (!nrow(k)) return(data.table(floor = f, n_groups = 0L, n_tests = 0L,
                                            n_BH_sig = 0L, min_BH = NA_real_))
            bh <- p.adjust(k$p_value, method = "BH")
            data.table(floor = f, n_groups = uniqueN(k$ccc_group), n_tests = nrow(k),
                       n_BH_sig = sum(bh < 0.05), min_BH = min(bh))
        }))
        write_tsv_safe(sweep, file.path(opt$outdir, "agt_coexpression_floor_sweep.tsv"))
        message("AGT detection floor sweep:"); print(sweep)
    }

    message(sprintf(paste("Co-expression: %d of %d rows tested (sender AGT",
                          "detection >= %.4g); %d excluded by the floor,",
                          "%d not estimable"),
                    sum(coex$tested), nrow(coex), opt$min_agt_detect,
                    sum(!coex$tested & is.finite(coex$p_value)),
                    sum(!is.finite(coex$p_value))))
    message("Top AGT/partner co-expression (tested rows only):")
    print(head(coex[tested == TRUE], 10))
}

## -------------------------------------------------- (C) target convergence ----
## THE HYPERGEOMETRIC NULL HERE IS WRONG, and it is kept only so the published
## value stays traceable. It assumes each ligand draws its targets uniformly
## from the universe. That universe is the shortlisted link table -- 24 genes
## over 30 ligands -- and the targets are wildly unequal in popularity:
## SERPINE1, CCN2, FN1, MMP2, NFKBIA and THBS1 are each hit by 26-28 of the 30
## ligands. Any two ligands share those almost automatically, and six of the ten
## genes AGT shares with CCN2 are exactly those six. Under the uniform null the
## overlap looks extraordinary (P = 1.7e-3); against a null that preserves how
## often each target is used it is ordinary (P ~ 0.20).
##
## So two calibrated readouts are reported alongside it:
##
##   pair_pctile  -- where this pair's Jaccard falls among ALL ligand pairs in
##                   the same table. Model-free, and the honest framing: the
##                   median pair here already shares Jaccard ~0.58, so "shares
##                   most of its targets" is the norm, not a finding.
##   degree_p     -- permutation P against a degree-preserving null (targets
##                   drawn with probability proportional to their usage across
##                   ligands, at the observed set sizes).
##
## This is (C) in the header's caveat: entirely inherited from the curated
## prior. No donor data enters it, so it cannot be a line of evidence
## independent of the NicheNet run -- only a consistency check on that prior.
if (file.exists(opt$links)) {
    lk <- fread(opt$links)
    lig_col <- intersect(c("ligand", "from"), names(lk))[1]
    tgt_col <- intersect(c("target", "to"), names(lk))[1]
    tg <- function(l) unique(lk[[tgt_col]][lk[[lig_col]] == l])
    all_targets <- unique(lk[[tgt_col]])
    universe <- length(all_targets)
    agt_t <- tg("AGT")

    ## Background: the Jaccard of every ligand pair in this table.
    all_ligs <- unique(lk[[lig_col]])
    tsets <- setNames(lapply(all_ligs, tg), all_ligs)
    jac_of <- function(a, b)
        length(intersect(a, b)) / max(1L, length(union(a, b)))
    pair_bg <- if (length(all_ligs) >= 3) combn(all_ligs, 2, function(p)
        jac_of(tsets[[p[1]]], tsets[[p[2]]])) else numeric(0)

    ## Degree-preserving null: target usage across ligands sets the draw weights.
    deg <- as.numeric(table(factor(lk[[tgt_col]], levels = all_targets)))
    NPERM <- 20000L
    degree_p <- function(na, nb, obs) {
        if (!na || !nb) return(NA_real_)
        hits <- sum(vapply(seq_len(NPERM), function(i) {
            a <- sample(all_targets, na, prob = deg)
            b <- sample(all_targets, nb, prob = deg)
            length(intersect(a, b)) >= obs
        }, logical(1)))
        (hits + 1) / (NPERM + 1)
    }

    conv <- rbindlist(lapply(c("TGFB1", "TGFB2", "TGFB3", "CCN2"), function(l) {
        o <- tg(l)
        if (!length(agt_t) || !length(o)) return(NULL)
        inter <- length(intersect(agt_t, o))
        jac <- inter / length(union(agt_t, o))
        ph <- phyper(inter - 1, length(o), universe - length(o),
                     length(agt_t), lower.tail = FALSE)
        data.table(ligand = l, n_agt_targets = length(agt_t),
                   n_other_targets = length(o), n_shared = inter,
                   jaccard = jac, universe = universe,
                   n_ligand_pairs = length(pair_bg),
                   jaccard_median_all_pairs = median(pair_bg),
                   pair_pctile = 100 * mean(pair_bg < jac),
                   degree_p = degree_p(length(agt_t), length(o), inter),
                   hyper_p_MISCALIBRATED = ph)
    }), fill = TRUE)
    if (nrow(conv)) {
        conv[, degree_p_BH := p.adjust(degree_p, method = "BH")]
        conv[, hyper_p_MISCALIBRATED_BH := p.adjust(hyper_p_MISCALIBRATED,
                                                    method = "BH")]
        write_tsv_safe(conv, file.path(opt$outdir, "agt_target_overlap.tsv"))
        message("AGT vs TGF-beta target overlap (prior-derived, not independent):")
        print(conv)
    }
}

readme <- c(
    "AGT versus other pericyte ligands -- generated summary",
    sprintf("AGT point-estimate rank in the frozen NicheNet run: %s",
            ifelse(length(agt_rank), agt_rank, "NA")),
    sprintf("Rank bootstrap run: %s", boot_done),
    "",
    "NOTE: (A) rank and (C) target overlap inherit the NicheNet prior network and",
    "are hypothesis-generating. (B) co-expression is the only block measured in",
    "this dataset, so (B) and (C) are NOT two independent lines of evidence.",
    "In (C), read pair_pctile and degree_p; hyper_p_MISCALIBRATED assumes uniform",
    "target draws, which this 24-gene shortlist badly violates.",
    "",
    if (exists("conv") && is.data.frame(conv) && nrow(conv))
        paste(utils::capture.output(print(conv)), collapse = "\n") else "",
    "",
    if (nrow(coex))
        paste(utils::capture.output(print(head(coex, 10))), collapse = "\n") else "")
writeLines(readme, file.path(opt$outdir, "agt_vs_ligands_README.txt"))

sessioninfo::session_info()
