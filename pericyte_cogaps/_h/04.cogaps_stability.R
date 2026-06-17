## CoGAPS reproducibility: stability across random seeds and across nPatterns.
##
## Reviewers of NMF/CoGAPS results reasonably ask whether the learned patterns are
## a stable property of the data or a single-run artifact. We quantify two things:
##
##   (A) Seed stability (fixed nP): re-run CoGAPS at the chosen nP with several
##       random seeds (step_4.sh writes feature_loadings_np5_seed<NN>.tsv.gz), then
##       greedily match each reference (seed-13) pattern to its best-correlated
##       pattern in every other seed (Pearson over shared genes). The matched
##       correlation is the stability of that pattern; a pattern that reappears
##       run-to-run has matched r ~ 1.
##
##   (B) nP correspondence: for each reference (nP=5) pattern, find its best match
##       among the nP=7 and nP=9 patterns. High correspondence means the coarser
##       solution's programs persist as the resolution increases (they split, not
##       vanish) -- i.e. nP=5 is not over-/under-fit out of existence.
##
## Outputs (to --outdir):
##   cogaps_seed_stability.tsv          per ref-pattern x seed matched correlation
##   cogaps_seed_stability_summary.tsv  per ref-pattern mean/min matched r over seeds
##   cogaps_np_correspondence.tsv       per ref(nP=5) pattern best match in nP=7,9
##   figures/cogaps_seed_stability.{pdf,png}
##   figures/cogaps_np_correspondence.{pdf,png}
suppressPackageStartupMessages({
    library(data.table); library(ggplot2)
})
.libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
            .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INDIR   <- parse_arg("--indir", "../_m")
REF_NP  <- as.integer(parse_arg("--ref-np", "5"))
SEEDS   <- strsplit(parse_arg("--seed-tags", "1,42,2024"), ",")[[1]]  # tags of re-runs
NP_SET  <- as.integer(strsplit(parse_arg("--np-set", "7,9"), ",")[[1]])
OUTDIR  <- parse_arg("--outdir", "../_m")
fig_dir <- file.path(OUTDIR, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
save_gg <- function(fn, p, w, h)
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), p, width = w, height = h)

load_loadings <- function(np, tag = "") {
    f <- file.path(INDIR, sprintf("feature_loadings_np%d%s.tsv.gz", np, tag))
    if (!file.exists(f)) { message("missing: ", f); return(NULL) }
    d <- as.data.frame(fread(f)); rownames(d) <- d[[1]]; d[[1]] <- NULL
    as.matrix(d)
}

## greedy 1:1 match of reference columns to candidate columns by |Pearson r|
best_match <- function(ref, cand) {
    g <- intersect(rownames(ref), rownames(cand))
    ref <- ref[g, , drop = FALSE]; cand <- cand[g, , drop = FALSE]
    cm <- abs(cor(ref, cand))                       # ref-patterns x cand-patterns
    out <- data.table(ref_pattern = colnames(ref), match_pattern = NA_character_,
                      r = NA_real_)
    avail <- colnames(cand)
    ## assign in order of strongest available correlation (greedy, 1:1)
    pairs <- which(cm == cm, arr.ind = TRUE)
    ord <- order(-cm[pairs])
    used_ref <- character(0)
    for (k in ord) {
        i <- pairs[k, 1]; j <- pairs[k, 2]
        rp <- colnames(ref)[i]; cp <- colnames(cand)[j]
        if (rp %in% used_ref || !(cp %in% avail)) next
        out[ref_pattern == rp, `:=`(match_pattern = cp, r = cm[i, j])]
        used_ref <- c(used_ref, rp); avail <- setdiff(avail, cp)
    }
    out[]
}

ref <- load_loadings(REF_NP, "")
stopifnot(!is.null(ref))

## ---- (A) seed stability -----------------------------------------------------
seed_rows <- rbindlist(lapply(SEEDS, function(tg) {
    cand <- load_loadings(REF_NP, paste0("_seed", tg))
    if (is.null(cand)) return(NULL)
    m <- best_match(ref, cand); m[, seed := tg]; m[]
}), fill = TRUE)

if (nrow(seed_rows)) {
    fwrite(seed_rows, file.path(OUTDIR, "cogaps_seed_stability.tsv"), sep = "\t")
    summ <- seed_rows[, .(n_seeds = .N, mean_r = mean(r, na.rm = TRUE),
                          min_r = min(r, na.rm = TRUE)), by = ref_pattern]
    fwrite(summ, file.path(OUTDIR, "cogaps_seed_stability_summary.tsv"), sep = "\t")
    cat("\n== Seed stability (nP=", REF_NP, "), matched |Pearson r| ==\n", sep = "")
    print(summ)
    cat(sprintf("Overall mean matched r across patterns x seeds: %.3f\n",
                mean(seed_rows$r, na.rm = TRUE)))

    pA <- ggplot(seed_rows, aes(ref_pattern, r, group = ref_pattern)) +
        geom_hline(yintercept = c(0.7, 0.9), linetype = c(3, 2), colour = "grey60") +
        geom_jitter(width = 0.12, height = 0, size = 2.4, alpha = 0.8,
                    aes(colour = seed)) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.5, colour = "black") +
        coord_cartesian(ylim = c(0, 1)) +
        labs(x = sprintf("Reference pattern (nP=%d, seed 13)", REF_NP),
             y = "matched |Pearson r| to re-run", colour = "seed",
             title = "CoGAPS pattern stability across random seeds") +
        theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    save_gg(file.path(fig_dir, "cogaps_seed_stability"), pA, 7, 5)
} else {
    cat("No seed re-runs found (run step_4.sh seed sweep first).\n")
}

## ---- (B) nP correspondence --------------------------------------------------
np_rows <- rbindlist(lapply(NP_SET, function(np) {
    cand <- load_loadings(np, "")
    if (is.null(cand)) return(NULL)
    m <- best_match(ref, cand); m[, target_np := np]; m[]
}), fill = TRUE)

if (nrow(np_rows)) {
    fwrite(np_rows, file.path(OUTDIR, "cogaps_np_correspondence.tsv"), sep = "\t")
    cat("\n== nP correspondence (ref nP=", REF_NP, "), matched |Pearson r| ==\n", sep = "")
    print(np_rows[, .(ref_pattern, target_np, match_pattern, r = round(r, 3))])

    pB <- ggplot(np_rows, aes(factor(target_np), ref_pattern, fill = r)) +
        geom_tile() + geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        labs(x = "target nPatterns", y = sprintf("reference pattern (nP=%d)", REF_NP),
             fill = "|r|", title = "Persistence of nP=5 programs at higher resolution") +
        theme_bw(base_size = 12)
    save_gg(file.path(fig_dir, "cogaps_np_correspondence"), pB, 6, 5)
}

cat("\n---- sessionInfo ----\n")
sessioninfo::session_info()
