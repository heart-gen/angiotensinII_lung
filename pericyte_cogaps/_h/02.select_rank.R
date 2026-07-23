## De-novo CoGAPS rank (nPatterns) selection + reproducibility.
##
## nP was previously fixed a priori to the curated program count. Here we let the
## data choose it, using two standard NMF rank-selection signals computed across a
## sweep of nP x random seeds (fits produced by step_1.sh):
##
##   (A) Cross-seed pattern robustness. At each nP, take the canonical (seed-13)
##       run as reference and greedily 1:1-match every other seed's patterns to it
##       by |Pearson r| over shared genes. A pattern that reappears run-to-run has
##       matched r ~ 1. We summarise each nP by:
##         mean_r  = mean over patterns of (mean matched r across seeds)
##         min_r   = min  over patterns of (mean matched r across seeds)  <- the
##                   weakest program; the first to fall as nP is pushed too high.
##
##   (B) Reconstruction error. meanChiSq per fit (from cogaps_meta_*.tsv),
##       averaged over seeds at each nP; lower = better fit, read for its elbow.
##
## Recommendation rule (a suggestion, not a hard gate): the largest nP whose
## min_r stays >= --stability-threshold. Rationale: keep adding patterns while
## every program is still reproducible; stop once the least-stable one degrades.
## The curated-program overlap (03.cogaps_validate.R) then *validates* the chosen
## nP rather than choosing it.
##
## Outputs (to --outdir):
##   cogaps_seed_stability.tsv       per (np, ref_pattern, seed) matched r
##   cogaps_nP_selection.tsv         per np: mean_r, min_r, meanChiSq, n_seeds
##   figures/cogaps_nP_selection.{pdf,png}   robustness + chiSq vs nP, rec marked
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(patchwork)
})
.libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
            .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INDIR     <- parse_arg("--indir", "../_m")
OUTDIR    <- parse_arg("--outdir", "../_m")
NP_SWEEP  <- as.integer(strsplit(parse_arg("--np-sweep", "4,5,6,7,8,9,10"), ",")[[1]])
SEEDS     <- strsplit(parse_arg("--seed-tags", "1,42,2024"), ",")[[1]]  # non-canonical
THRESH    <- as.numeric(parse_arg("--stability-threshold", "0.8"))
fig_dir   <- file.path(OUTDIR, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
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

## ---- (A) cross-seed robustness at each nP ----------------------------------
seed_rows <- rbindlist(lapply(NP_SWEEP, function(np) {
    ref <- load_loadings(np, "")
    if (is.null(ref)) return(NULL)
    rbindlist(lapply(SEEDS, function(tg) {
        cand <- load_loadings(np, paste0("_seed", tg))
        if (is.null(cand)) return(NULL)
        m <- best_match(ref, cand); m[, `:=`(np = np, seed = tg)]; m[]
    }), fill = TRUE)
}), fill = TRUE)

if (!nrow(seed_rows))
    stop("No matched seed runs found under ", INDIR,
         " -- run step_1.sh (nP x seed array) first.")
fwrite(seed_rows, file.path(OUTDIR, "cogaps_seed_stability.tsv"), sep = "\t")

## per (np, pattern): mean matched r across seeds; then per np: mean/min over patterns
per_pat <- seed_rows[, .(pat_mean_r = mean(r, na.rm = TRUE)),
                     by = .(np, ref_pattern)]
robust  <- per_pat[, .(mean_r = mean(pat_mean_r), min_r = min(pat_mean_r)), by = np]

## ---- (B) reconstruction error at each nP -----------------------------------
## Distributed CoGAPS does NOT populate meanChiSq on the aggregated result (it is
## 0; chi-sq is computed per subset, not on the consensus), so we compute the
## reconstruction error directly from the saved factors: D ~= A %*% t(P), with
## A = feature loadings (genes x patterns) and P = sampleFactors (cells x
## patterns). We report the fraction of total sum-of-squares left unexplained
## (lower = better fit), which is comparable across nP.
load_patterns <- function(np) {
    f <- file.path(INDIR, sprintf("patterns_np%d.tsv.gz", np))
    if (!file.exists(f)) { message("missing: ", f); return(NULL) }
    d <- as.data.frame(fread(f)); rownames(d) <- d[[1]]; d[[1]] <- NULL
    as.matrix(d)
}
recon <- data.table(np = NP_SWEEP, frac_unexplained = NA_real_)
mtx_f <- file.path(INDIR, "cogaps_input.mtx")
if (file.exists(mtx_f)) {
    suppressPackageStartupMessages(library(Matrix))
    D <- as.matrix(readMM(mtx_f))
    rownames(D) <- readLines(file.path(INDIR, "cogaps_genes.tsv"))
    colnames(D) <- readLines(file.path(INDIR, "cogaps_barcodes.tsv"))
    for (k in NP_SWEEP) {
        A <- load_loadings(k, ""); P <- load_patterns(k)
        if (is.null(A) || is.null(P)) next
        g   <- intersect(rownames(D), rownames(A))
        cel <- intersect(colnames(D), rownames(P))
        Dsub <- D[g, cel, drop = FALSE]
        Dhat <- A[g, , drop = FALSE] %*% t(P[cel, , drop = FALSE])
        recon[np == k, frac_unexplained := sum((Dsub - Dhat)^2) / sum(Dsub^2)]
    }
    rm(D); gc()
} else message("cogaps_input.mtx not found; reconstruction error will be NA")

n_fit_dt <- seed_rows[, .(n_seeds = uniqueN(seed)), by = np]
sel <- Reduce(function(a, b) merge(a, b, by = "np", all = TRUE),
              list(robust, recon, n_fit_dt))
setorder(sel, np)
fwrite(sel, file.path(OUTDIR, "cogaps_nP_selection.tsv"), sep = "\t")

## recommendation: largest nP whose weakest program is still reproducible
ok  <- sel[!is.na(min_r) & min_r >= THRESH]
rec <- if (nrow(ok)) max(ok$np) else sel[which.max(min_r), np]

cat("\n== De-novo nP selection ==\n")
print(sel[, .(np, mean_r = round(mean_r, 3), min_r = round(min_r, 3),
              frac_unexplained = round(frac_unexplained, 4), n_seeds)])
cat(sprintf("\nRecommended nP = %d (largest nP with min_r >= %.2f).\n", rec, THRESH))
cat("Cross-check against the reconstruction-error elbow and 03.cogaps_validate.R overlap.\n")

## ---- selection figure (manuscript style: no in-panel titles, direct labels) -
theme_sel <- theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = 10, colour = "grey20"),
          axis.text  = element_text(size = 9, colour = "grey25"),
          plot.tag   = element_text(face = "bold", size = 12),
          legend.position = "none",
          plot.margin = margin(6, 12, 6, 6))

## A -- cross-seed pattern robustness; lines direct-labelled instead of a legend
rob     <- melt(sel[, .(np, mean_r, min_r)], id.vars = "np")
labmap  <- c(mean_r = "mean over\nprograms", min_r = "weakest\nprogram")
lab_pos <- rob[np == max(np)]
pr <- ggplot(rob, aes(np, value, colour = variable)) +
    geom_hline(yintercept = THRESH, linetype = 3, colour = "grey60") +
    annotate("text", x = min(NP_SWEEP), y = THRESH, vjust = -0.6, hjust = 0,
             label = sprintf("reproducibility floor (%.1f)", THRESH),
             size = 2.9, colour = "grey55") +
    geom_vline(xintercept = rec, linetype = 2, colour = "#B2182B") +
    annotate("text", x = rec, y = 0.99, hjust = 1.08, vjust = 0, size = 3.1,
             colour = "#B2182B", label = sprintf("selected nP = %d", rec)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.2) +
    geom_text(data = lab_pos, hjust = 0, nudge_x = 0.15, size = 3, lineheight = 0.9,
              aes(label = labmap[as.character(variable)])) +
    scale_colour_manual(values = c(mean_r = "#2166AC", min_r = "#E08214")) +
    scale_x_continuous(breaks = NP_SWEEP, expand = expansion(mult = c(0.03, 0.24))) +
    coord_cartesian(ylim = c(0.4, 1.0)) +
    labs(x = "number of patterns (nP)", y = "cross-seed matched |r|", tag = "A")

## B -- reconstruction error (fraction of variance left unexplained)
pc <- ggplot(sel, aes(np, frac_unexplained)) +
    geom_vline(xintercept = rec, linetype = 2, colour = "#B2182B") +
    geom_line(linewidth = 0.7, colour = "grey35") +
    geom_point(size = 2.2, colour = "grey20") +
    scale_x_continuous(breaks = NP_SWEEP) +
    labs(x = "number of patterns (nP)", y = "fraction of variance unexplained",
         tag = "B")

pp <- (pr | pc) & theme_sel
save_gg(file.path(fig_dir, "cogaps_nP_selection"), pp, 9, 3.6)

cat("\n---- sessionInfo ----\n")
sessioninfo::session_info()
