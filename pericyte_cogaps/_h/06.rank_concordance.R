## Does the CoGAPS story depend on the rank we picked?
##
## step_2 chose nP on cross-seed reproducibility and carried TWO ranks forward:
## nP=8 (main; min_r 0.978) and nP=9 (sensitivity; the largest nP with
## min_r >= 0.80). 02.select_rank.R already writes a np-correspondence table, but
## it is anchored on nP=5 -- the rank that was canonical before the sweep -- so
## nothing on disk compares the two ranks we actually report.
##
## This script does exactly that, and nothing else:
##   1. greedily 1:1-match every nP=8 pattern to an nP=9 pattern by |Pearson r|
##      over the shared HVG loadings (same algorithm as 02.select_rank.R, so the
##      numbers are on the same footing as the cross-seed ones);
##   2. for each pattern at each rank, read its curated-program identity from the
##      cell-level Spearman table written by 03.cogaps_validate.R -- the argmax
##      program and its rho;
##   3. correlate the two 6-program rho PROFILES of a matched pair. This is the
##      number that answers the question: a pair with profile_r ~ 1 stands for the
##      same biology at both ranks, whatever we call it.
##
## Interpretation caveat carried into the legend: argmax over rho is a weak label
## when the winning rho is small. nP=8 Pattern_2 and Pattern_8 correlate NEGATIVELY
## with five of the six panels and are labelled basement_membrane only because that
## is their one positive correlation (rho 0.25 / 0.39); the extra pattern nP=9 buys
## tops out at rho 0.16. top_rho is emitted alongside the name so the figure can
## show that rather than hide it.
##
## Outputs (to --outdir):
##   cogaps_np8_vs_np9_concordance.tsv   one row per nP=8 pattern + any unmatched
##                                       nP=9 pattern (matched = FALSE)
suppressPackageStartupMessages({ library(data.table) })
.libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
            .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INDIR  <- parse_arg("--indir", "../_m")
OUTDIR <- parse_arg("--outdir", "../_m")
NP_A   <- as.integer(parse_arg("--np-main", "8"))
NP_B   <- as.integer(parse_arg("--np-sens", "9"))

## The six curated NVU panels, in the order 03.cogaps_validate.R scores them.
PROGRAMS <- c("vascular_stabilizing", "inflammatory", "synthetic_contractile",
              "activated_migratory", "fibroblast_like", "basement_membrane")

load_loadings <- function(np) {
    f <- file.path(INDIR, sprintf("feature_loadings_np%d.tsv.gz", np))
    if (!file.exists(f)) stop("missing: ", f)
    d <- as.data.frame(fread(f)); rownames(d) <- d[[1]]; d[[1]] <- NULL
    as.matrix(d)
}

## Program-identity profile: pattern x 6 curated programs, cell-level Spearman.
load_profile <- function(np) {
    f <- file.path(INDIR, sprintf("validation_np%d", np), "pattern_score_spearman.tsv.gz")
    if (!file.exists(f)) stop("missing: ", f, " -- run 03.cogaps_validate.R first")
    d <- fread(f)
    d[, program := sub("_score$", "", score)]
    d <- d[program %in% PROGRAMS & is.finite(rho)]
    dcast(d, pattern ~ factor(program, levels = PROGRAMS), value.var = "rho")
}

## greedy 1:1 match of reference columns to candidate columns by |Pearson r|
## (verbatim from 02.select_rank.R -- kept identical on purpose)
best_match <- function(ref, cand) {
    g <- intersect(rownames(ref), rownames(cand))
    ref <- ref[g, , drop = FALSE]; cand <- cand[g, , drop = FALSE]
    cm <- abs(cor(ref, cand))
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
    list(matches = out[], unmatched = avail)
}

La <- load_loadings(NP_A); Lb <- load_loadings(NP_B)
cat(sprintf("loadings: nP=%d %d x %d ; nP=%d %d x %d ; %d shared genes\n",
            NP_A, nrow(La), ncol(La), NP_B, nrow(Lb), ncol(Lb),
            length(intersect(rownames(La), rownames(Lb)))))

mm <- best_match(La, Lb)
res <- mm$matches
setnames(res, c("ref_pattern", "match_pattern", "r"),
         c("pattern_np_main", "pattern_np_sens", "r_loadings"))

pa <- load_profile(NP_A); pb <- load_profile(NP_B)
prof_mat <- function(p) { m <- as.matrix(p[, ..PROGRAMS]); rownames(m) <- p$pattern; m }
Ma <- prof_mat(pa); Mb <- prof_mat(pb)

top_of <- function(M, pat) {
    v <- M[pat, ]
    list(program = names(v)[which.max(v)], rho = max(v))
}

res[, `:=`(
    top_program_main = vapply(pattern_np_main, function(p) top_of(Ma, p)$program, ""),
    top_rho_main     = vapply(pattern_np_main, function(p) top_of(Ma, p)$rho, 0),
    top_program_sens = vapply(pattern_np_sens, function(p) top_of(Mb, p)$program, ""),
    top_rho_sens     = vapply(pattern_np_sens, function(p) top_of(Mb, p)$rho, 0),
    profile_r        = mapply(function(a, b) cor(Ma[a, ], Mb[b, ]),
                              pattern_np_main, pattern_np_sens))]
res[, `:=`(program_agree = top_program_main == top_program_sens,
           matched = TRUE, np_main = NP_A, np_sens = NP_B)]

## The extra pattern(s) the higher rank buys, carried as explicit rows so the
## figure cannot quietly drop them.
if (length(mm$unmatched)) {
    extra <- data.table(
        pattern_np_main = NA_character_, pattern_np_sens = mm$unmatched,
        r_loadings = NA_real_, top_program_main = NA_character_,
        top_rho_main = NA_real_,
        top_program_sens = vapply(mm$unmatched, function(p) top_of(Mb, p)$program, ""),
        top_rho_sens = vapply(mm$unmatched, function(p) top_of(Mb, p)$rho, 0),
        profile_r = NA_real_, program_agree = NA, matched = FALSE,
        np_main = NP_A, np_sens = NP_B)
    res <- rbind(res, extra)
}

fwrite(res, file.path(OUTDIR, sprintf("cogaps_np%d_vs_np%d_concordance.tsv",
                                      NP_A, NP_B)), sep = "\t")

m <- res[matched == TRUE]
cat(sprintf("matched %d/%d patterns; median r_loadings = %.3f (min %.3f); ",
            nrow(m), ncol(La), median(m$r_loadings), min(m$r_loadings)))
cat(sprintf("median profile_r = %.3f; argmax program agrees for %d/%d\n",
            median(m$profile_r), sum(m$program_agree), nrow(m)))
if (length(mm$unmatched))
    cat("unmatched nP=", NP_B, " pattern(s): ", paste(mm$unmatched, collapse = ", "),
        "\n", sep = "")
print(res)

cat("\n---- sessionInfo ----\n")
sessioninfo::session_info()
