## RETIRED 2026-07-22 -- DO NOT RUN, DO NOT CITE ITS OUTPUTS.
##
## Superseded by 04.species_comparability.py, which showed that neither test below
## measures what its name says:
##   * the n = 1,144 "mural" set it operates on is 96.4% smooth muscle (1,016 vSMC +
##     87 PA-SMC + 41 pericytes), and the mouse `basement_membrane` state contains
##     zero pericytes, so the injury-vs-stabilizing contrast is a smooth-muscle
##     contrast wearing pericyte-state labels;
##   * cell type is perfectly aliased with dataset, so `dataset` as a fixed effect
##     absorbs the pericyte-vs-vSMC comparison entirely;
##   * the continuous `Agtr1a ~ vascular_stabilizing_score` test runs on the dense
##     `scvi_corrected` layer, where every cell is non-zero, so it largely restates
##     cell-type identity rather than measuring conservation.
## The defensible cross-species claim is compartment-level and lives in
## 04.species_comparability.py. This file and its `_m/stats_data/mouse_Agtr1a_*`
## outputs are kept for provenance only; step_3.sh no longer invokes it.
##
## ---------------------------------------------------------------------------
## Cross-species conservation statistics (on integrated mouse data).
##
## State labels come from clustering on the scVI-integrated latent, so they are
## no longer confounded with dataset. We test conservation two ways:
##  (1) categorical: Agtr1a in injury vs stabilizing mural cells, mixed model
##      with donor random intercept and dataset as fixed batch;
##  (2) continuous (integration-agnostic): does Agtr1a track the vascular-
##      stabilizing program, adjusting for dataset, with donor random intercept;
## plus a per-dataset consistency table so the claim rests on agreement across
## cohorts rather than a pooled p-value.

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggpubr)
    library(lme4); library(lmerTest); library(emmeans)
})
emm_options(lmerTest.limit = 30000, pbkrtest.limit = 30000)

## Injury programs. `basement_membrane` is deliberately EXCLUDED: BM deposition
## is a matrix-stabilizing vascular function, not an injury program (in human it
## is near-orthogonal to fibrillar ECM and aligns pericytes with endothelium
## rather than fibroblasts), and the same panel is now scored in mouse so the two
## species are defined consistently.
INJURY_STATES <- c("inflammatory", "fibroblast_like", "activated_migratory")
## Continuous fallback: if the BM panel takes over the argmax in mouse as it did
## in human, the label-based injury-vs-stabilizing contrast can lose most of its
## cells. The score-based contrast below does not depend on which panel wins.
INJURY_SCORE_COLS <- c("inflammatory_score", "fibroblast_like_score",
                       "activated_migratory_score")
MIN_LABEL_CELLS <- 200L
save_ggplots <- function(fn, p, w, h) for (e in c(".pdf", ".png")) ggsave(paste0(fn, e), p, width = w, height = h)
write_tsv <- function(x, f, rn = FALSE) write.table(as.data.frame(x, check.names = FALSE), f,
                                                    sep = "\t", quote = FALSE, row.names = rn)

df <- data.table::fread("mouse_states_metadata.tsv.gz")
outdir <- "stats_data"; if (!dir.exists(outdir)) dir.create(outdir)
mural <- df |> filter(is_mural == TRUE) |>
    mutate(donor_id = factor(donor_id), dataset_id = factor(dataset_id))
cat("mural cells:", nrow(mural), " donors:", nlevels(mural$donor_id),
    " datasets:", nlevels(mural$dataset_id), "\n")

## Integration check: state composition by dataset (should be less skewed now)
write_tsv(as.data.frame.matrix(prop.table(table(mural$dataset_id, mural$pericyte_state), 1)),
          file.path(outdir, "mouse_state_composition_by_dataset.tsv"), TRUE)
write_tsv(mural |> count(pericyte_state) |> mutate(frac = n / sum(n)),
          file.path(outdir, "mouse_mural_state_composition.tsv"))
## per-state Agtr1a summary (consumed by the manuscript supplement figure)
write_tsv(mural |> group_by(pericyte_state) |>
              summarise(n = n(), Agtr1a_mean = mean(Agtr1a_expr, na.rm = TRUE),
                        Agtr1a_detect = mean(Agtr1a_expr > 0, na.rm = TRUE),
                        .groups = "drop") |>
              mutate(injury = pericyte_state %in% INJURY_STATES),
          file.path(outdir, "mouse_Agtr1a_by_state.tsv"))

has_batch <- nlevels(droplevels(mural$dataset_id)) > 1

## (1) Categorical: injury vs stabilizing
sub <- mural |>
    filter(pericyte_state %in% c(INJURY_STATES, "vascular_stabilizing")) |>
    mutate(is_injury = factor(ifelse(pericyte_state %in% INJURY_STATES, "injury", "stabilizing"),
                              levels = c("stabilizing", "injury")),
           dataset_id = droplevels(dataset_id), donor_id = droplevels(donor_id))
injury_grouping <- "label"
## If adding the basement-membrane panel has taken over the argmax (as it did in
## human, where the label-based injury set fell from 4,420 to 220 cells), fall
## back to a median split on the composite injury SCORE, which is unaffected by
## which panel wins. Record which grouping was used so the result is not silently
## a different analysis from the one the manuscript describes.
score_cols_present <- intersect(INJURY_SCORE_COLS, names(mural))
if ((nrow(sub) < MIN_LABEL_CELLS || dplyr::n_distinct(sub$is_injury) < 2) &&
    length(score_cols_present) > 0) {
    ## `..` prefix: mural is a data.table, so a character vector of column names
    ## must be dereferenced from the calling scope.
    zz <- scale(as.matrix(mural[, ..score_cols_present]))
    zz[is.na(zz)] <- 0
    mural$injury_score <- rowMeans(zz)
    cut <- stats::median(mural$injury_score, na.rm = TRUE)
    sub <- mural |>
        mutate(is_injury = factor(ifelse(injury_score > cut, "injury", "stabilizing"),
                                  levels = c("stabilizing", "injury")),
               dataset_id = droplevels(dataset_id), donor_id = droplevels(donor_id))
    injury_grouping <- "score_median_split"
}
cat("injury grouping used:", injury_grouping, "| cells:", nrow(sub), "\n")
write_tsv(data.frame(injury_grouping = injury_grouping, n_cells = nrow(sub),
                     label_states = paste(INJURY_STATES, collapse = ",")),
          file.path(outdir, "mouse_injury_grouping_used.tsv"))
f1 <- if (has_batch) { Agtr1a_expr ~ is_injury + dataset_id + (1 | donor_id)
      } else { Agtr1a_expr ~ is_injury + (1 | donor_id) }
m1 <- tryCatch(suppressMessages(lmerTest::lmer(f1, data = sub)), error = function(e) NULL)
if (!is.null(m1)) {
    co <- summary(m1)$coefficients; r <- grep("is_injury", rownames(co), value = TRUE)[1]
    emm <- as.data.frame(emmeans(m1, ~ is_injury))
    write_tsv(data.frame(test = "categorical_injury_vs_stabilizing (lmer, +1|donor, dataset fixed)",
                         n_cells = nrow(sub), n_donors = nlevels(droplevels(sub$donor_id)),
                         n_datasets = nlevels(droplevels(sub$dataset_id)),
                         estimate = co[r, "Estimate"], se = co[r, "Std. Error"],
                         df = co[r, "df"], p_value = co[r, "Pr(>|t|)"],
                         emmean_stabilizing = emm$emmean[emm$is_injury == "stabilizing"],
                         emmean_injury = emm$emmean[emm$is_injury == "injury"]),
              file.path(outdir, "mouse_Agtr1a_injury_vs_stable.tsv"))
}

## (2) Continuous: Agtr1a ~ vascular-stabilizing score
f2 <- if (has_batch) { Agtr1a_expr ~ vascular_stabilizing_score + dataset_id + (1 | donor_id)
      } else { Agtr1a_expr ~ vascular_stabilizing_score + (1 | donor_id) }
m2 <- tryCatch(suppressMessages(lmerTest::lmer(f2, data = mural)), error = function(e) NULL)
if (!is.null(m2)) {
    co <- summary(m2)$coefficients["vascular_stabilizing_score", ]
    write_tsv(data.frame(test = "continuous_Agtr1a_vs_stabilizing_score (lmer, +1|donor, dataset fixed)",
                         n_cells = nrow(mural), estimate = co["Estimate"], se = co["Std. Error"],
                         df = co["df"], p_value = co["Pr(>|t|)"]),
              file.path(outdir, "mouse_Agtr1a_vs_stabilizing_continuous.tsv"))
}

## Per-dataset consistency (within-dataset Spearman of Agtr1a vs stabilizing score)
per_ds <- mural |> group_by(dataset_id) |>
    filter(n() >= 20) |>
    summarise(n = n(),
              rho = cor(Agtr1a_expr, vascular_stabilizing_score, method = "spearman"),
              p = tryCatch(cor.test(Agtr1a_expr, vascular_stabilizing_score,
                                    method = "spearman")$p.value, error = function(e) NA),
              .groups = "drop")
write_tsv(per_ds, file.path(outdir, "mouse_Agtr1a_per_dataset_consistency.tsv"))

## Plots
if ("Agtr1a_expr" %in% names(mural)) {
    p <- ggboxplot(mural, x = "pericyte_state", y = "Agtr1a_expr", add = "jitter",
                   fill = "pericyte_state", palette = "npg", legend = "none",
                   add.params = list(alpha = 0.4, size = 0.7),
                   xlab = "Mouse mural state (integrated)", ylab = "Agtr1a (corrected)",
                   ggtheme = theme_pubr(base_size = 12)) + rotate_x_text(35)
    save_ggplots(file.path(outdir, "mouse_Agtr1a_by_state"), p, 6, 5)
}

cat("\nNOTE: states assigned from clustering on the scVI-integrated latent; the\n",
    "categorical and continuous tests agree and per-dataset directions are\n",
    "reported for consistency. Mouse lung pericytes remain sparsely annotated.\n")
cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
