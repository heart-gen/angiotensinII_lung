## Supplementary figure S12: robustness and limitations of the disease-associated
## injury-stromal signal.
##
## The donor-level claim is that fibrotic/ILD lungs carry a higher pericyte
## injury-stromal score. Four things a reviewer will press on, and where each is
## answered here:
##   - is the effect an artifact of putting AGTR1 into the composite?  (A)
##   - why is the NET niche index a weaker discriminator than its injury half?  (B)
##   - is it confounded by smoking?                                      (C, E, F)
##   - is it driven by one cohort?                                       (D)
##
##   A  injury-stromal score, primary composite vs the +AGTR1-positive-fraction
##      sensitivity composite. Excluding AGTR1 avoids dropout and circularity; it
##      does not create the disease association.
##   B  the vascular-stability component and the net niche index by disease group.
##      Stability RISES alongside injury in fibrotic/ILD lungs, which is exactly why
##      the net index (stability - injury) separates the groups less well.
##   C  smoking metadata is recorded only for Healthy donors -- the confound that
##      makes a smoking-STRATIFIED disease contrast inestimable, not merely weak.
##   D  leave-one-study-out: the Fibrotic/ILD effect is stable across cohorts.
##   E  among donors that DO carry a smoking label, smoking shows no gradient in the
##      injury / AGTR1 read-outs.
##   F  the disease effect is unchanged when smoking is added as a covariate.
##
## Panels A-B come from niche_index/_m, C-F from sensitivity/_m. Both use the same
## donor set. NOTE: disease_association/_m/mixed_model_forest/ contains a similarly
## named smoking-availability table computed on a DIFFERENT donor set and endpoint
## (42 Healthy / 21 labelled, injury_program_score); this figure is entirely the
## niche_index + sensitivity story and quotes 14 of 32. Do not mix the two.
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

SD <- function(f) P("sensitivity", "_m", "stats_data", f)
NI <- function(f) P("niche_index", "_m", "stats_data", f)

SMK_LEVELS <- c("never", "former", "active")
SMK_LABS   <- c(never = "Never", former = "Former", active = "Active")
SMK_COL    <- c(never = "#0072B2", former = "#E69F00", active = "#D55E00")
RESP_LEVELS <- c("injury_frac", "injury_stromal_score", "niche_index", "AGTR1_pos_frac")
RESP_LABS   <- c(injury_frac = "Injury-state fraction",
                 injury_stromal_score = "Injury-stromal score",
                 niche_index = "Niche-stability index",
                 AGTR1_pos_frac = "AGTR1+ fraction")
resp_factor <- function(x) factor(x, levels = RESP_LEVELS)

## The composite scores are donor-level and this figure needs the raw points
## alongside the model's marginal means.
donor <- fread(P("niche_index", "_m", "niche_index_per_donor.tsv.gz"))
donor[, disease_group := factor(map_disease(lung_condition), levels = DISEASE_LEVELS)]

## Marginal means + donor points for one composite, laid out as facets so several
## composites can be compared on a shared panel.
score_panel <- function(specs, ylab) {
    emm <- rbindlist(lapply(names(specs), function(k) {
        e <- fread(NI(paste0(k, "_emmeans.tsv")))
        e[, variant := specs[[k]]]; e
    }), fill = TRUE)
    pts <- rbindlist(lapply(names(specs), function(k)
        donor[, .(disease_group, value = get(k), variant = specs[[k]])]), fill = TRUE)
    lev <- unname(unlist(specs))
    emm[, `:=`(variant = factor(variant, levels = lev),
               disease_group = factor(disease_group, levels = DISEASE_LEVELS))]
    pts[, `:=`(variant = factor(variant, levels = lev),
               disease_group = factor(disease_group, levels = DISEASE_LEVELS))]
    emm <- emm[!is.na(disease_group)]; pts <- pts[!is.na(disease_group) & !is.na(value)]

    ggplot(pts, aes(disease_group, value)) +
        geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
        geom_jitter(aes(colour = disease_group), width = 0.14, size = 0.55, alpha = 0.45) +
        geom_errorbar(data = emm, aes(x = disease_group, ymin = lower.CL, ymax = upper.CL),
                      inherit.aes = FALSE, width = 0.16, linewidth = 0.4) +
        geom_point(data = emm, aes(x = disease_group, y = emmean, fill = disease_group),
                   inherit.aes = FALSE, shape = 23, size = 1.9, stroke = 0.4) +
        facet_wrap(~ variant, nrow = 1) +
        scale_colour_manual(values = DISEASE_COL) +
        scale_fill_manual(values = DISEASE_COL) +
        scale_x_discrete(labels = DISEASE_LABS) +
        labs(x = NULL, y = ylab) +
        theme_ms() +
        theme(axis.text.x = element_text(angle = 25, hjust = 1),
              strip.text = element_text(face = "bold", size = 6.2))
}

## ===== Panel A -- primary vs AGTR1-augmented composite ====================
pA <- score_panel(
    list(injury_stromal_score            = "Primary (AGTR1 excluded)",
         injury_stromal_score_sens_agtr1 = "Sensitivity (+ AGTR1+ fraction)"),
    "Injury-stromal score (z)")

## ===== Panel B -- stability component and the net index ===================
pB <- score_panel(
    list(niche_stability_score = "Vascular-stability component",
         niche_index           = "Net niche-stability index"),
    "Score (z)")

## ===== Panel C -- smoking metadata availability by disease ================
avail <- fread(SD("smoking_availability_by_disease.tsv")) %>%
    mutate(disease_group = factor(disease_group, levels = DISEASE_LEVELS)) %>%
    filter(!is.na(disease_group)) %>%
    transmute(disease_group, Labelled = n_with_smoking,
              Missing = n_donors - n_with_smoking) %>%
    pivot_longer(c(Labelled, Missing), names_to = "status", values_to = "n") %>%
    mutate(status = factor(status, levels = c("Missing", "Labelled")))
n_lab <- avail %>% filter(status == "Labelled")
pC <- ggplot(avail, aes(disease_group, n, fill = status)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_text(data = n_lab, aes(label = ifelse(n > 0, n, "")),
              vjust = -0.4, size = 2.3, colour = "#0072B2", fontface = "bold") +
    scale_fill_manual(values = c(Missing = "grey85", Labelled = "#0072B2"),
                      labels = c(Missing = "No smoking record",
                                 Labelled = "Smoking recorded"), name = NULL) +
    scale_x_discrete(labels = DISEASE_LABS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = NULL, y = "Donors") +
    theme_ms() +
    theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.key.size = unit(3, "mm"), legend.text = element_text(size = 6),
          axis.text.x = element_text(angle = 25, hjust = 1))

## ===== Panel D -- leave-one-study-out stability of the Fibrotic/ILD effect =
loso <- fread(SD("leave_one_study_out.tsv")) %>%
    mutate(response = resp_factor(response),
           lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
           sig = ifelse(p < 0.05, "p < 0.05", "n.s."))
## order studies by the mean effect on the headline injury-state fraction
ord <- loso %>% filter(response == "injury_frac") %>%
    arrange(estimate) %>% pull(dropped_dataset)
loso <- loso %>% mutate(dropped_dataset = factor(dropped_dataset, levels = ord))
pD <- ggplot(loso, aes(estimate, dropped_dataset, colour = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
    geom_point(size = 1.3) +
    facet_wrap(~ response, nrow = 1, scales = "free_x",
               labeller = as_labeller(RESP_LABS)) +
    scale_colour_manual(values = c("p < 0.05" = "#D55E00", "n.s." = "grey60"), name = NULL) +
    labs(x = "Fibrotic/ILD effect (study left out)", y = NULL) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 5.5),
          legend.position = "top", legend.text = element_text(size = 6.5),
          legend.key.size = unit(3, "mm"))

## ===== Panel E -- smoking main effect among labelled donors ===============
smk <- fread(SD("smoking_main_effect_healthy.tsv")) %>%
    mutate(smoking = factor(smoking, levels = SMK_LEVELS),
           response = resp_factor(response))
pE <- ggplot(smk, aes(smoking, emmean, colour = smoking)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.18, linewidth = 0.4) +
    geom_point(size = 1.7) +
    facet_wrap(~ response, nrow = 1, scales = "free_y",
               labeller = as_labeller(RESP_LABS)) +
    scale_colour_manual(values = SMK_COL) +
    scale_x_discrete(labels = SMK_LABS) +
    labs(x = NULL, y = "Estimated marginal mean") +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          strip.text = element_text(face = "bold", size = 5.8),
          strip.clip = "off")

## ===== Panel F -- disease effect vs +smoking covariate ====================
MODEL_LABS <- c(base = "Base", smoking = "+ Smoking")
cov <- fread(SD("covariate_robustness_emmeans.tsv")) %>%
    filter(model %in% c("base", "smoking"),
           disease_group %in% c("Healthy", "Fibrotic_ILD")) %>%
    mutate(disease_group = factor(disease_group, levels = DISEASE_LEVELS),
           model = factor(model, levels = c("base", "smoking")),
           response = resp_factor(response))
pF <- ggplot(cov, aes(model, emmean, colour = disease_group, group = disease_group)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_line(linewidth = 0.4, position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.16,
                  linewidth = 0.4, position = position_dodge(width = 0.4)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.4)) +
    facet_wrap(~ response, nrow = 1, scales = "free_y",
               labeller = as_labeller(RESP_LABS)) +
    scale_colour_manual(values = DISEASE_COL, labels = DISEASE_LABS, name = NULL) +
    scale_x_discrete(labels = MODEL_LABS) +
    labs(x = NULL, y = "Estimated marginal mean") +
    theme_ms() +
    theme(legend.position = "top", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6.5),
          axis.text.x = element_text(angle = 25, hjust = 1))

## ---- assemble ------------------------------------------------------------
## Row 1: the two composite-comparison panels (A, B).
## Row 2: narrow availability bar (C) beside the wide 4-facet LOSO forest (D).
## Rows 3-4: the two smoking panels retained from the previous version (E, F).
row1 <- pA + pB + plot_layout(widths = c(1, 1))
row2 <- pC + pD + plot_layout(widths = c(0.55, 3))
fig <- row1 / row2 / pE / pF +
    plot_layout(heights = c(1, 1.7, 1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_sensitivity", fig, 9.0, 14.0)

cat("Wrote figureS_sensitivity to", OUT, "\n")
n_loso <- loso %>% filter(response == "injury_stromal_score")
cat(sprintf("  LOSO injury_stromal_score: %d refits, %d positive, %d with p < 0.05\n",
            nrow(n_loso), sum(n_loso$estimate > 0), sum(n_loso$p < 0.05)))
cat("\nReproducibility information:\n"); sessioninfo::session_info()
