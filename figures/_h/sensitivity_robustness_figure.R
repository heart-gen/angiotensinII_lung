## Supplementary figure: smoking + cohort robustness of the donor-level disease
## associations (sensitivity module). Four panels, left-to-right / top-to-bottom:
##   A  smoking metadata is recorded only for Healthy donors (the confound that
##      makes a smoking-STRATIFIED disease contrast inestimable);
##   B  among donors that DO carry a smoking label, smoking shows no gradient in
##      the injury / AGTR1 read-outs (the phenotype is not a smoking artifact);
##   C  the Healthy-vs-Fibrotic/ILD disease effect is unchanged when smoking is
##      added as a covariate (marginal means +/- 95% CI, base vs +smoking);
##   D  leave-one-study-out: the Fibrotic/ILD effect is stable across cohorts.
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P <- function(...) file.path(ROOT, ...)
OUT <- P("figures", "mechanism"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SD <- function(f) P("sensitivity", "_m", "stats_data", f)

## ---- shared visual language (matches manuscript_mechanism_figure.R) -----
DISEASE_LEVELS <- c("Healthy", "COPD", "Fibrotic_ILD", "Other")
DISEASE_LABS   <- c(Healthy = "Healthy", COPD = "COPD",
                    Fibrotic_ILD = "Fibrotic/ILD", Other = "Other")
DISEASE_COL <- c(Healthy = "#0072B2", COPD = "#E69F00",
                 Fibrotic_ILD = "#D55E00", Other = "#999999")
SMK_LEVELS <- c("never", "former", "active")
SMK_LABS   <- c(never = "Never", former = "Former", active = "Active")
SMK_COL    <- c(never = "#0072B2", former = "#E69F00", active = "#D55E00")
RESP_LEVELS <- c("injury_frac", "injury_stromal_score", "niche_index", "AGTR1_pos_frac")
RESP_LABS   <- c(injury_frac = "Injury-state fraction",
                 injury_stromal_score = "Injury-stromal score",
                 niche_index = "Niche-stability index",
                 AGTR1_pos_frac = "AGTR1+ fraction")
resp_factor <- function(x) factor(x, levels = RESP_LEVELS)

theme_ms <- function(base = 8) {
    theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              legend.position = "none",
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", size = 7))
}
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}

## ===== Panel A — smoking metadata availability by disease ================
avail <- fread(SD("smoking_availability_by_disease.tsv")) %>%
    mutate(disease_group = factor(disease_group, levels = DISEASE_LEVELS)) %>%
    filter(!is.na(disease_group)) %>%
    transmute(disease_group, Labelled = n_with_smoking,
              Missing = n_donors - n_with_smoking) %>%
    pivot_longer(c(Labelled, Missing), names_to = "status", values_to = "n") %>%
    mutate(status = factor(status, levels = c("Missing", "Labelled")))
n_lab <- avail %>% filter(status == "Labelled")
pA <- ggplot(avail, aes(disease_group, n, fill = status)) +
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

## ===== Panel B — smoking main effect among labelled donors ===============
smk <- fread(SD("smoking_main_effect_healthy.tsv")) %>%
    mutate(smoking = factor(smoking, levels = SMK_LEVELS),
           response = resp_factor(response))
pB <- ggplot(smk, aes(smoking, emmean, colour = smoking)) +
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

## ===== Panel C — disease effect vs +smoking covariate ====================
MODEL_LABS <- c(base = "Base", smoking = "+ Smoking")
cov <- fread(SD("covariate_robustness_emmeans.tsv")) %>%
    filter(model %in% c("base", "smoking"),
           disease_group %in% c("Healthy", "Fibrotic_ILD")) %>%
    mutate(disease_group = factor(disease_group, levels = DISEASE_LEVELS),
           model = factor(model, levels = c("base", "smoking")),
           response = resp_factor(response))
pC <- ggplot(cov, aes(model, emmean, colour = disease_group, group = disease_group)) +
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

## ===== Panel D — leave-one-study-out stability of the Fibrotic/ILD effect =
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

## ---- assemble ----------------------------------------------------------
## top row: narrow availability bar (A) + wide 4-facet smoking panel (B)
row1 <- pA + pB + plot_layout(widths = c(0.7, 3))
fig <- row1 / pC / pD +
    plot_layout(heights = c(1, 1, 1.7)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_sensitivity", fig, 9.0, 9.5)

cat("Wrote figureS_sensitivity to", OUT, "\n")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
