## Disease main figure (Circ Research revision): the reviewer-defensible
## Healthy-vs-Fibrotic/ILD pericyte injury-program result from
## disease_association/03. Three panels, no in-panel titles (interpretation in
## the caption), shared manuscript visual language.
##
##   A  within-study random-effects meta-analysis forest (the headline: the
##      Fibrotic-Healthy effect is consistent WITHIN studies, not a between-study
##      batch artifact).
##   B  program specificity -- which injury programs drive the composite (fibrillar
##      fibroblast-like + activated/migratory, not a global shift).
##   C  donor-level endpoint by disease group with sex+study-adjusted marginal
##      means (the pooled contrast the LMM tests).
suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P    <- function(...) file.path(ROOT, ...)
OUT  <- P("figures", "mechanism"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SD   <- P("disease_association", "_m", "mixed_model_forest")

DISEASE_COL <- c(Healthy = "#0072B2", Fibrotic_ILD = "#D55E00")
DISEASE_LAB <- c(Healthy = "Healthy", Fibrotic_ILD = "Fibrotic/ILD")

theme_ms <- function(base = 8) {
    theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              legend.position = "none", strip.background = element_blank())
}
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}
fmt_p <- function(p) if (p < 1e-3) sprintf("P = %.1e", p) else sprintf("P = %.3f", p)

## ---- load 03 outputs -------------------------------------------------------
study  <- fread(file.path(SD, "forest_per_study.tsv"))
pool   <- fread(file.path(SD, "forest_pooled_RE.tsv"))
comp   <- fread(file.path(SD, "component_effects.tsv"))
prim   <- fread(file.path(SD, "primary_effect.tsv"))
emm    <- fread(file.path(SD, "primary_emmeans.tsv"))
donor  <- fread(file.path(SD, "donor_endpoint_table.tsv"))

## ===========================================================================
## Panel A -- within-study meta forest
## ===========================================================================
study_lab <- sprintf("%s  (%d H / %d F)", study$dataset, study$nH, study$nF)
pool_lab  <- sprintf("RE pooled  (I2 = %.0f%%)", pool$I2)
fp <- rbind(
    data.table(label = study_lab, y = study$yi, lo = study$ci_lo, hi = study$ci_hi,
               w = study$weight_pct, kind = "study"),
    data.table(label = pool_lab, y = pool$estimate, lo = pool$ci_lo, hi = pool$ci_hi,
               w = max(study$weight_pct), kind = "pooled"))
fp[, label := factor(label, levels = rev(c(study_lab[order(study$yi)], pool_lab)))]

pA <- ggplot(fp, aes(y = label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = kind), height = 0.2, linewidth = 0.5) +
    geom_point(aes(x = y, colour = kind, size = w, shape = kind)) +
    annotate("text", x = pool$estimate, y = 0.45,
             label = sprintf("+%.2f SD   %s", pool$estimate, fmt_p(pool$p_pooled)),
             size = 2.5, fontface = "bold", colour = "#D55E00") +
    scale_colour_manual(values = c(study = "#0072B2", pooled = "#D55E00"), guide = "none") +
    scale_shape_manual(values = c(study = 16, pooled = 18), guide = "none") +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    coord_cartesian(clip = "off", ylim = c(0.8, length(levels(fp$label)))) +
    labs(x = "Fibrotic/ILD minus Healthy\ninjury-program score (SD units)", y = NULL) +
    theme_ms() + theme(panel.grid.major.y = element_blank())

## ===========================================================================
## Panel B -- program specificity (component effects, +/-95% CI)
## ===========================================================================
nice_prog <- c(z_inflammatory = "Inflammatory", z_activated_migratory = "Activated/migratory",
               z_fibroblast_like = "Fibrillar\nfibroblast-like",
               z_vascular_stabilizing = "Vascular-\nstabilizing")
prog_col  <- c(z_inflammatory = "#E69F00", z_activated_migratory = "#CC79A7",
               z_fibroblast_like = "#D55E00", z_vascular_stabilizing = "#0072B2")
comp <- comp[response %in% names(nice_prog)]
comp[, `:=`(lo = beta_fib_minus_healthy - 1.96 * se,
            hi = beta_fib_minus_healthy + 1.96 * se,
            prog = factor(nice_prog[response], levels = nice_prog[order(comp$beta_fib_minus_healthy)]),
            sig = ifelse(p < 0.05, "*", ""))]

pB <- ggplot(comp, aes(beta_fib_minus_healthy, prog, colour = response)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.6) +
    geom_text(aes(x = hi, label = sig), hjust = -0.4, vjust = 0.75, size = 4, show.legend = FALSE) +
    scale_colour_manual(values = prog_col, guide = "none") +
    labs(x = "Fibrotic/ILD minus Healthy (SD units)", y = NULL) +
    theme_ms() + theme(panel.grid.major.y = element_blank())

## ===========================================================================
## Panel C -- donor endpoint by disease group + sex/study-adjusted marginal means
## ===========================================================================
donor <- donor[disease_group %in% names(DISEASE_COL)]
donor[, dx := factor(disease_group, levels = names(DISEASE_COL))]
emm[, dx := factor(disease_group, levels = names(DISEASE_COL))]
set.seed(13)
p_prim <- prim$p[1]

pC <- ggplot(donor, aes(dx, injury_program_score)) +
    geom_hline(yintercept = 0, linetype = 3, colour = "grey70", linewidth = 0.3) +
    geom_jitter(aes(colour = dx), width = 0.13, height = 0, alpha = 0.5, size = 1.3) +
    geom_errorbar(data = emm, aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  width = 0.16, linewidth = 0.6, colour = "black", inherit.aes = TRUE) +
    geom_point(data = emm, aes(y = emmean), size = 2.4, colour = "black") +
    annotate("segment", x = 1, xend = 2,
             y = max(donor$injury_program_score) * 1.05,
             yend = max(donor$injury_program_score) * 1.05, linewidth = 0.3) +
    annotate("text", x = 1.5, y = max(donor$injury_program_score) * 1.12,
             label = sprintf("+%.2f SD   %s", prim$beta_fib_minus_healthy[1], fmt_p(p_prim)),
             size = 2.5, fontface = "bold") +
    scale_colour_manual(values = DISEASE_COL, guide = "none") +
    scale_x_discrete(labels = DISEASE_LAB) +
    labs(x = NULL, y = "Injury-program score (SD units)") +
    theme_ms()

## ---- assemble --------------------------------------------------------------
fig <- pA / (pB | pC) +
    plot_layout(heights = c(1, 1.05)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figure_disease_main", fig, 7.2, 6.0)

cat("wrote", file.path(OUT, "figure_disease_main.{pdf,svg,png}"), "\n")
cat("\nReproducibility information:\n"); Sys.time(); options(width = 120); sessioninfo::session_info()
