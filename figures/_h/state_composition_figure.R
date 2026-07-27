## Supplementary figure S11: discrete pericyte-state composition does not differ
## across disease groups.
##
## This is a deliberate null. The disease signal in this study is carried by the
## CONTINUOUS injury-program score (see the disease main figure), not by a shift in
## how many pericytes fall into each discrete cluster. Showing the composition null
## explicitly forecloses the reading that the continuous result is a relabelled
## abundance change, and it bounds what the discrete state model can be asked to do.
##
##   A  donor fractions of the six stable pericyte clusters, by disease group
##   B  donor fractions of the three dominant programs
##   C  grouped injury-associated state fraction
##   D  forest of the disease contrasts across every cluster, program and the
##      grouped fraction, with BH-adjusted significance
##
## Panel D is the efficient statement of the null; A-C show the donor-level spread
## behind it. Models are donor-level ANCOVA (frac ~ disease_group + age + sex) fit by
## pericyte_states/_h/01.state_stats.R; only donors with >= 20 pericytes are included.
##
## No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork)
})

source("../_h/_fig_common.R")

SD <- P("pericyte_states", "_m", "stats_data")
GRP_LEVELS <- c("Healthy", "Fibrotic_ILD", "Other")

PROG_NICE <- c(vascular_stabilizing = "Vascular-stabilizing",
               basement_membrane = "Basement-membrane",
               activated_migratory = "Activated/migratory")
PROG_ORDER <- names(PROG_NICE)

## Cluster -> program family, so panel A can be read against panel B
## (pericyte_states/_m/annotations/state_program_map.tsv).
CLUST_PROG <- c(`0` = "vascular_stabilizing", `2` = "vascular_stabilizing",
                `1` = "basement_membrane", `3` = "basement_membrane",
                `5` = "basement_membrane", `4` = "activated_migratory")
CLUST_ORDER <- c("0", "2", "1", "3", "5", "4")

dx <- function(x) factor(x, levels = GRP_LEVELS)

## ---- donor-level composition (written by 01.state_stats.R) ---------------
st <- fread(file.path(SD, "composition_state_by_donor.tsv"))[, level := as.character(level)]
pg <- fread(file.path(SD, "composition_program_by_donor.tsv"))[, level := as.character(level)]
ij <- fread(file.path(SD, "injury_fraction_by_donor.tsv"))
for (d in list(st, pg, ij)) d[, disease_group := dx(disease_group)]

st <- st[level %chin% CLUST_ORDER]
st[, level := factor(level, levels = CLUST_ORDER)]
st[, facet := sprintf("Cluster %s\n(%s)", level, PROG_NICE[CLUST_PROG[as.character(level)]])]
st[, facet := factor(facet, levels = unique(facet[order(level)]))]

pg <- pg[level %chin% PROG_ORDER]
pg[, facet := factor(PROG_NICE[level], levels = unname(PROG_NICE))]

## ---- shared box + jitter + marginal mean ---------------------------------
## Boxes show the donor spread; the white diamond is the age/sex-adjusted marginal
## mean with its 95% CI, i.e. the quantity the ANCOVA actually compares.
comp_panel <- function(d, yvar, emm, ylab, ncol) {
    ggplot(d, aes(disease_group, .data[[yvar]])) +
        geom_boxplot(aes(fill = disease_group), width = 0.62, outlier.shape = NA,
                     alpha = 0.75, linewidth = 0.25) +
        geom_jitter(width = 0.13, size = 0.4, alpha = 0.45, colour = "grey25") +
        {if (!is.null(emm))
            list(geom_errorbar(data = emm,
                               aes(x = disease_group, ymin = lower.CL, ymax = upper.CL),
                               inherit.aes = FALSE, width = 0.16, linewidth = 0.35),
                 geom_point(data = emm, aes(x = disease_group, y = emmean),
                            inherit.aes = FALSE, shape = 23, size = 1.4,
                            fill = "white", stroke = 0.35))} +
        facet_wrap(~ facet, ncol = ncol, scales = "free_y") +
        scale_fill_manual(values = DISEASE_COL) +
        scale_x_discrete(labels = DISEASE_LABS) +
        labs(x = NULL, y = ylab) +
        theme_ms() +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.8),
              strip.text = element_text(size = 5.6, face = "bold", lineheight = 0.95))
}

## Marginal means live in one file per cluster / program; stitch them back into a
## single facet-keyed table.
read_emmeans <- function(prefix, keys, facet_of) {
    rows <- lapply(keys, function(k) {
        f <- file.path(SD, sprintf("%s_%s_emmeans.tsv", prefix, k))
        if (!file.exists(f)) return(NULL)
        e <- fread(f); e[, key := k]; e
    })
    e <- rbindlist(rows, fill = TRUE)
    if (!nrow(e)) return(NULL)
    e[, disease_group := dx(disease_group)]
    e <- e[!is.na(disease_group)]
    e[, facet := facet_of(key)]
    e
}

emm_st <- read_emmeans("composition_state", CLUST_ORDER, function(k)
    factor(sprintf("Cluster %s\n(%s)", k, PROG_NICE[CLUST_PROG[k]]),
           levels = levels(st$facet)))
emm_pg <- read_emmeans("composition_program", PROG_ORDER, function(k)
    factor(PROG_NICE[k], levels = unname(PROG_NICE)))

pA <- comp_panel(st, "frac", emm_st, "Fraction of donor's pericytes", 6)
pB <- comp_panel(pg, "frac", emm_pg, "Fraction of donor's pericytes", 3)

## ===== Panel C -- grouped injury-associated state fraction ================
emm_ij <- fread(file.path(SD, "injury_fraction_emmeans.tsv"))
emm_ij[, disease_group := dx(disease_group)]
emm_ij <- emm_ij[!is.na(disease_group)]

pC <- ggplot(ij, aes(disease_group, injury_frac)) +
    geom_boxplot(aes(fill = disease_group), width = 0.58, outlier.shape = NA,
                 alpha = 0.75, linewidth = 0.25) +
    geom_jitter(width = 0.12, size = 0.55, alpha = 0.5, colour = "grey25") +
    geom_errorbar(data = emm_ij, aes(x = disease_group, ymin = lower.CL, ymax = upper.CL),
                  inherit.aes = FALSE, width = 0.15, linewidth = 0.35) +
    geom_point(data = emm_ij, aes(x = disease_group, y = emmean), inherit.aes = FALSE,
               shape = 23, size = 1.6, fill = "white", stroke = 0.4) +
    scale_fill_manual(values = DISEASE_COL) +
    scale_x_discrete(labels = DISEASE_LABS) +
    labs(x = NULL, y = "Injury-associated state fraction\n(per donor)") +
    theme_ms() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 6.5))

## ===== Panel D -- forest of disease contrasts =============================
## lower.CL/upper.CL are NOMINAL 95% intervals while p.value is BH-adjusted across
## contrasts within each response (see posthoc_with_ci in 01.state_stats.R), so an
## interval may exclude zero while the adjusted p does not.
read_posthoc <- function(prefix, keys, labeller, family) {
    rows <- lapply(keys, function(k) {
        f <- file.path(SD, sprintf("%s_%s_posthoc.tsv", prefix, k))
        if (!file.exists(f)) return(NULL)
        d <- fread(f); d[, `:=`(key = k, response = labeller(k), family = family)]; d
    })
    rbindlist(rows, fill = TRUE)
}

fo <- rbindlist(list(
    read_posthoc("composition_state", CLUST_ORDER,
                 function(k) sprintf("Cluster %s", k), "Stable cluster"),
    read_posthoc("composition_program", PROG_ORDER,
                 function(k) unname(PROG_NICE[k]), "Dominant program"),
    {d <- fread(file.path(SD, "injury_fraction_posthoc.tsv"))
     d[, `:=`(key = "injury", response = "Injury-associated fraction",
              family = "Grouped")]; d}
), fill = TRUE)

## Only the contrasts against Healthy; the Fibrotic-vs-Other contrast is not a
## question this figure asks.
fo <- fo[grepl("^Healthy - ", contrast)]
## Report the effect as (group - Healthy) so a positive value means "higher in
## disease", which is how the caption reads.
fo[, `:=`(estimate = -estimate, lo = -upper.CL, hi = -lower.CL,
          comparison = sub("^Healthy - ", "", contrast))]
fo[, comparison := factor(DISEASE_LABS[comparison], levels = DISEASE_LABS[c("Fibrotic_ILD", "Other")])]
fo[, family := factor(family, levels = c("Stable cluster", "Dominant program", "Grouped"))]
resp_order <- unique(fo[order(family, match(key, c(CLUST_ORDER, PROG_ORDER, "injury"))), response])
fo[, response := factor(response, levels = rev(resp_order))]
fo[, sig := fifelse(p.value < 0.05, "BH p < 0.05", "n.s.")]

pD <- ggplot(fo, aes(estimate, response, colour = comparison, shape = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35, linetype = 2) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4,
                   position = position_dodge(width = 0.55)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.55)) +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = unname(DISEASE_COL[c("Fibrotic_ILD", "Other")]), name = NULL) +
    scale_shape_manual(values = c(`BH p < 0.05` = 16, `n.s.` = 1), name = NULL) +
    labs(x = "Difference in donor fraction vs Healthy\n(age/sex-adjusted, 95% CI)", y = NULL) +
    theme_ms() +
    theme(axis.text.y = element_text(size = 6),
          strip.text.y = element_text(size = 5.6, face = "bold", angle = 0),
          legend.position = "top", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 6), legend.box = "horizontal")

## ---- assemble ------------------------------------------------------------
fig <- pA / (pB | pC) / pD +
    plot_layout(heights = c(1, 1, 1.45)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figureS_state_composition", fig, 8.5, 10.0)

cat("Wrote figureS_state_composition to", OUT, "\n")
cat(sprintf("  donors: %d; contrasts in forest: %d; any BH p < 0.05: %s\n",
            uniqueN(ij$donor_id), nrow(fo), any(fo$p.value < 0.05)))
cat("\nReproducibility information:\n"); sessioninfo::session_info()
