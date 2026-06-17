## Manuscript-quality mechanistic figures (Circ Research revision).
##
## Main figure: donor-level disease effects (niche-stability index, injury-stromal
## score, AT1R-AT2R balance) + state composition + state-resolved balance +
## continuum trends. Supplements: alluvial (cluster->category->role) and mouse
## cross-species. No in-panel titles; interpretation belongs in the caption.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(tidyr)
    library(ggplot2); library(ggpubr); library(patchwork); library(ggalluvial)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P <- function(...) file.path(ROOT, ...)
OUT <- P("figures", "mechanism"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- shared visual language --------------------------------------------
DISEASE_LEVELS <- c("Healthy", "COPD", "Fibrotic_ILD", "Other")
DISEASE_LABS   <- c(Healthy = "Healthy", COPD = "COPD",
                    Fibrotic_ILD = "Fibrotic/ILD", Other = "Other")
DISEASE_COL <- c(Healthy = "#0072B2", COPD = "#E69F00",
                 Fibrotic_ILD = "#D55E00", Other = "#999999")
STATE_LEVELS <- c("vascular_stabilizing", "synthetic_contractile",
                  "activated_migratory", "inflammatory", "fibroblast_like")
STATE_LABS <- c(vascular_stabilizing = "Vascular-\nstabilizing",
                synthetic_contractile = "Synthetic/\ncontractile",
                activated_migratory = "Activated/\nmigratory",
                inflammatory = "Inflammatory",
                fibroblast_like = "Fibroblast-\nlike")
STATE_LABS1 <- c(vascular_stabilizing = "Vascular-stabilizing",
                 synthetic_contractile = "Synthetic/contractile",
                 activated_migratory = "Activated/migratory",
                 inflammatory = "Inflammatory",
                 fibroblast_like = "Fibroblast-like")
STATE_COL <- c(vascular_stabilizing = "#0072B2", synthetic_contractile = "#009E73",
               activated_migratory = "#CC79A7", inflammatory = "#E69F00",
               fibroblast_like = "#D55E00")

theme_ms <- function(base = 8) {
    theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              legend.position = "none",
              strip.background = element_blank())
}

map_disease <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc) ~ "Healthy", lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc,
              ignore.case = TRUE) ~ "Fibrotic_ILD", TRUE ~ "Other")
}
dx_factor <- function(x) {
    f <- factor(map_disease(x), levels = DISEASE_LEVELS)
    droplevels(f)
}
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}

## key contrast vs Healthy for the box panels
dx_comparisons <- function(levs)
    Filter(function(p) all(p %in% levs),
           list(c("Healthy", "Fibrotic_ILD"), c("Healthy", "COPD")))

box_by_disease <- function(df, yvar, ylab) {
    df <- df %>% filter(!is.na(.data[[yvar]]), !is.na(disease_group))
    levs <- levels(droplevels(df$disease_group))
    ggplot(df, aes(disease_group, .data[[yvar]], fill = disease_group)) +
        geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
        geom_jitter(width = 0.12, size = 0.5, alpha = 0.5, colour = "grey20") +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6,
                     fill = "white", stroke = 0.4) +
        stat_compare_means(comparisons = dx_comparisons(levs), method = "wilcox.test",
                           size = 2.4, tip.length = 0.01, bracket.size = 0.25) +
        scale_fill_manual(values = DISEASE_COL) +
        scale_x_discrete(labels = DISEASE_LABS) +
        labs(x = NULL, y = ylab) + theme_ms()
}

## ---- load donor-level data ---------------------------------------------
niche <- fread(P("niche_index", "_m", "niche_index_per_donor.tsv.gz")) %>%
    mutate(disease_group = dx_factor(lung_condition))

## NVU-pattern model: the interpretable program is in `state_program` (the six
## stable `pericyte_state` clusters collapse onto three programs); balance panels
## key on `state_program`.
bal_cell <- fread(P("pathway_balance", "_m", "pathway_balance_metadata.tsv.gz"))
INJURY <- c("inflammatory", "fibroblast_like", "activated_migratory")
bal_donor <- bal_cell %>%
    filter(state_program %in% INJURY) %>%
    group_by(donor_id, lung_condition) %>%
    summarise(balance = mean(AT1R_AT2R_balance, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 10) %>% mutate(disease_group = dx_factor(lung_condition))
bal_state <- bal_cell %>%
    group_by(donor_id, state_program) %>%
    summarise(balance = mean(AT1R_AT2R_balance, na.rm = TRUE), n = n(), .groups = "drop") %>%
    filter(n >= 5) %>%
    mutate(pericyte_state = factor(state_program, levels = STATE_LEVELS)) %>%
    filter(!is.na(pericyte_state))

## ---- Panels -------------------------------------------------------------
pA <- box_by_disease(niche, "niche_index", "Niche-stability index")
pB <- box_by_disease(niche, "injury_stromal_score", "Injury-stromal score")
## AT1R-AT2R balance by disease is retained as a COROLLARY of the injury-stromal
## program (it is redundant with injury intensity: disease effect collapses when
## adjusting for the injury-stromal score; see MECHANISM_ANALYSES). The by-state
## balance panel (NS) is moved to the supplement (figureS_balance_by_state).
pC <- box_by_disease(bal_donor, "balance", "AT1R-AT2R balance (corollary)")

## D: state-program composition by disease (mean donor fraction, stacked).
## Computed from the cell-level metadata (donor -> state_program fractions), since
## the per-donor niche file carries only the collapsed injury_frac under the
## NVU-pattern model.
comp_donor <- bal_cell %>%
    filter(!is.na(state_program)) %>%
    group_by(donor_id, lung_condition, state_program) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(donor_id) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
    mutate(disease_group = dx_factor(lung_condition),
           state = factor(state_program, levels = STATE_LEVELS))
comp <- comp_donor %>%
    filter(!is.na(disease_group), !is.na(state)) %>%
    group_by(disease_group, state) %>%
    summarise(frac = mean(frac, na.rm = TRUE), .groups = "drop")
pD <- ggplot(comp, aes(disease_group, frac, fill = state)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = STATE_COL, labels = STATE_LABS, name = NULL) +
    scale_x_discrete(labels = DISEASE_LABS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Mean state fraction") +
    theme_ms() + theme(legend.position = "right",
                       legend.text = element_text(size = 5.5),
                       legend.key.size = unit(2.6, "mm"),
                       axis.text.x = element_text(angle = 30, hjust = 1))

## (SUPPLEMENT) AT1R-AT2R balance across pericyte states -- program differences are
## NOT significant (smallest pairwise p = 0.067); demoted out of the main figure.
pE <- ggplot(bal_state, aes(pericyte_state, balance, fill = pericyte_state)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 0.4, alpha = 0.4, colour = "grey20") +
    scale_fill_manual(values = STATE_COL) +
    scale_x_discrete(limits = rev(STATE_LEVELS), labels = STATE_LABS1) +
    labs(x = NULL, y = "AT1R-AT2R balance",
         title = "AT1R-AT2R balance by pericyte program (n.s., smallest p = 0.067)") +
    coord_flip() + theme_ms() + theme(plot.title = element_text(size = 7))
save_fig("figureS_balance_by_state", pE, 4.2, 3.2)

## ---- Supplement: ACTA2 contractile-identity control --------------------
## Benchmark (not a pillar): is AGTR1 simply canonical ACTA2+ contractile mural
## identity? (A) Centered donor-aware program emmeans for ACTA2 (raw) vs AGTR1 (raw)
## vs AGTR1 (denoised): the raw markers share a vascular-stabilizing-high pattern
## that the denoised AGTR1 does NOT. (B) Donor x program pseudobulk Pearson
## correlation (95% CI) of AGTR1 (raw vs denoised) against ACTA2 and the
## leave-ACTA2-out contractile program: the raw coupling collapses to ~0 / n.s.
## after denoising, so AGTR1 is not reducible to contractile identity.
SD <- P("pericyte_states", "_m", "stats_data")
acta2_f <- file.path(SD, "acta2_by_program_emmeans.tsv")
agtr1_f <- file.path(SD, "agtr1_lenses_by_program_emmeans.tsv")
cor_f   <- file.path(SD, "acta2_control_correlations.tsv")
if (all(file.exists(c(acta2_f, agtr1_f, cor_f)))) {
    PROG_LEVELS <- c("vascular_stabilizing", "fibroblast_like", "activated_migratory")
    READOUT_LABS <- c(ACTA2_expr = "ACTA2 (raw)", AGTR1_expr = "AGTR1 (raw)",
                      AGTR1_scvi = "AGTR1 (denoised)")
    READOUT_COL  <- c("ACTA2 (raw)" = "#009E73", "AGTR1 (raw)" = "#56B4E9",
                      "AGTR1 (denoised)" = "#D55E00")

    emm <- rbind(fread(acta2_f)[lens == "ACTA2_expr"],
                 fread(agtr1_f)[lens %in% c("AGTR1_expr", "AGTR1_scvi")]) %>%
        mutate(readout = factor(READOUT_LABS[lens], levels = READOUT_LABS),
               program = factor(state_program, levels = PROG_LEVELS)) %>%
        filter(!is.na(program)) %>%
        group_by(readout) %>% mutate(centered = emmean - mean(emmean)) %>% ungroup()
    sA <- ggplot(emm, aes(program, centered, colour = readout, group = readout)) +
        geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
        geom_line(linewidth = 0.6) +
        geom_errorbar(aes(ymin = centered - SE, ymax = centered + SE),
                      width = 0.12, linewidth = 0.4) +
        geom_point(size = 1.9) +
        scale_colour_manual(values = READOUT_COL, name = NULL) +
        scale_x_discrete(labels = STATE_LABS[PROG_LEVELS]) +
        labs(x = NULL, y = "Centered donor-aware emmean") +
        theme_ms() + theme(legend.position = c(0.5, 0.12),
                           legend.background = element_blank(),
                           legend.text = element_text(size = 6),
                           legend.key.size = unit(3, "mm"))

    ## (B) pseudobulk correlation forest with Fisher-z 95% CI
    TARGET_LABS <- c(ACTA2_expr = "ACTA2", synth_contr_noACTA2 = "Contractile\nprogram (-ACTA2)")
    cor_dt <- fread(cor_f)[x %in% c("AGTR1_expr", "AGTR1_scvi") &
                           y %in% names(TARGET_LABS)] %>%
        mutate(lens = factor(READOUT_LABS[x], levels = READOUT_LABS[c("AGTR1_expr","AGTR1_scvi")]),
               target = factor(TARGET_LABS[y], levels = rev(TARGET_LABS)),
               z = atanh(pearson_r), se = 1 / sqrt(n - 3),
               lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
               sig = ifelse(pearson_p < 0.05, "p < 0.05", "n.s."))
    sB <- ggplot(cor_dt, aes(pearson_r, target, colour = lens, shape = sig)) +
        geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
        geom_errorbarh(aes(xmin = lo, xmax = hi),
                       position = position_dodge(width = 0.55), height = 0.18, linewidth = 0.4) +
        geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
        scale_colour_manual(values = READOUT_COL, name = NULL) +
        scale_shape_manual(values = c("p < 0.05" = 16, "n.s." = 1), name = NULL) +
        labs(x = "Pearson r (donor x program pseudobulk, 95% CI)", y = NULL) +
        theme_ms() + theme(legend.position = "right",
                           legend.text = element_text(size = 6),
                           legend.key.size = unit(3, "mm"),
                           legend.spacing.y = unit(1, "mm"))

    pACTA2 <- (sA | sB) + plot_layout(widths = c(1, 1.25)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(face = "bold", size = 10))
    save_fig("figureS_acta2_control", pACTA2, 7.4, 3.2)

    ## source-data table (figure underlying values)
    tab_emm <- emm %>%
        transmute(panel = "A", readout, program = state_program,
                  emmean, SE, centered)
    tab_cor <- cor_dt %>%
        transmute(panel = "B", readout = lens,
                  target = gsub("\n", " ", as.character(target)),
                  n, pearson_r, ci_low = lo, ci_high = hi, pearson_p)
    fwrite(rbindlist(list(tab_emm, tab_cor), fill = TRUE),
           file.path(OUT, "tableS_acta2_control.tsv"), sep = "\t")
}

## F: continuum donor-level trends (Spearman rho of feature vs pseudotime)
nice_feat <- c(vascular_stabilizing = "Vascular-stabilizing",
               synthetic_contractile = "Synthetic/contractile",
               activated_migratory = "Activated/migratory",
               inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
               AGTR1 = "AGTR1")
trend <- fread(P("pericyte_states", "_m", "pseudotime_trend_correlations.tsv")) %>%
    filter(level == "donor") %>%
    mutate(feature = sub("_score$", "", feature),
           feature = recode(feature, AGTR1_expr = "AGTR1"),
           feature = recode(feature, !!!nice_feat),
           sig = ifelse(p_value < 0.05, "p < 0.05", "n.s.")) %>%
    arrange(spearman_rho)
trend$feature <- factor(trend$feature, levels = trend$feature)
pF <- ggplot(trend, aes(spearman_rho, feature, colour = sig)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = spearman_rho, yend = feature), linewidth = 0.5) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("p < 0.05" = "#D55E00", "n.s." = "grey60"), name = NULL) +
    labs(x = "Spearman correlation", y = NULL) +
    theme_ms() + theme(legend.position = c(0.72, 0.18),
                       legend.background = element_blank(),
                       legend.text = element_text(size = 6),
                       plot.margin = margin(3, 9, 3, 3))

## ---- assemble main figure ----------------------------------------------
## Five panels (the by-state balance panel is now in the supplement): the empty
## bottom-right cell keeps pD/pF aligned to the top-row panel widths.
main <- (pA | pB | pC) / (pD | pF | plot_spacer()) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 10))
save_fig("figure_mechanism_main", main, 7.2, 6.0)

## ---- Supplement: program x protein-category enrichment (dot-heatmap) -----
## The science figure the collaborator asked for: how the five pericyte programs
## relate to the eight marker protein categories. Dot size = prevalence (mean
## fraction of the program's cells expressing the category); fill = relative
## enrichment (z of mean detection across programs, within category) so the
## program-specific signal -- otherwise swamped by category baseline detectability
## -- is visible. Cells are assigned to their dominant program by z-scored argmax
## (03b.program_category_enrichment.py), keeping all five programs.
CAT_LEVELS <- c("Signaling ligands", "Chemokines", "Cytokines", "Adhesion molecules",
                "ECM structural", "Matrix-remodeling", "Fibrotic mediators",
                "Mural identity")
pce_f <- P("cell_communication", "_m", "program_category_enrichment.tsv")
gpd_f <- P("cell_communication", "_m", "gene_program_detection.tsv")
if (file.exists(pce_f)) {
    pce <- fread(pce_f) %>%
        group_by(category) %>%
        mutate(enrich = as.numeric(scale(mean_detect))) %>% ungroup() %>%
        mutate(program = factor(program, levels = rev(STATE_LEVELS)),
               category = factor(category, levels = CAT_LEVELS)) %>%
        filter(!is.na(program), !is.na(category))
    pHeat <- ggplot(pce, aes(category, program)) +
        geom_point(aes(size = mean_detect, fill = enrich), shape = 21, colour = "grey30",
                   stroke = 0.2) +
        scale_size_area(max_size = 8, name = "Fraction\nexpressing",
                        breaks = c(0.05, 0.2, 0.4, 0.6)) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "Relative\nenrichment (z)") +
        scale_y_discrete(labels = STATE_LABS1) +
        labs(x = NULL, y = NULL) + theme_ms() +
        theme(axis.text.x = element_text(angle = 35, hjust = 1),
              panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
              legend.position = "right", legend.box = "vertical",
              legend.key.size = unit(3.5, "mm"), legend.title = element_text(size = 6.5),
              legend.text = element_text(size = 6))
    save_fig("figureS_program_category", pHeat, 6.6, 3.0)
}

## ---- Supplement: alluvial  stable cluster -> program -> effector category --
## A single coherent flow: the six stable pericyte subclusters collapse onto their
## dominant program (vascular-stabilizing / fibroblast-like / activated-migratory;
## the discrete state model), and each program's effector output is decomposed into
## six functional molecule classes. Left->middle ribbon width = cluster cell count;
## middle->right width = that program's effector-expression composition (mean
## fraction expressing per class, normalised within program so flow is conserved).
## Effector classes only (no role axis, mural-identity dropped): the classes are
## kept homogeneous per the collaborator's request.
EFF_LEVELS <- c("Signaling ligands", "Chemokines/cytokines", "Adhesion molecules",
                "ECM structural", "Matrix-remodeling", "Fibrotic mediators")
PROG_ORDER <- c("vascular_stabilizing", "fibroblast_like", "activated_migratory")
mf_f  <- P("cell_communication", "_m", "marker_fractions_by_cluster.tsv.gz")
map_f <- P("pericyte_states", "_m", "annotations", "state_program_map.tsv")
cnt_f <- P("pericyte_states", "_m", "state_counts.tsv")
if (all(file.exists(c(mf_f, map_f, cnt_f)))) {
    prog_map <- fread(map_f) %>%
        transmute(cluster = as.character(pericyte_state), program = state_program)
    counts <- fread(cnt_f) %>%
        transmute(cluster = as.character(pericyte_state), n = as.numeric(count))
    ## per cluster x effector-class detection (mean over genes in the class)
    cl_eff <- fread(mf_f) %>%
        filter(grouping == "pericyte_state", category != "Mural identity") %>%
        mutate(eff = recode(category, Chemokines = "Chemokines/cytokines",
                            Cytokines = "Chemokines/cytokines"),
               cluster = as.character(cluster)) %>%
        group_by(cluster, eff) %>% summarise(detect = mean(frac_expr), .groups = "drop") %>%
        left_join(prog_map, by = "cluster") %>% left_join(counts, by = "cluster")
    ## program effector composition (cell-weighted across its clusters, normalised)
    prog_eff <- cl_eff %>% group_by(program, eff) %>%
        summarise(detect = weighted.mean(detect, n), .groups = "drop") %>%
        group_by(program) %>% mutate(share = detect / sum(detect)) %>% ungroup()
    ## cluster ordering: grouped by program (matching the middle axis), then id
    cl_order <- counts %>% left_join(prog_map, by = "cluster") %>%
        mutate(program = factor(program, levels = PROG_ORDER)) %>%
        arrange(program, as.integer(cluster)) %>% pull(cluster)
    ## Colour alluvia by STABLE CLUSTER (left axis) so the left->middle merge is
    ## traceable; cluster hues are grouped into program families (blues = vascular-
    ## stabilizing, oranges = fibroblast-like, pink = activated-migratory) so the
    ## program grouping still reads at the middle axis.
    CLUST_COL <- c(P0 = "#08519C", P2 = "#6BAED6",                 # vascular_stabilizing
                   P1 = "#D94801", P3 = "#FD8D3C", P5 = "#FDD0A2",  # fibroblast_like
                   P4 = "#CC79A7")                                  # activated_migratory
    alluv <- counts %>% left_join(prog_map, by = "cluster") %>%
        left_join(prog_eff %>% select(program, eff, share), by = "program") %>%
        mutate(freq = n * share,
               cluster = factor(paste0("P", cluster), levels = paste0("P", cl_order)),
               program = factor(STATE_LABS[program], levels = STATE_LABS[PROG_ORDER]),
               eff = factor(eff, levels = EFF_LEVELS))
    pAll <- ggplot(alluv, aes(axis1 = cluster, axis2 = program, axis3 = eff, y = freq)) +
        geom_alluvium(aes(fill = cluster), alpha = 0.85, width = 0.28) +
        geom_stratum(width = 0.28, fill = "grey96", colour = "grey55", linewidth = 0.3) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.4,
                  lineheight = 0.85) +
        scale_x_discrete(limits = c("Stable cluster", "Dominant program", "Effector class"),
                         expand = expansion(add = c(0.35, 0.35))) +
        scale_fill_manual(values = CLUST_COL, guide = "none") +
        labs(x = NULL, y = NULL) +
        theme_ms(9) +
        theme(legend.position = "none", panel.border = element_blank(),
              panel.grid = element_blank(),
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_text(face = "bold"))
    save_fig("figureS_alluvial", pAll, 7.2, 5.0)
}

## ---- Supplement: mouse cross-species (continuous, integration-agnostic) --
## The sparse mouse mural set does not support the 5-state clustering, so the
## conserved signal is shown as the continuous relationship between Agtr1a and
## the vascular-stabilizing program -- positive and consistent across all four
## mouse datasets. Agtr1a marks the homeostatic contractile mural compartment,
## matching the wet-lab pericyte-loss / losartan-rescue phenotype.
mouse_cell_f <- P("cross_species", "_m", "mouse_states_metadata.tsv.gz")
mouse_rho_f  <- P("cross_species", "_m", "stats_data", "mouse_Agtr1a_per_dataset_consistency.tsv")
if (file.exists(mouse_cell_f) && file.exists(mouse_rho_f)) {
    mc <- fread(mouse_cell_f)[is_mural == TRUE]
    rho <- fread(mouse_rho_f)
    ## compact, stable dataset labels (D1..Dk) shared across both sub-panels
    ds_levels <- rho[order(-rho), dataset_id]
    ds_labs <- setNames(paste0("D", seq_along(ds_levels)), ds_levels)
    mc[, ds := factor(ds_labs[dataset_id], levels = unname(ds_labs))]
    rho[, ds := factor(ds_labs[dataset_id], levels = unname(ds_labs))]
    ds_col <- setNames(c("#0072B2", "#009E73", "#E69F00", "#CC79A7",
                         "#56B4E9", "#D55E00")[seq_along(ds_levels)], unname(ds_labs))

    pScatter <- ggplot(mc, aes(vascular_stabilizing_score, Agtr1a_expr, colour = ds)) +
        geom_point(size = 0.5, alpha = 0.35) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
        scale_colour_manual(values = ds_col, name = NULL) +
        labs(x = "Vascular-stabilizing score",
             y = expression(italic("Agtr1a") ~ "(scVI-corrected)")) +
        theme_ms() +
        theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
              legend.justification = c(1, 0),
              legend.key.size = unit(7, "pt"), legend.text = element_text(size = 6),
              legend.background = element_blank())

    pRho <- ggplot(rho, aes(ds, rho, fill = ds)) +
        geom_col(width = 0.7) +
        geom_hline(yintercept = 0, linewidth = 0.3) +
        geom_text(aes(label = paste0("n=", n)), vjust = -0.4, size = 1.9) +
        scale_fill_manual(values = ds_col) +
        scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
        labs(x = NULL, y = expression("Spearman " * rho * " (Agtr1a vs stabilizing)")) +
        theme_ms() + theme(axis.title.y = element_text(size = 6.5))

    pMouse <- (pScatter | pRho) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(face = "bold", size = 10))
    save_fig("figureS_crossspecies_mouse", pMouse, 6.4, 3.0)
}

## ===== Figure: CCC + NicheNet into pericytes ============================
RECV <- "Pericytes"
la_f <- P("cell_communication", "_m", "nichenet", paste0("ligand_activities_", RECV, ".tsv"))
lt_f <- P("cell_communication", "_m", "nichenet", paste0("ligand_target_links_", RECV, ".tsv"))
fr_f <- P("cell_communication", "_m", "expressed_fraction_per_group.tsv.gz")
if (all(file.exists(c(la_f, lt_f, fr_f)))) {
    TGF <- c("TGFB1", "TGFB2", "TGFB3", "CCN1", "CCN2")
    la <- fread(la_f) %>% slice_max(aupr_corrected, n = 15)
    la$test_ligand <- factor(la$test_ligand, levels = rev(la$test_ligand))
    keep_lig <- rev(levels(la$test_ligand))

    ## cA: NicheNet ligand activity (TGF-beta family highlighted)
    cA <- ggplot(la, aes(aupr_corrected, test_ligand,
                         colour = test_ligand %in% TGF)) +
        geom_segment(aes(x = 0, xend = aupr_corrected, yend = test_ligand), linewidth = 0.5) +
        geom_point(size = 1.9) +
        scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
        labs(x = "Ligand activity (AUPR corrected)", y = NULL) + theme_ms()

    ## cB: ligand -> target regulatory potential heatmap
    lt <- fread(lt_f) %>% filter(ligand %in% keep_lig)
    top_tgt <- lt %>% group_by(target) %>% summarise(w = sum(weight), .groups = "drop") %>%
        slice_max(w, n = 24) %>% arrange(desc(w)) %>% pull(target)
    lt <- lt %>% filter(target %in% top_tgt) %>%
        mutate(ligand = factor(ligand, levels = keep_lig),
               target = factor(target, levels = top_tgt))
    cB <- ggplot(lt, aes(target, ligand, fill = weight)) +
        geom_tile(colour = "white", linewidth = 0.2) +
        scale_fill_gradient(low = "#FDE0DD", high = "#7A0177", name = "Regulatory\npotential") +
        labs(x = "Predicted target gene", y = NULL) + theme_ms() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5.5),
              panel.grid = element_blank(), legend.position = "right",
              legend.key.size = unit(3, "mm"), legend.title = element_text(size = 6))

    ## cC: which niche cells express the prioritized ligands
    frac <- fread(fr_f); setnames(frac, 1, "gene")
    fl <- frac %>% filter(gene %in% keep_lig) %>%
        pivot_longer(-gene, names_to = "group", values_to = "frac") %>%
        filter(!group %in% c(RECV, "Pericyte_AGTR1neg"))
    ord <- fl %>% group_by(group) %>% summarise(s = sum(frac), .groups = "drop") %>%
        slice_max(s, n = 12) %>% arrange(desc(s)) %>% pull(group)
    fl <- fl %>% filter(group %in% ord) %>%
        mutate(group = factor(group, levels = ord),
               gene = factor(gene, levels = keep_lig))
    cC <- ggplot(fl, aes(group, gene, size = frac, colour = frac)) +
        geom_point() +
        scale_size_area(max_size = 3.2, name = "Fraction\nexpressing") +
        scale_colour_viridis_c(option = "magma", direction = -1, guide = "none") +
        labs(x = NULL, y = NULL) + theme_ms() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
              legend.position = "right", legend.key.size = unit(3, "mm"),
              legend.title = element_text(size = 6))

    ## cD: projectR transfer -- pericyte-learned CoGAPS patterns projected onto
    ## the niche (cell-type x donor pseudobulk). Shows the fibroblast-like injury
    ## program is mirrored most strongly in bona-fide fibroblasts, while the
    ## stabilizing / contractile programs stay mural-compartment-specific.
    pp_f  <- P("pericyte_cogaps", "_m", "projected_pattern_by_celltype_np5.tsv")
    ann_f <- P("cell_communication", "_m", "cogaps_receiver_annotation_np5.tsv")
    cD <- NULL
    if (file.exists(pp_f)) {
        pp <- fread(pp_f)
        nice <- c(vascular_stabilizing = "Vascular-stabilizing",
                  synthetic_contractile = "Synthetic/contractile",
                  inflammatory = "Inflammatory", fibroblast_like = "Fibroblast-like",
                  activated_migratory = "Activated/migratory")
        prog <- c(Pattern_4 = "Vascular-stabilizing", Pattern_5 = "Synthetic/contractile",
                  Pattern_1 = "Inflammatory", Pattern_3 = "Fibroblast-like")
        if (file.exists(ann_f)) {
            a <- fread(ann_f)
            m <- nice[a$assigned_program]; names(m) <- a$pattern
            m <- m[!is.na(m)]                       # drop unassigned/mixed patterns
            if (length(m)) prog <- m
        }
        pats <- names(prog)
        mat <- as.matrix(pp[, ..pats]); rownames(mat) <- pp$cell_type
        z <- scale(mat)                             # z within pattern (across cell types)
        prog_levels <- intersect(c("Vascular-stabilizing", "Synthetic/contractile",
                                    "Inflammatory", "Fibroblast-like",
                                    "Activated/migratory"), prog)
        fl_pat <- names(prog)[prog == "Fibroblast-like"][1]
        ct_order <- rownames(z)[order(z[, fl_pat], decreasing = TRUE)]
        long <- as.data.table(as.table(z)); setnames(long, c("cell_type", "pattern", "z"))
        long[, program := factor(prog[as.character(pattern)], levels = rev(prog_levels))]
        long[, cell_type := factor(cell_type, levels = ct_order)]
        cD <- ggplot(long, aes(cell_type, program, fill = z)) +
            geom_tile(colour = "white", linewidth = 0.2) +
            scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                 midpoint = 0, name = "z (within\nprogram)") +
            labs(x = NULL, y = NULL) + theme_ms() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
                  axis.text.y = element_text(size = 6.5), panel.grid = element_blank(),
                  legend.position = "right", legend.key.size = unit(3, "mm"),
                  legend.title = element_text(size = 6))
    }

    figA <- (cA + cB + plot_layout(widths = c(1, 2.1))) / cC
    if (!is.null(cD)) figA <- figA / cD
    figA <- figA +
        plot_layout(heights = if (!is.null(cD)) c(1, 0.9, 0.75) else c(1, 0.9)) +
        plot_annotation(tag_levels = "A") &
        theme(plot.tag = element_text(face = "bold", size = 10))
    save_fig("figure_ccc_nichenet", figA, 9.2, if (!is.null(cD)) 9.2 else 7.0)
}

cat("Wrote manuscript figures to", OUT, "\n")
list.files(OUT, pattern = "\\.(pdf|svg)$")
cat("\nReproducibility information:\n"); sessioninfo::session_info()
