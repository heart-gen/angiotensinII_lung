## Assemble the two main mechanistic figures from module outputs.
##
## Figure A (CCC + NicheNet + alluvial): liana edges into AGTR1+/AGTR2+
##   receivers, NicheNet ligand->target heatmap, cluster->category->role alluvial.
## Figure B (states + continuum + niche index): state UMAP, state-by-disease,
##   DPT continuum trend, donor-level niche index by disease.
##
## Run AFTER the production module jobs (cell_communication, pericyte_states,
## niche_index) have written their _m outputs. Panels are combined from the
## per-module PDFs/PNGs; this script stitches the statistical panels it can
## rebuild directly and records the panel manifest for the rest.

suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P <- function(...) file.path(ROOT, ...)
outdir <- file.path(ROOT, "figures", "mechanism"); dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_if <- function(path) if (file.exists(path)) data.table::fread(path) else NULL

## ---- Figure B statistical panels (rebuildable from light tables) --------
niche <- read_if(P("niche_index", "_m", "stats_data", "niche_index_emmeans.tsv"))
injury <- read_if(P("pericyte_states", "_m", "stats_data", "injury_fraction_emmeans.tsv"))
trend <- read_if(P("pericyte_states", "_m", "pseudotime_trend_correlations.tsv"))

panels <- list()
if (!is.null(niche)) {
    panels$niche <- ggplot(niche, aes(disease_group, emmean)) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL)) +
        labs(x = "", y = "Niche-stability index", title = "Niche index by disease") +
        theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
if (!is.null(injury)) {
    panels$injury <- ggplot(injury, aes(disease_group, emmean)) +
        geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), color = "#B2182B") +
        labs(x = "", y = "Injury-state fraction", title = "Injury states by disease") +
        theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}
if (!is.null(trend)) {
    td <- trend |> filter(level == "donor")
    panels$trend <- ggplot(td, aes(reorder(feature, spearman_rho), spearman_rho)) +
        geom_col(fill = "#2166AC") + coord_flip() +
        labs(x = "", y = "Spearman rho (donor)", title = "Continuum trends") +
        theme_bw(base_size = 11)
}

if (length(panels)) {
    figB <- wrap_plots(panels, ncol = length(panels)) +
        plot_annotation(title = "States, continuum, and niche-stability index",
                        tag_levels = "A")
    ggsave(file.path(outdir, "figureB_states_continuum_niche.pdf"), figB,
           width = 5 * length(panels), height = 4.5)
    ggsave(file.path(outdir, "figureB_states_continuum_niche.png"), figB,
           width = 5 * length(panels), height = 4.5, dpi = 300)
}

## ---- Panel manifest (image panels assembled in vector editor) -----------
manifest <- tibble::tribble(
    ~figure, ~panel, ~source,
    "A", "liana edges -> pericytes",        P("cell_communication","_m","figures","dotplot_into_Pericytes.pdf"),
    "A", "liana edges -> AGTR2-det AT2",    P("cell_communication","_m","figures","dotplot_into_AT2_AGTR2det.pdf"),
    "A", "NicheNet ligand->target (pericyte)", P("cell_communication","_m","nichenet","ligand_target_heatmap_Pericytes.pdf"),
    "A", "alluvial cluster->category->role", P("cell_communication","_m","figures","alluvial_state.pdf"),
    "B", "state UMAP",                      P("pericyte_states","_m","figures","umap_pericyte_state.pdf"),
    "B", "DPT pseudotime UMAP",             P("pericyte_states","_m","figures","umap_pseudotime.pdf"),
    "B", "PAGA states",                     P("pericyte_states","_m","figures","paga_states.pdf"),
    "Supp", "AT1R/AT2R balance by disease", P("pathway_balance","_m","stats_data","balance_by_disease.pdf"),
    "Supp", "smoking + cohort robustness", P("figures","mechanism","figureS_sensitivity.pdf"),
    "Supp", "mouse Agtr1a by state",        P("cross_species","_m","stats_data","mouse_Agtr1a_by_state.pdf")
)
manifest$exists <- file.exists(manifest$source)
write.table(manifest, file.path(outdir, "figure_panel_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Panels present:", sum(manifest$exists), "/", nrow(manifest), "\n")
print(manifest[, c("figure", "panel", "exists")])
