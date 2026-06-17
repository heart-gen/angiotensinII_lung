## Alluvial plots: pericyte cluster/state -> protein category -> candidate role.
##
## Ribbon width = mean fraction of cells in a cluster/state expressing the
## markers of each category (from 03a.marker_fractions.py). Visualizes which
## pericyte populations carry signaling-ligand / chemokine / cytokine /
## adhesion / ECM / matrix-remodeling / fibrotic programs and their roles.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ggalluvial)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
INPUT  <- parse_arg("--input", "../_m/marker_fractions_by_cluster.tsv.gz")
OUTDIR <- parse_arg("--outdir", "../_m/figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

df <- data.table::fread(INPUT)

save_ggplots <- function(fn, p, w, h) {
    for (ext in c(".pdf", ".png")) ggsave(paste0(fn, ext), plot = p, width = w, height = h)
}

make_alluvial <- function(d, grouping, title, fn) {
    flows <- d |>
        filter(grouping == !!grouping) |>
        group_by(cluster, category, role) |>
        summarise(freq = mean(frac_expr, na.rm = TRUE), .groups = "drop") |>
        filter(freq > 0)

    p <- ggplot(flows,
                aes(axis1 = cluster, axis2 = category, axis3 = role, y = freq)) +
        geom_alluvium(aes(fill = category), alpha = 0.75, width = 1/4) +
        geom_stratum(width = 1/4, fill = "grey95", color = "grey40") +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                  size = 2.8) +
        scale_x_discrete(limits = c("Cluster", "Category", "Role"),
                         expand = c(0.10, 0.10)) +
        labs(title = title, y = "Mean fraction expressing", fill = "Category") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom",
              panel.grid.major.x = element_blank())
    save_ggplots(fn, p, 12, 8)
}

if ("leiden_pericytes" %in% unique(df$grouping)) {
    make_alluvial(df, "leiden_pericytes",
                  "Pericyte subcluster -> protein category -> role",
                  file.path(OUTDIR, "alluvial_leiden"))
}
if ("pericyte_state" %in% unique(df$grouping)) {
    make_alluvial(df, "pericyte_state",
                  "Pericyte state -> protein category -> role",
                  file.path(OUTDIR, "alluvial_state"))
}

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
