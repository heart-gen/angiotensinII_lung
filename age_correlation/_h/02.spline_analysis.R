suppressPackageStartupMessages({
    library(mgcv)
    library(dplyr)
    library(splines)
})

prepare_age_agtr1_donor_table <- function(outdir = "mean_expr") {
    return(data.table::fread(file.path(outdir, "donor_metadata.tsv")))
}

get_agtr1_enriched_celltypes <- function(
    donor_celltype, metric = c("mean_expr", "frac_expr"), top_n = 5) {
    metric <- match.arg(metric)
    summ <- donor_celltype |>
        group_by(cell_type) |>
        summarise(
            n_donors = n(), mean_expr = mean(AGTR1_mean, na.rm = TRUE),
            frac_expr = mean(frac_AGTR1_pos, na.rm = TRUE),
            .groups = "drop"
        ) |>
        arrange(desc(.data[[metric]])) |>
        slice_head(n = top_n)
    return(summ)
}

fit_age_spline_tests <- function(donor_celltype, keep_celltypes, df_spline = 3) {
    donor_celltype |>
        filter(cell_type %in% keep_celltypes) |>
        group_by(cell_type) |>
        group_modify(~{
            dat <- .x |> tidyr::drop_na(age, AGTR1_mean)
            if (nrow(dat) < (df_spline + 5)) return(tibble())

            m_lin <- lm(AGTR1_mean ~ age + sex + ethnicity, data = dat)
            m_spl <- lm(AGTR1_mean ~ ns(age, df = df_spline) + sex +
                      ethnicity, data = dat)

            a <- anova(m_lin, m_spl)  # nested model test
            tibble(
                n_donors = nrow(dat), df_spline = df_spline,
                p_spline_vs_linear = a$`Pr(>F)`[2],
                AIC_linear = AIC(m_lin), AIC_spline = AIC(m_spl)
            )
        }) |>
    ungroup() |>
    mutate(fdr = p.adjust(p_spline_vs_linear, method = "fdr")) |>
    arrange(fdr, p_spline_vs_linear)
}

fit_gam_by_celltype <- function(donor_celltype, keep_celltypes, k = 3) {
    donor_celltype |>
        filter(cell_type %in% keep_celltypes) |>
        group_by(cell_type) |>
        group_modify(~{
            dat <- .x |> tidyr::drop_na(age, AGTR1_mean)
            if (nrow(dat) < (k + 5)) return(tibble())

            m_gam <- gam(AGTR1_mean ~ s(age, k = k) + sex + ethnicity,
                         data = dat, method = "REML")
            sm <- summary(m_gam)

            tibble(
                n_donors = nrow(dat), edf = sm$s.table[1, "edf"],
                p_smooth = sm$s.table[1, "p-value"], r2 = sm$r.sq
            )
        }) |>
    ungroup() |>
    mutate(fdr = p.adjust(p_smooth, method = "fdr")) |>
    arrange(fdr, p_smooth)
}

age_agtr1_analysis <- function(
    top_n_celltypes = 5, enrichment_metric = c("mean_expr", "frac_expr"),
    dof = 3
) {
    enrichment_metric <- match.arg(enrichment_metric)
    outdir <- enrichment_metric

                                        # Load donor x cell_type table for AGTR1
    donor_celltype <- prepare_age_agtr1_donor_table(enrichment_metric)

                                        # Select enriched cell types descriptively
    enriched <- get_agtr1_enriched_celltypes(
        donor_celltype, metric = enrichment_metric, top_n = top_n_celltypes
    )
    keep_celltypes <- enriched$cell_type

                                        # Spline and GAM
    stats_spline <- fit_age_spline_tests(donor_celltype, keep_celltypes, dof)
    stats_gam    <- fit_gam_by_celltype(donor_celltype, keep_celltypes, dof)

                                        # Save outputs
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(stats_spline, file.path(outdir, "agtr1_spline_by_celltypes.tsv"),
                       sep = "\t")
    data.table::fwrite(stats_gam, file.path(outdir, "age_agtr1_gam_by_celltype.tsv"),
                       sep = "\t")
    return(NULL)
}

#### Main
                                        # Top 5 cell types by mean AGTR1 expression
age_agtr1_analysis(
    enrichment_metric = "mean_expr", top_n_celltypes = 5, dof = 3
)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
