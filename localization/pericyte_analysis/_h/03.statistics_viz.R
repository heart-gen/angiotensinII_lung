## Analyze pericyte subclusters
## Normality checked with shapiro (AGTR1 non-normal; ACTA2 okay)
## Heteroscedasticity okay for both genes

suppressPackageStartupMessages({
    library(dplyr)
    library(ggpubr)
    library(emmeans)
})

load_data <- function(min_donors = 2){
    df <- data.table::fread("pericytes_metadata.tsv.gz")
    donor_df <- group_by(df, donor_id) |>
        summarise(
            AGTR1_mean = mean(AGTR1_expr, na.rm=TRUE),
            ACTA2_mean = mean(ACTA2_expr, na.rm=TRUE),
            sex = first(sex),
            disease = first(disease),
            ethnicity = first(self_reported_ethnicity),
            age = mean(age_or_mean_of_age_range, na.rm=TRUE),
            leiden_pericytes = names(which.max(table(leiden_pericytes)))
        ) |>
        tidyr::drop_na() |>
        mutate(
            leiden_pericytes = factor(leiden_pericytes),
            sex = factor(sex),
            disease = factor(disease),
            ethnicity = factor(ethnicity),
            donor_id = factor(donor_id),
        ) |>
        group_by(leiden_pericytes) |>
        filter(n() >= min_donors) |>
        ungroup() |>
        mutate(
            leiden_pericytes = droplevels(leiden_pericytes),
            sex = droplevels(sex),
            disease = droplevels(disease),
            ethnicity = droplevels(ethnicity)
        )
    return(donor_df)
}

save_lm_diagnostics <- function(fit, file, width = 8, height = 8) {
    pdf(file, width = width, height = height)
    par(mfrow = c(2, 2))
    plot(fit)
    dev.off()
}

run_ancova <- function(df, response, plot_dir) {
    formula <- reformulate(
        c("leiden_pericytes", "sex", "disease", "ethnicity", "age"),
        response
    )

                                        # Fit model
    fit <- lm(formula, data = df)
    robust_test <- lmtest::coeftest(
        fit, vcov = sandwich::vcovHC(fit, type = "HC3")
    )

                                        # Plot QC
    if(!dir.exists(plot_dir)) { dir.create(plot_dir) }
    
    plot_file <- file.path(plot_dir, paste0(response, "_lm_diagnostics.pdf"))
    save_lm_diagnostics(fit, plot_file)

    list(
        model = fit,
        robust_coefs = robust_test,
        anova = car::Anova(fit, type = 2),
        emmeans = emmeans(fit, ~ leiden_pericytes, weights="proportional")
    )
}

save_tables <- function(res, outdir, gene) {
    if(!dir.exists(outdir)) { dir.create(outdir) }
    write.table(
        as.data.frame(res$anova),
        file = file.path(outdir, paste0(gene, "_anova.tsv")),
        sep = "\t", quote = FALSE
    )

    write.table(
        as.data.frame(res$emmeans),
        file = file.path(outdir, paste0(gene, "_emmeans.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
        as.data.frame(pairs(res$emmeans)),
        file = file.path(outdir, paste0(gene, "_posthoc.tsv")),
        sep = "\t", quote = FALSE, row.names = FALSE
    )

    write.table(
        as.data.frame(res$robust_coefs),
        file = file.path(outdir, paste0(gene, "_robust_coefs.tsv")),
        sep = "\t", quote = FALSE, row.names = TRUE
    )
}

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

boxplot_subcluster <- function(df, palette, label, gene_name, filename) {
    my_comp <- list(c("0", "1"), c("1", "2"), c("2", "3"), c("0", "2"),
                    c("1", "3"), c("0", "3"))
    bxp <- ggboxplot(
        df, x = "leiden_pericytes", y = label, add = "jitter",
        add.params = list(alpha=0.5, size=1),
        fill = "leiden_pericytes", palette = palette,
        xlab = "Pericyte Subclusters",
        ylab = paste0("Normalized Expression (", gene_name, ")"),
        legend = "none", ggtheme = theme_pubr(base_size = 15)
    ) +
        stat_summary(
            fun = mean, geom = "point", shape = 23, size = 3,
            fill = "white", color = "black"
        ) +
        stat_compare_means(comparisons = my_comp)
    save_ggplots(filename, bxp, 4, 5)
}

### Main script section
df <- load_data()
print(dim(df))

                                        # Statistics
outdir  <- "stats_data"
if(!dir.exists(outdir)) { dir.create(outdir) }

res_agtr1 <- run_ancova(df, "AGTR1_mean", outdir)
res_acta2 <- run_ancova(df, "ACTA2_mean", outdir)

save_tables(res_agtr1, outdir, "AGTR1")
save_tables(res_acta2, outdir, "ACTA2")

tab20 <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
  "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
  "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
  "#17becf", "#9edae5" 
)
boxplot_subcluster(df, tab20, "AGTR1_mean", "AGTR1",
                   file.path(outdir, "boxplot_agtr1.norm_expression"))
boxplot_subcluster(df, tab20, "ACTA2_mean", "ACTA2",
                   file.path(outdir, "boxplot_acta2.norm_expression"))

                                        # Gene comparisons
corr_by_cluster <- df |>
    group_by(leiden_pericytes) |>
    summarise(
        n_donors = n(),
        spearman_rho = cor(
            AGTR1_mean, ACTA2_mean, method = "spearman",
            use = "pairwise.complete.obs"
        ),
        p_value = cor.test(
            AGTR1_mean, ACTA2_mean, method = "spearman"
        )$p.value,
        .groups = "drop"
  )
corr_by_cluster |> as.data.frame() |>
    write.table(
        file = file.path(outdir, "correlation_by_subcluster.tsv"),
        sep = "\t", quote = FALSE, row.names = TRUE
    )

sca <- ggscatter(df, x = "AGTR1_mean", y = "ACTA2_mean", add = "reg.line",
                 add.params = list( color = "blue", fill = "lightgray"),
                 facet.by = "leiden_pericytes", conf.int = TRUE, cor.coef = TRUE,
                 cor.coeff.args = list(method = "spearman", label.sep="\n",
                                       label.x=1.5, label.y=4),
                 xlab = "Normalized Expression (AGTR1)",
                 ylab = "Normalized Expression (ACTA2)", ncol=4,
                 legend = "none", ggtheme = theme_pubr(base_size = 15))
save_ggplots(file.path(outdir, "correlation_by_subcluster"), sca, 12, 4)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
