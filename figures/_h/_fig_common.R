## Shared visual language for every script in figures/_h/.
##
## Source this from the top of a figure script, AFTER the library() block:
##     source(file.path(dirname(sys.frame(1)$ofile %||% "."), "_fig_common.R"))
## or, as the scripts here are run with `Rscript ../_h/<name>.R` from figures/_m:
##     source("../_h/_fig_common.R")
##
## Provides ROOT/P/OUT, save_fig(), the Okabe-Ito palette, the disease vocabulary,
## and a parameterised theme_ms().
##
## Deliberately NOT shared: STATE_LEVELS / STATE_LABS / PROG_* / CLUST_*. Those
## legitimately differ between scripts -- e.g. manuscript_mechanism_figure.R orders
## STATE_LEVELS with basement_membrane before fibroblast_like while
## pericyte_layer_figure.R does the reverse, and that ordering drives plotting
## order. Unifying them would silently reorder published panels.

ROOT <- normalizePath(file.path(getwd(), "..", ".."))
P    <- function(...) file.path(ROOT, ...)
OUT  <- P("figures", "mechanism")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## geom_jitter draws its offsets at RENDER time, so without a seed the same script
## produces a visibly different figure on every run -- a published panel could not
## be regenerated. Seed once here, for every figure script.
set.seed(20260727)

## ---- palette ------------------------------------------------------------
OKABE <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9",
           "#F0E442", "#999999")

DISEASE_LEVELS <- c("Healthy", "COPD", "Fibrotic_ILD", "Other")
DISEASE_LABS   <- c(Healthy = "Healthy", COPD = "COPD",
                    Fibrotic_ILD = "Fibrotic/ILD", Other = "Other")
DISEASE_COL    <- c(Healthy = "#0072B2", COPD = "#E69F00",
                    Fibrotic_ILD = "#D55E00", Other = "#999999")

## ---- theme --------------------------------------------------------------
## `.theme_ms` is the parameterised base; `theme_ms` is the majority variant
## (pericyte_layer / disease_main / sensitivity) and is what a script gets for free.
## The two scripts that deviate re-bind theme_ms after sourcing this file:
##   manuscript_mechanism_figure.R  theme_ms <- \(base = 8) .theme_ms(base, grid_x = FALSE)
##   basement_membrane_figure.R     theme_ms <- \(base = 8) .theme_ms(base, legend = NULL)
.theme_ms <- function(base = 8, grid_x = TRUE, legend = "none", strip_size = NULL) {
    th <- theme_bw(base_size = base) +
        theme(plot.title = element_blank(), plot.subtitle = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text = element_text(colour = "black"),
              axis.title = element_text(colour = "black"),
              strip.background = element_blank())
    if (!grid_x) th <- th + theme(panel.grid.major.x = element_blank())
    if (!is.null(legend)) th <- th + theme(legend.position = legend)
    if (!is.null(strip_size))
        th <- th + theme(strip.text = element_text(face = "bold", size = strip_size))
    th
}
theme_ms <- function(base = 8) .theme_ms(base = base)

## ---- export -------------------------------------------------------------
## Every supplement/main figure goes out as vector PDF (cairo, for the journal),
## SVG (for vector-editor assembly; gitignored) and 350 dpi PNG (for quick review).
save_fig <- function(fn, p, w, h) {
    ggsave(file.path(OUT, paste0(fn, ".pdf")), p, width = w, height = h, device = cairo_pdf)
    ggsave(file.path(OUT, paste0(fn, ".svg")), p, width = w, height = h)
    ggsave(file.path(OUT, paste0(fn, ".png")), p, width = w, height = h, dpi = 350)
}

## Panel tag applied to a standalone plot (patchwork's plot_annotation handles the
## assembled case).
tag <- function(p, lab) p + labs(tag = lab) +
    theme(plot.tag = element_text(face = "bold", size = 10))

## ---- shared helpers -----------------------------------------------------
map_disease <- function(lc) {
    lc <- as.character(lc)
    dplyr::case_when(
        grepl("^Healthy", lc) ~ "Healthy", lc %in% c("COPD") ~ "COPD",
        grepl("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis", lc,
              ignore.case = TRUE) ~ "Fibrotic_ILD", TRUE ~ "Other")
}
dx_factor <- function(x) droplevels(factor(map_disease(x), levels = DISEASE_LEVELS))

## ---- stale-label normalisation ------------------------------------------
## The 2026-07-21 gate (basement_membrane/_h/01.state_gate.py) relabelled stable
## clusters 1/3/5 -- 4,200 cells, 35.96% of pericytes, i.e. the ENTIRE former
## `fibroblast_like` set -- to `basement_membrane`. On any table keyed by
## state_program the rename is therefore one-to-one and leaves the estimates
## untouched. Tables regenerated after that date already carry the new level;
## a table that predates it carries the old one and, because every figure script
## filters state_program against a PROG_* vector containing `basement_membrane`,
## silently DROPS that program instead of erroring. That is exactly how panel D
## of figure_pericyte_layer shipped plotting 2 of 3 programs (and, because the
## panel centres within lens, with every remaining point shifted too).
##
## Call this on any state_program-keyed table before filtering. It is a stopgap
## for tables on disk, not a substitute for re-running the producing script --
## if both levels are ever present the 1:1 mapping no longer holds and we stop.
bm_relabel <- function(dt, col = "state_program", src = "table") {
    lv <- unique(as.character(dt[[col]]))
    if (!"fibroblast_like" %in% lv) return(dt)
    if ("basement_membrane" %in% lv)
        stop(src, " carries BOTH fibroblast_like and basement_membrane in `", col,
             "`; the 1:1 relabel no longer holds -- regenerate it rather than merging ",
             "two distinct programs.")
    warning(src, " predates the 2026-07-21 basement-membrane relabel; mapping ",
            "fibroblast_like -> basement_membrane. Re-run the producing script.",
            call. = FALSE)
    dt[[col]] <- sub("^fibroblast_like$", "basement_membrane", as.character(dt[[col]]))
    dt
}

## emmeans names its interval columns by the degrees of freedom it used:
## `lower.CL`/`upper.CL` for a t-based interval, `asymp.LCL`/`asymp.UCL` when df
## is infinite (asymptotic). Two tables written by different fits therefore carry
## different column names for the same quantity, and rbind()-ing them fails with
## "Column 5 ['lower.CL'] of item 2 is missing in item 1" -- which is what stopped
## figureS_acta2_control panel A. Normalise before combining; the fits themselves
## are unchanged and the panel uses only emmean and SE.
norm_ci <- function(dt) {
    setnames(dt, c("asymp.LCL", "asymp.UCL"), c("lower.CL", "upper.CL"),
             skip_absent = TRUE)
    dt
}

## Guard the other half of the same failure: a program that vanished during
## filtering must never be plotted around silently.
require_programs <- function(present, expected, panel) {
    missing <- setdiff(expected, unique(as.character(present)))
    if (length(missing))
        stop(panel, " is missing program(s): ", paste(missing, collapse = ", "),
             ". Refusing to plot a partial panel.")
    invisible(TRUE)
}

fmt_p <- function(p) if (p < 1e-3) sprintf("P = %.1e", p) else sprintf("P = %.3f", p)

## Spearman rho annotation, e.g. "rho = 0.54, P = 4.4e-09"
fmt_rho <- function(rho, p) sprintf("rho = %.2f, %s", rho, tolower(fmt_p(p)))
