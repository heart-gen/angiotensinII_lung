## Shared plumbing for every supplementary-table script in tables/_h/.
##
## Source from the top of a table script, AFTER the library() block. The scripts
## here are run with `Rscript ../_h/<name>.R` from tables/_m:
##     source("../_h/_tab_common.R")
##
## Mirrors figures/_h/_fig_common.R: same P() path idiom, same disease vocabulary.
## Provides ROOT/P/OUT/TSVDIR, write_part(), and the numbering registry.

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

ROOT   <- normalizePath(file.path(getwd(), "..", ".."))
P      <- function(...) file.path(ROOT, ...)
OUT    <- P("tables", "_m")
TSVDIR <- file.path(OUT, "tsv")
dir.create(TSVDIR, showWarnings = FALSE, recursive = TRUE)

## Disease vocabulary -- keep identical to figures/_h/_fig_common.R so a table and
## its figure never disagree on a group label.
DISEASE_LEVELS <- c("Healthy", "COPD", "Fibrotic_ILD", "Other")
DISEASE_LABS   <- c(Healthy = "Healthy", COPD = "COPD",
                    Fibrotic_ILD = "Fibrotic/ILD", Other = "Other")

## ---- numbering registry -------------------------------------------------
## write_part() appends here; 08.assemble_tables.R reads it to build the workbook.
## This file is the single authority for supplementary-TABLE numbering, exactly as
## the tribble in figures/_h/assemble_mechanism_figures.R is for FIGURE numbering.
## Table scripts must never hard-code a sheet order.
REGISTRY <- file.path(OUT, "manifest_parts.tsv")

## `status`:
##   complete         -- every requested column is populated from a source file
##   pending_upstream -- an upstream re-run has not happened yet, so the part is
##                       written with what exists and flagged. 08 refuses to build
##                       a final workbook while any part is pending.
write_part <- function(df, part, title, supports = NA_character_,
                       sources = character(), status = "complete",
                       notes = NA_character_) {
    stopifnot(is.data.frame(df), nchar(part) > 0)
    df <- as.data.frame(df, check.names = FALSE)
    f  <- file.path(TSVDIR, paste0("tableS", part, ".tsv"))
    write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

    row <- data.frame(
        part = part, title = title, supports = supports,
        n_rows = nrow(df), n_cols = ncol(df), status = status,
        notes = notes, file = basename(f),
        sources = paste(sources, collapse = "; "),
        built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE)

    ## Re-registering a part replaces its row rather than duplicating it, so a
    ## script can be re-run standalone without corrupting the manifest.
    if (file.exists(REGISTRY)) {
        ## as.data.frame FIRST, deliberately. Inside `[.data.table` the bare `part`
        ## on the right-hand side resolves to the COLUMN, not this function's
        ## argument, so `old[old$part != part]` compares the column with itself,
        ## matches nothing, and silently discards every previously registered part.
        old <- as.data.frame(fread(REGISTRY, colClasses = "character"))
        old <- old[old$part != part, , drop = FALSE]
        row <- rbind(old, row)
    }
    row <- row[order(nchar(row$part), row$part), ]
    write.table(row, REGISTRY, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    cat(sprintf("  [S%s] %-58s %4d x %2d  %s\n", part, substr(title, 1, 58),
                nrow(df), ncol(df), status))
    invisible(f)
}

## ---- helpers ------------------------------------------------------------
## Several source files carry p_value == 0 from underflow. Printing "0" in a
## supplementary table is a claim no test can support, so floor it.
fmt_p <- function(p, floor_at = 2.2e-16) {
    ifelse(is.na(p), NA_character_,
           ifelse(p < floor_at, paste0("<", format(floor_at, scientific = TRUE)),
                  format(p, digits = 3, scientific = TRUE)))
}

bh <- function(p) p.adjust(p, method = "BH")

## Fisher z 95% CI for a Pearson/Spearman correlation. Used where a source file
## stored a bare coefficient with no interval (state_gate_axis_correlation.tsv).
fisher_ci <- function(r, n, conf = 0.95) {
    z  <- atanh(pmin(pmax(r, -0.999999), 0.999999))
    se <- 1 / sqrt(pmax(n - 3, 1))
    k  <- qnorm(1 - (1 - conf) / 2)
    list(lo = tanh(z - k * se), hi = tanh(z + k * se))
}

## Read a source file, or return NULL with a warning. Table scripts use this so a
## single missing upstream output degrades one part instead of killing the run.
read_src <- function(path, ...) {
    if (!file.exists(path)) { warning("MISSING source: ", path, call. = FALSE); return(NULL) }
    fread(path, ...)
}

have <- function(...) all(vapply(list(...), function(x) !is.null(x) && nrow(x) > 0, logical(1)))
