## Assemble the supplementary-table workbook and verify it against the manuscript.
##
## Two jobs:
##   1. ANCHOR CHECKS -- a set of values quoted in the manuscript are re-read from
##      the built tables and compared. A mismatch fails the job. A supplementary
##      table that silently disagrees with the text is worse than no table, so this
##      is a hard gate rather than a warning.
##   2. WORKBOOK -- every registered part becomes one worksheet of
##      supplementary_tables.xlsx, ordered by the manifest, with a README sheet.
##
## The manifest (tables/_m/manifest_parts.tsv, written by write_part()) is the
## single authority for table numbering, mirroring the role the tribble in
## figures/_h/assemble_mechanism_figures.R plays for figure numbering.

suppressPackageStartupMessages({
    library(data.table); library(dplyr); library(writexl)
})
source("../_h/_tab_common.R")

if (!file.exists(REGISTRY)) stop("no manifest at ", REGISTRY, " -- run the table scripts first")
man <- fread(REGISTRY, colClasses = "character")
man[, n_rows := as.integer(n_rows)][, n_cols := as.integer(n_cols)]

## Order parts as S1, S2A, S2B, ... S14C (numeric part first, then letter suffix).
man[, tbl_num := as.integer(sub("^([0-9]+).*$", "\\1", part))]
man[, tbl_sfx := sub("^[0-9]+", "", part)]
setorder(man, tbl_num, tbl_sfx)

## =========================================================================
## Anchor checks
## =========================================================================
tsv <- function(part) {
    f <- file.path(TSVDIR, paste0("tableS", part, ".tsv"))
    if (!file.exists(f)) return(NULL)
    fread(f)
}
CHECKS <- list()
chk <- function(label, got, want, tol = 1e-4, threshold_dependent = FALSE) {
    ok <- length(got) == 1 && !is.na(got) && abs(got - want) <= tol
    CHECKS[[length(CHECKS) + 1L]] <<- data.table(
        check = label, expected = want, observed = if (length(got) == 1) got else NA_real_,
        passed = ok, threshold_dependent = threshold_dependent)
    invisible(ok)
}
maybe <- function(x, expr) if (is.null(x)) NA_real_ else tryCatch(expr, error = function(e) NA_real_)

## -- values that do not depend on the donor cell-count threshold ----------
x <- tsv("10A")
chk("S10A TGFB2 corrected AUPR", maybe(x, x[test_ligand == "TGFB2", aupr_corrected]), 0.208170)
chk("S10A TGFB1 corrected AUPR", maybe(x, x[test_ligand == "TGFB1", aupr_corrected]), 0.203057)
chk("S10A AGT rank",             maybe(x, as.numeric(x[test_ligand == "AGT", rank])), 11, tol = 0)

x <- tsv("11A")
chk("S11A MMP14 permutation z",  maybe(x, x[test_ligand == "MMP14", perm_z]), 27.1049, tol = 1e-3)
chk("S11A MMP14 BH",             maybe(x, x[test_ligand == "MMP14", perm_p_BH]), 0.0321, tol = 1e-3)
chk("S11A n FDR-significant",    maybe(x, sum(x$fdr_significant)), 1, tol = 0)
chk("S11A TGFB2 BM rank",        maybe(x, as.numeric(x[test_ligand == "TGFB2", rank])), 5, tol = 0)
chk("S11A TGFB1 BM rank",        maybe(x, as.numeric(x[test_ligand == "TGFB1", rank])), 12, tol = 0)

x <- tsv("05A")
chk("S05A dropout obs/exp ratio", maybe(x, x$ratio_observed_to_expected[1]), 0.99112, tol = 1e-4)

x <- tsv("08C")
chk("S08C BM genes",      maybe(x, uniqueN(x$gene)), 13, tol = 0)
chk("S08C cell types",    maybe(x, uniqueN(x$ccc_group)), 22, tol = 0)
chk("S08C donor profiles", maybe(x, sum(x[gene == x$gene[1], n_profiles], na.rm = TRUE)), 2329, tol = 0)

x <- tsv("12A")
chk("S12A cell types",     maybe(x, uniqueN(x$ccc_group)), 22, tol = 0)
chk("S12A donor profiles", maybe(x, sum(x[gene == "AGT", n_profiles], na.rm = TRUE)), 4376, tol = 0)
chk("S12A max REN detection", maybe(x, max(x[gene == "REN", detect], na.rm = TRUE)), 0.0229, tol = 1e-3)

x <- tsv("12B")
chk("S12B AGTR1 detection in pericytes",
    maybe(x, x[ccc_group == "Pericytes", AGTR1_detection_fraction]), 0.3418, tol = 1e-3)
chk("S12B complete autonomous circuits", maybe(x, sum(x$autonomous_circuit)), 0, tol = 0)
chk("S12B max core roles", maybe(x, max(x$n_core_roles)), 1, tol = 0)

x <- tsv("09B")
chk("S09B BM-score contrasts at BH<0.05",
    maybe(x, sum(x[score == "basement_membrane_score_z", p_BH_within_score] < 0.05)), 12, tol = 0)

x <- tsv("02C1")
chk("S02C1 detection-vs-depth rho (pooled)",
    maybe(x, x[dataset_id == "ALL DATASETS POOLED" &
                   grepl("^Agtr1a detection", measure), spearman_rho]), 0.432, tol = 1e-3)

## -- threshold-dependent: assert against the >=20 SENSITIVITY rows, where the
##    previously published values still apply. The >=10 primary fit is new and has
##    no prior value to check, so it is validated structurally instead (below).
ge20 <- function(d) if (is.null(d) || !"analysis_role" %in% names(d)) NULL else
    d[grepl(">=20", analysis_role)]

x <- ge20(tsv("13E"))
chk("S13E (>=20) F statistic",
    maybe(x, x[block == "anova" & term == "disease_group", `F value`]), 5.96555, tol = 1e-3, TRUE)
chk("S13E (>=20) Healthy marginal mean",
    maybe(x, x[block == "emmeans" & disease_group == "Healthy", emmean]), -0.164118, tol = 1e-4, TRUE)
chk("S13E (>=20) Fibrotic_ILD marginal mean",
    maybe(x, x[block == "emmeans" & disease_group == "Fibrotic_ILD", emmean]), 0.564861, tol = 1e-4, TRUE)

x <- ge20(tsv("13D"))
chk("S13D (>=20) Healthy injury fraction",
    maybe(x, x[block == "marginal means" & disease_group == "Healthy", emmean]), 0.019157, tol = 1e-4, TRUE)

x <- ge20(tsv("14B"))
chk("S14B (>=20) niche index F",
    maybe(x, x[response == "niche_index" & block == "anova" &
                   term == "disease_group", `F value`]), 3.41048, tol = 1e-3, TRUE)

checks <- rbindlist(CHECKS)
cat("\n==== anchor checks ====\n"); print(checks)

## Structural check on the new primary fit: it must include strictly more donors
## than the sensitivity fit and keep the direction of the Healthy->Fibrotic effect.
struct_ok <- TRUE
e13 <- tsv("13E")
if (!is.null(e13) && "analysis_role" %in% names(e13) &&
    any(grepl("PRIMARY", e13$analysis_role)) && any(grepl(">=20", e13$analysis_role))) {
    nd <- e13[block == "emmeans", .(n = max(n_donors, na.rm = TRUE)),
              by = .(primary = grepl("PRIMARY", analysis_role))]
    if (nrow(nd) == 2 && nd[primary == TRUE, n] <= nd[primary == FALSE, n]) {
        struct_ok <- FALSE
        cat("FAIL: the >=10 primary fit does not have more donors than the >=20 fit\n")
    }
    mm <- e13[block == "emmeans" & disease_group %in% c("Healthy", "Fibrotic_ILD")]
    d <- mm[, .(diff = emmean[disease_group == "Fibrotic_ILD"] -
                    emmean[disease_group == "Healthy"]),
            by = .(primary = grepl("PRIMARY", analysis_role))]
    if (nrow(d) == 2 && prod(sign(d$diff)) < 0) {
        struct_ok <- FALSE
        cat("FAIL: Healthy->Fibrotic effect changes sign between thresholds\n")
    }
}

failed <- checks[passed == FALSE]
if (nrow(failed) || !struct_ok) {
    print(failed)
    stop("anchor checks failed -- the tables disagree with the manuscript. ",
         "Fix the source or the text before shipping.")
}
cat("all", nrow(checks), "anchor checks passed\n")

## =========================================================================
## Workbook
## =========================================================================
pending <- man[status != "complete"]
if (nrow(pending)) {
    cat("\n==== PENDING PARTS (upstream re-runs outstanding) ====\n")
    print(pending[, .(part, title, status)])
}

sheets <- list()
readme <- man[, .(Table = paste0("S", part), Title = title,
                  Supports = supports, Rows = n_rows, Columns = n_cols,
                  Status = status, Sources = sources, Notes = notes)]
sheets[["README"]] <- as.data.frame(readme)
for (i in seq_len(nrow(man))) {
    d <- tsv(man$part[i])
    if (is.null(d)) { warning("missing TSV for part ", man$part[i]); next }
    ## Excel sheet names are capped at 31 characters and must be unique.
    nm <- paste0("S", man$part[i])
    sheets[[nm]] <- as.data.frame(d)
}
xlsx <- file.path(OUT, "supplementary_tables.xlsx")
write_xlsx(sheets, xlsx)
cat("\nwrote", xlsx, "with", length(sheets), "sheets\n")

write.table(man[, .(part, title, supports, n_rows, n_cols, status, file, sources)],
            file.path(OUT, "table_manifest.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE, na = "")
write.table(checks, file.path(OUT, "anchor_checks.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
print(man[, .(part, n_rows, n_cols, status)])

cat("\nReproducibility information:\n")
Sys.time(); options(width = 120); sessioninfo::session_info()
