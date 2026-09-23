## Supplementary Table S15 (parts A-F): the distributed-RAS circuit behind Figure 5
## and its robustness supplement, Figure S18. Source: ras_circuit/_m.
##
## S14 is vacant (moved with the disease layer on 2026-09-10); S15 is the next free
## table number, taken here per the append-only convention used for figures.
##
## A part whose source is missing is written with status `pending_upstream`, which
## makes 08.assemble_tables.R refuse the final workbook -- a missing module run is
## loud, never a silently absent sheet. Model descriptions are DERIVED from the
## `model` column of each source (the P2-35 rule), never asserted in a note.

suppressPackageStartupMessages({ library(data.table); library(dplyr) })
source("../_h/_tab_common.R")

RC <- function(...) P("ras_circuit", "_m", ...)
SD <- function(f) RC("stats_data", f)
models_of <- function(x) {
    if (is.null(x) || !"model" %in% names(x)) return("Model: see source script header.")
    m <- unique(na.omit(as.character(x$model))); m <- m[nzchar(m)]
    paste0("Models (from the `model` column): ", paste(m, collapse = "; "), ".")
}
part <- function(files, id, title, supports, notes_fun, combine = NULL) {
    xs <- lapply(files, read_src)
    names(xs) <- basename(files)
    if (any(vapply(xs, is.null, logical(1)))) {
        write_part(data.frame(missing_source = files[vapply(xs, is.null, logical(1))]),
                   id, title, supports = supports, sources = sub(paste0(ROOT, "/"), "", files),
                   status = "pending_upstream",
                   notes = "Source missing -- run the ras_circuit module (steps 1-4) first.")
        return(invisible(NULL))
    }
    df <- if (is.null(combine)) rbindlist(Map(function(x, n) x[, source_table := n], xs, names(xs)),
                                          fill = TRUE) else combine(xs)
    write_part(df, id, title, supports = supports,
               sources = sub(paste0(ROOT, "/"), "", files), notes = notes_fun(df, xs))
}

## ---- S15A niche affinity ------------------------------------------------------------
part(c(SD("niche_affinity_global.tsv"), SD("niche_affinity_agtr1_models.tsv"),
       SD("niche_affinity_contrasts.tsv"), SD("niche_affinity_per_axis_models.tsv"),
       SD("niche_affinity_cell_models.tsv"), SD("niche_affinity_legacy_comparison.tsv")),
     "15A", "Niche-affinity decomposition of the pericyte airspace score against AGTR1",
     "Figure 5B; Figure S18A",
     function(df, xs) {
         g <- xs[["niche_affinity_global.tsv"]][arm == "primary"]
         paste0("Four cosine similarities in X_pca_harmony (AT1, AT2, aerocyte EC, general ",
                "capillary EC) whose mean reproduces the published airspace_score. ",
                "Primary: ", g$n_donors, " donors / ", g$n_studies, " studies, no age term ",
                "(age is a study filter here; see `_ageadj` for the restricted arm and its ",
                "`studies_dropped`). `p_emp` refers each primary statistic to ",
                g$n_null, " detection-matched genes substituted for AGTR1. ", models_of(df))
     })

## ---- S15B RAS network -----------------------------------------------------------------------
part(c(SD("ras_network_edges.tsv"), SD("ras_network_chain_models.tsv"),
       SD("ras_dag_nodes.tsv"), SD("ras_dag_edges.tsv"), SD("ras_dag_scalars.tsv"),
       SD("ras_dag_annotations.tsv"), SD("ras_dag_notes.tsv")),
     "15B", "Donor-level RAS covariance restricted to the proposed DAG",
     "Figure 5A, 5C; Figure S18D",
     function(df, xs) {
         e <- xs[["ras_network_edges.tsv"]]
         paste0("Partial Spearman correlations over donors, each node residualized on its ",
                "own unit depth with (1 | study). `family`: allowed = DAG edges (BH within); ",
                "forbidden = implied conditional independences, expected null; ",
                "marginal_descriptive = individual RAS nodes, not tested as a family. ",
                "Bootstrap: donors resampled with residuals held fixed. Leave-one-dataset-out ",
                "over ", e$n_datasets_lodo[1], " datasets. Consistency with a pre-specified DAG, ",
                "not causal mediation; Ang I / Ang II are latent. `ras_dag_edges` carries the ",
                "PRIOR sign each drawn edge is given in Figure 5A (arrow head vs bar head); ",
                "`ras_dag_annotations` carries the estimate that annotates it, with the panel ",
                "it came from, so a prior and its estimate can disagree on the page. ",
                models_of(df))
     })

## ---- S15C signature provenance -----------------------------------------------------------
part(c(P("ras_circuit", "_h", "signatures", "angii_response_signature.tsv"),
       SD("at1r_signature_audit.tsv")),
     "15C", "AngII response signature: provenance and disjointness audit",
     "Figure 5D; Figure S18F",
     function(df, xs) {
         a <- xs[["at1r_signature_audit.tsv"]]
         s <- xs[["angii_response_signature.tsv"]]
         paste0("Derived from ", s$source_id[1], "; contrast ", s$contrast[1], "; cell type ",
                s$cell_type[1], ". Rule pre-specified in ras_circuit/_h/signatures/README.md. ",
                sum(a$kept), " of ", nrow(a), " genes kept after pruning genes shared with ",
                "state, matrix, tracer, TGF-beta or RAS panels (`pruned_reason`).")
     })

## ---- S15D AT1R response -------------------------------------------------------------------
part(c(SD("at1r_vs_agtr1_count.tsv"), SD("at1r_vs_agtr1_scvi.tsv"),
       SD("at1r_continuum_summary.tsv"), SD("at1r_by_subcluster_emmeans.tsv"),
       SD("at1r_by_subcluster_posthoc.tsv"), SD("at1r_vs_programs.tsv")),
     "15D", "Pericyte AT1R-response programme vs AGTR1, the continuum and subclusters",
     "Figure 5D; Figure S18B",
     function(df, xs) {
         a <- xs[["at1r_vs_agtr1_count.tsv"]][score == "at1r_response_score" &
                                                 model %in% c("NB GLMM", "Poisson+OLRE")]
         paste0("Arbiter: AGTR1 integer counts as the response of an NB GLMM with a ",
                "library-size offset (", a$model[1], "). The signature's `p_emp` refers its ",
                "estimate to ", a$n_null[1], " detection-matched null panels (null centre ",
                sprintf("%.3f", a$null_mean[1]), "); PROGENy and CollecTRI rows are tested ",
                "against zero and are secondary. Pseudotime is DPT rooted at the ",
                "vascular-stabilizing pole. ", models_of(df))
     })

## ---- S15E matrix consequence -----------------------------------------------------------------
part(c(SD("matrix_vs_at1r.tsv"), SD("matrix_vs_at1r_donor_rho_summary.tsv"),
       SD("matrix_context_rows.tsv")),
     "15E", "Basement-membrane versus fibrillar matrix along the AT1R-response axis",
     "Figure 5E; Figure S18C",
     function(df, xs) {
         paste0("The estimand is the BM - fibrillar CONTRAST; rows with claimable = FALSE ",
                "(either score alone) are reported for completeness only. `p_emp` against the ",
                "detection-matched null. `matrix_context_rows` are copied from ",
                "basement_membrane/_m/stats_data for the legend caveats, not recomputed. ",
                models_of(df))
     })

## ---- S15F outgoing communication -----------------------------------------------------------
part(c(RC("liana_direction_summary.tsv"),
       RC("nichenet_outgoing", "outgoing_nichenet_summary.tsv"),
       RC("nichenet_outgoing", "program_genes_dropped.tsv"),
       SD("outgoing_donor_validation.tsv")),
     "15F", "Outgoing pericyte communication: LIANA edges, NicheNet (pericyte as sender), donor validation",
     "Figure 5F; Figure S18E",
     function(df, xs) {
         paste0("LIANA edges re-filtered from the existing all-pairs run (no re-run); ",
                "COPA, MMP14, SIRPB2 flagged `not_a_ligand` and AGT `substrate_not_ligand`. ",
                "NicheNet ligands ordered by permutation z. Programme genes that pericytes ",
                "could themselves be scored for as ligands were dropped (listed). Donor ",
                "validation `p_emp` against detection-matched non-ligand composites. ",
                models_of(df))
     })
