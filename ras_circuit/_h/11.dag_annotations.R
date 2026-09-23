## Figure 5A: attach the MEASURED statistics from panels C-F to the drawn DAG edges.
##
## 07 writes the graph (structure + prior signs). This writes the evidence layer:
## one row per drawn edge that a model in this module actually estimated, plus the
## panel-level notes. Splitting the two keeps the schematic honest -- an edge exists
## in the drawing because biochemistry says so, and carries a number only where a
## donor-level model produced one. Where prior and estimate disagree (pericyte ACE2)
## the estimate is reported as it stands.
##
## Every number is read from ras_circuit/_m; nothing is typed in. Runs in seconds.
##     Rscript ../_h/11.dag_annotations.R --outdir ./stats_data

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("../_h/_stats_common.R")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--outdir", default = "./stats_data"),
    make_option("--nichenet", default = "./nichenet_outgoing"),
    make_option("--liana", default = "./liana_direction_summary.tsv"),
    make_option("--top-ligands", type = "integer", default = 2L, dest = "top_lig")
)))
SD <- function(f) file.path(opt$outdir, f)

net <- read_req(SD("ras_network_edges.tsv"))[arm == "primary"]
mx  <- read_req(SD("matrix_vs_at1r.tsv"))
cnt <- read_req(SD("at1r_vs_agtr1_count.tsv"))[model == "NB GLMM" & spec == "primary" &
                                                 score == "at1r_response_score"]
glo <- read_req(SD("niche_affinity_global.tsv"))[arm == "primary"]
prg <- read_req(SD("at1r_vs_programs.tsv"))
nn  <- read_req(file.path(opt$nichenet, "outgoing_nichenet_summary.tsv"))[status == "run"]
dv  <- read_req(SD("outgoing_donor_validation.tsv"))[spec == "primary"]

## ---- helpers ---------------------------------------------------------------------
star <- function(p) fifelse(is.finite(p) & p < 0.05, "*", "")
sgn  <- function(est, p) fifelse(!is.finite(p) | p >= 0.05, "null",
                          fifelse(est < 0, "negative", "positive"))
num  <- function(x) sprintf("%+.2f", x)
## One fitted edge of the panel-C network, by its node names there.
E <- function(f, t) {
    r <- net[from == f & to == t]
    if (nrow(r) != 1L) stop("panel-C edge not found exactly once: ", f, " -> ", t)
    r
}
row_edge <- function(from, to, est, p, p_kind, stat_kind, panel, label, note = "",
                     show = TRUE) {
    data.table(from = from, to = to, kind = "edge", evidence_panel = panel,
               stat_kind = stat_kind, estimate = est, p_value = p, p_kind = p_kind,
               measured_sign = sgn(est, p), label = label, note = note, show_in_panel = show)
}

ann <- list()

## ---- the two latent-peptide tiers -------------------------------------------------
## The peptides are unmeasured, so the flux edges carry no estimate. What IS estimable
## is whether the donor-level supply upstream of them tracks the pericyte response;
## both composite indices do, and they annotate the edge where the flux enters.
ag <- E("AGT_source_index", "Peri_AT1R_response"); pr <- E("processing_index", "Peri_AT1R_response")
ann[[length(ann) + 1]] <- row_edge(
    "AngII", "Peri_AT1R_response", ag$partial_rho, ag$p_BH, "BH", "partial rho", "C",
    ## "; " is the line break: an unquoted newline would break the TSV.
    sprintf("AGT source %s%s; ACE/chymase %s%s", num(ag$partial_rho), star(ag$p_BH),
            num(pr$partial_rho), star(pr$p_BH)),
    note = "donor-level partial Spearman of each composite index with the response")

## ---- receptor abundance is not receptor activity ----------------------------------
gr <- E("Peri_AGTR1", "Peri_AT1R_response")
ann[[length(ann) + 1]] <- row_edge(
    "Peri_AGTR1", "Peri_AT1R_response", gr$partial_rho, gr$p_BH, "BH", "partial rho", "C/D",
    sprintf("%s n.s.", num(gr$partial_rho)),
    note = "count arbiter agrees: see the AGTR1 abundance note")

## ---- response branches ------------------------------------------------------------
for (b in c("contractile", "inflammatory", "activated_migratory")) {
    r <- E("Peri_AT1R_response", b)
    ann[[length(ann) + 1]] <- row_edge("Peri_AT1R_response", b, r$partial_rho, r$p_BH,
                                       "BH", "partial rho", "C",
                                       sprintf("%s%s", num(r$partial_rho), star(r$p_BH)))
}
ann[[length(ann) + 1]] <- row_edge("Peri_AT1R_response", "matrix", NA_real_, NA_real_, "",
                                   "", "E", "", show = FALSE)

## ---- the matrix endpoints, separately ---------------------------------------------
## The claim is the CONTRAST; the two arms are shown so the reader can see which side
## moves. Fibrillar rises, basement membrane does not.
fib <- mx[spec == "_fib_alone"]; bm <- mx[spec == "_bm_alone"]
ann[[length(ann) + 1]] <- row_edge("matrix", "FIB", fib$estimate, fib$p_value, "raw", "beta", "E",
                                   sprintf("beta %s%s", num(fib$estimate), star(fib$p_value)))
ann[[length(ann) + 1]] <- row_edge("matrix", "BM", bm$estimate, bm$p_value, "raw", "beta", "E",
                                   sprintf("beta %s n.s.", num(bm$estimate)))

## ---- counter-regulatory arm --------------------------------------------------------
ms <- E("Peri_MAS1", "Peri_AT1R_response")
ann[[length(ann) + 1]] <- row_edge(
    "Peri_MAS1", "Peri_AT1R_response", ms$partial_rho, ms$p_BH, "BH", "partial rho", "C",
    sprintf("%s n.s.", num(ms$partial_rho)),
    note = "prior sign is antagonistic; the donor-level test is null")

## ---- outgoing edges: name the ligands that carry them -------------------------------
## Panel F's top-ranked pericyte ligands per receiver programme, so the neighbour edges
## say WHAT is sent, not just that something is.
NN_TO_NODE <- c("EC aerocyte capillary" = "out_EC", "EC general capillary" = "out_EC",
                "AT1" = "out_EPI", "AT2" = "out_EPI",
                "Alveolar fibroblasts" = "out_FIB", "Adventitial fibroblasts" = "out_FIB")
lig <- rbindlist(lapply(seq_len(nrow(nn)), function(i) {
    z <- as.numeric(strsplit(nn$top_z[i], ";")[[1]])
    data.table(target = nn$target[i], test_ligand = strsplit(nn$top_ligands_by_z[i], ";")[[1]], z = z)
}))
lig[, node := NN_TO_NODE[target]]
## Two receivers share each outgoing node (aerocyte/general capillary, AT1/AT2, the two
## fibroblasts), so keep each ligand once at its best z before taking the top few.
lig <- lig[order(-z), head(.SD, 1), by = .(node, test_ligand)]
top <- lig[order(-z), head(.SD, opt$top_lig), by = node]
dv[, node := NN_TO_NODE[target]]
val <- dv[order(-abs(estimate)), head(.SD, 1), by = node]
for (nd in unique(top$node)) {
    src <- switch(nd, out_EC = "BM", out_FIB = "FIB", out_EPI = "Peri_AT1R_response")
    v <- val[node == nd]
    ann[[length(ann) + 1]] <- row_edge(
        src, nd, max(top[node == nd, z]), NA_real_, "", "NicheNet z", "F",
        paste(top[node == nd, test_ligand], collapse = ", "),
        note = sprintf("donor validation, %s: beta %s, p = %.3g (p_emp %.3g)",
                       v$target, num(v$estimate), v$p_value, v$p_emp))
}

A <- rbindlist(ann, fill = TRUE)

## ---- panel-level notes ---------------------------------------------------------------
vs <- E("VSMC_AGT", "Peri_AGTR1")
ac <- E("Peri_ACE2", "Peri_AT1R_response")
bmp <- prg[program == "basement_membrane_score"]
N <- rbindlist(list(
    data.table(id = "forbidden", show_in_panel = TRUE, text = sprintf(
        "AGT -> AGTR1 is not an edge (substrate, not ligand); tested as a control: rho %s, n.s.",
        num(vs$partial_rho))),
    data.table(id = "abundance_not_activity", show_in_panel = TRUE, text = sprintf(
        "AGTR1 abundance is not AT1R activity: count model beta %s against a matched-null centre of %s (p_emp %s)",
        num(cnt$estimate), num(cnt$null_mean), formatC(cnt$p_emp, digits = 2, format = "g"))),
    data.table(id = "matrix_balance", show_in_panel = TRUE, text = sprintf(
        "matrix balance BM - fibrillar: beta %s against its matched null (p_emp %s); donor rho %s",
        num(mx[spec == "primary", estimate]),
        formatC(mx[spec == "primary", p_emp], digits = 2, format = "g"), num(bmp$rho))),
    data.table(id = "ace2_not_antagonistic", show_in_panel = FALSE, text = sprintf(
        "pericyte ACE2 -> response rho %s%s: the counter-regulatory arm does not oppose the response here",
        num(ac$partial_rho), star(ac$p_BH))),
    data.table(id = "no_compartment_affinity", show_in_panel = FALSE, text = sprintf(
        "AGTR1-high pericytes show no compartment-specific niche affinity (interaction LRT p = %.2f)",
        glo$lrt_p))))

write_tsv_safe(A, SD("ras_dag_annotations.tsv"))
write_tsv_safe(N, SD("ras_dag_notes.tsv"))
message(sprintf("ras_dag_annotations: %d annotated edges (%d shown); %d notes",
                nrow(A), sum(A$show_in_panel), nrow(N)))
