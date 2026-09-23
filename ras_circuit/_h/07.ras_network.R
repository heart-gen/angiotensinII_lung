## Figures 5A and 5C -- the multicellular RAS as a biology-constrained DAG.
##
## (A) Data-annotated DAG tables for the schematic: every drawn node that is a
##     (cell type, gene) pair carries its detection and model estimate from
##     agt_axis (ras_celltype_profile.tsv). Ang I, Ang II and Ang 1-7 are LATENT
##     (peptides; invisible to scRNA-seq). AGT -> AGTR1 is a FORBIDDEN edge and is
##     never drawn: angiotensinogen is renin's substrate, not the AT1R ligand.
##
## (C) Donor-level covariance restricted to that DAG. One row per donor. Every node
##     is residualized on its own unit's sequencing depth with (1 | study) (arms:
##     + disease_group as a nuisance; Healthy donors only). Because the peptides are
##     latent, their parents are summarised as two indices (brief: "AGT-source
##     index", "processing index"):
##       AGT_source_index = mean z of AGT in VSMC, alveolar and adventitial fibroblasts
##       processing_index = mean z of ACE in aerocyte/general capillary EC and
##                          alveolar/interstitial macrophages, CMA1 and CTSG in mast cells
##     Edge statistic = Spearman partial correlation of child and parent given the
##     child's OTHER DAG parents (complete cases for that parent set). BH over the
##     allowed-edge family. Implied conditional independencies are tested as
##     "forbidden" controls and should be null.
##     Robustness: 1,000 donor bootstraps (residuals held fixed), leave-one-dataset-
##     out over the datasets actually present, a sequential linear chain, and
##     dagitty::localTests where installed.
##
## This is CONSISTENCY WITH A PROPOSED DAG, not causal mediation: the DAG is
## specified by biochemistry, not learned (no PC / NOTEARS), and Ang II is unmeasured.

suppressPackageStartupMessages({
    library(optparse); library(data.table); library(lme4)
})
source("../_h/_stats_common.R")

opt <- parse_args(OptionParser(option_list = list(
    make_option("--pseudobulk", default = "./ras_circuit_pseudobulk.tsv.gz"),
    make_option("--programs", default = "./receiver_programs.tsv"),
    make_option("--meta", default = "./at1r_response_metadata.tsv.gz"),
    make_option("--agt", default = "../../agt_axis/_m/stats_data"),
    make_option("--outdir", default = "./stats_data"),
    make_option("--min-cells", type = "integer", default = 5L, dest = "min_cells"),
    make_option("--min-pericytes", type = "integer", default = 10L, dest = "min_peri"),
    make_option("--nboot", type = "integer", default = 1000L),
    make_option("--seed", type = "integer", default = 13L)
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(opt$seed)

## =============================== (A) DAG tables ==================================
prof <- read_req(file.path(opt$agt, "ras_celltype_profile.tsv"))
comp <- read_req(file.path(opt$agt, "ras_circuit_completeness.tsv"))
look <- function(ct, g) {
    r <- prof[ccc_group == ct & gene == g]
    if (!nrow(r)) return(list(detect = NA_real_, emmean = NA_real_, rank = NA_integer_))
    rk <- prof[gene == g & !is.na(emmean)][, rk := frank(-emmean, ties.method = "min")][
        ccc_group == ct, rk]
    list(detect = r$detect[1], emmean = r$emmean[1],
         rank = if (length(rk)) as.integer(rk) else NA_integer_)
}
N <- list(
    ## tier 1: AGT sources
    list("VSMC_AGT", "AGT, VSMC", "Vascular smooth muscle", "AGT", "substrate", 1, 3),
    list("AlvFib_AGT", "AGT, alveolar fibroblast", "Alveolar fibroblasts", "AGT", "substrate", 1, 2),
    list("AdvFib_AGT", "AGT, adventitial fibroblast", "Adventitial fibroblasts", "AGT", "substrate", 1, 1),
    list("AngI", "Ang I (latent)", NA, NA, "latent", 2, 2),
    ## tier 3: processing
    list("ECaero_ACE", "ACE, aerocyte capillary EC", "EC aerocyte capillary", "ACE", "processing", 3, 4),
    list("ECgcap_ACE", "ACE, general capillary EC", "EC general capillary", "ACE", "processing", 3, 3),
    list("AlvMac_ACE", "ACE, alveolar macrophage", "Alveolar macrophages", "ACE", "processing", 3, 2),
    list("Mast_CMA1", "CMA1, mast cell", "Mast cells", "CMA1", "processing", 3, 1),
    list("Mast_CTSG", "CTSG, mast cell", "Mast cells", "CTSG", "processing", 3, 0),
    list("AngII", "Ang II (latent)", NA, NA, "latent", 4, 2),
    ## tier 5: receptor + counter-regulatory
    list("Peri_AGTR1", "AGTR1, pericyte", "Pericytes", "AGTR1", "receptor", 5, 2),
    list("Peri_ACE2", "ACE2, pericyte", "Pericytes", "ACE2", "counter_regulatory", 5, 0.5),
    list("Ang17", "Ang 1-7 (latent)", NA, NA, "latent", 6, 0.5),
    list("Peri_MAS1", "MAS1, pericyte", "Pericytes", "MAS1", "counter_regulatory", 7, 0.5),
    ## tier 6-9: response and consequences
    list("Peri_AT1R_response", "AT1R-response programme", "Pericytes", NA, "response", 6, 2),
    list("contractile", "contractile", "Pericytes", NA, "response_branch", 7, 3.2),
    list("inflammatory", "inflammatory", "Pericytes", NA, "response_branch", 7, 2.4),
    list("matrix", "matrix remodelling", "Pericytes", NA, "response_branch", 7, 1.6),
    list("BM", "basement membrane", "Pericytes", NA, "matrix", 8, 2.2),
    list("FIB", "fibrillar ECM", "Pericytes", NA, "matrix", 8, 1.2),
    list("out_EC", "capillary EC", NA, NA, "neighbor", 9, 3),
    list("out_FIB", "fibroblasts", NA, NA, "neighbor", 9, 2),
    list("out_EPI", "AT1 / AT2", NA, NA, "neighbor", 9, 1))
nodes <- rbindlist(lapply(N, function(z) {
    ct <- z[[3]]; g <- z[[4]]
    v <- if (!is.na(ct) && !is.na(g)) look(ct, g) else list(detect = NA, emmean = NA, rank = NA)
    data.table(node = z[[1]], label = z[[2]], cell_type = ct, gene = g, node_class = z[[5]],
               x = z[[6]], y = z[[7]], detect = v$detect, emmean = v$emmean,
               rank_in_gene = v$rank)
}))
E <- rbind(
    data.table(from = c("VSMC_AGT", "AlvFib_AGT", "AdvFib_AGT"), to = "AngI", edge_class = "substrate"),
    data.table(from = "AngI", to = "AngII", edge_class = "processing"),
    data.table(from = c("ECaero_ACE", "ECgcap_ACE", "AlvMac_ACE", "Mast_CMA1", "Mast_CTSG"),
               to = "AngII", edge_class = "processing"),
    data.table(from = c("AngII", "Peri_AGTR1"), to = "Peri_AT1R_response", edge_class = "receptor"),
    data.table(from = "Peri_AT1R_response", to = c("contractile", "inflammatory", "matrix"),
               edge_class = "response"),
    data.table(from = "matrix", to = c("BM", "FIB"), edge_class = "matrix"),
    data.table(from = c("BM", "FIB"), to = c("out_EC", "out_FIB"), edge_class = "neighbor"),
    data.table(from = "Peri_AT1R_response", to = "out_EPI", edge_class = "neighbor"),
    data.table(from = c("AngII", "Peri_ACE2", "Ang17"), to = c("Ang17", "Ang17", "Peri_MAS1"),
               edge_class = "counter_regulatory"),
    data.table(from = "VSMC_AGT", to = "Peri_AGTR1", edge_class = "forbidden"))
E[, drawn := edge_class != "forbidden"]
stopifnot(all(c(E$from, E$to) %in% nodes$node))
write_tsv_safe(nodes, file.path(opt$outdir, "ras_dag_nodes.tsv"))
write_tsv_safe(E, file.path(opt$outdir, "ras_dag_edges.tsv"))
scal <- data.table(
    quantity = c("max_REN_detect", "max_REN_cell_type", "n_cell_types",
                 "n_autonomous_circuit", "autonomy_detect_threshold", "max_steps_in_one_cell_type"),
    value = c(sprintf("%.4f", max(prof[gene == "REN", detect], na.rm = TRUE)),
              prof[gene == "REN"][which.max(detect), ccc_group],
              as.character(uniqueN(comp$ccc_group)),
              as.character(sum(as.logical(comp$autonomous_circuit), na.rm = TRUE)),
              "0.05 (agt_axis/_h/01.ras_landscape_stats.R --detect-thr default)",
              as.character(max(comp$n_steps_present, na.rm = TRUE))))
write_tsv_safe(scal, file.path(opt$outdir, "ras_dag_scalars.tsv"))
sw <- file.path(opt$agt, "ras_autonomy_threshold_sweep.tsv")
if (file.exists(sw)) file.copy(sw, file.path(opt$outdir, "ras_autonomy_threshold_sweep.tsv"),
                               overwrite = TRUE)

## ============================ (C) donor node table ===============================
pb <- read_req(opt$pseudobulk)[n_cells >= opt$min_cells]
prog <- read_req(opt$programs)
unit_node <- function(ct, g, nm) {
    col <- paste0(g, "__expr")
    x <- pb[ccc_group == ct, c("donor_id", col, "mean_log10_total_counts"), with = FALSE]
    setnames(x, c("donor_id", nm, paste0("depth__", nm))); x
}
spec <- list(c("Vascular smooth muscle", "AGT", "VSMC_AGT"),
             c("Alveolar fibroblasts", "AGT", "AlvFib_AGT"),
             c("Adventitial fibroblasts", "AGT", "AdvFib_AGT"),
             c("EC aerocyte capillary", "ACE", "ECaero_ACE"),
             c("EC general capillary", "ACE", "ECgcap_ACE"),
             c("Alveolar macrophages", "ACE", "AlvMac_ACE"),
             c("Interstitial macrophages", "ACE", "IntMac_ACE"),
             c("Mast cells", "CMA1", "Mast_CMA1"),
             c("Mast cells", "CTSG", "Mast_CTSG"),
             c("Pericytes", "AGTR1", "Peri_AGTR1"),
             c("Pericytes", "ACE2", "Peri_ACE2"),
             c("Pericytes", "MAS1", "Peri_MAS1"))
tabs <- lapply(spec, function(s) unit_node(s[1], s[2], s[3]))
## EC readout: mean z of the barrier genes in each capillary EC type, then averaged
ecg <- prog[program == "EC_readout", gene]
ec <- rbindlist(lapply(c("EC aerocyte capillary", "EC general capillary"), function(ct) {
    x <- pb[ccc_group == ct]
    cols <- intersect(paste0(ecg, "__expr"), names(x))
    z <- sapply(cols, function(c) z_within_dataset(x[[c]], x$dataset))
    data.table(donor_id = x$donor_id, v = rowMeans(as.matrix(z), na.rm = TRUE),
               dep = x$mean_log10_total_counts)
}))
ec <- ec[, .(EC_readout = mean(v), depth__EC_readout = mean(dep)), by = donor_id]
tabs[[length(tabs) + 1]] <- ec

m <- read_req(opt$meta); setnames(m, 1, "index")
m[, bm_minus_fib := basement_membrane_score - fibrillar_collagen_score]
peri <- m[, .(n_peri = .N,
              Peri_AT1R_response = mean(at1r_response_score),
              contractile = mean(synthetic_contractile_score),
              inflammatory = mean(inflammatory_score),
              activated_migratory = mean(activated_migratory_score),
              vascular_stabilizing = mean(vascular_stabilizing_score),
              BM = mean(basement_membrane_score), FIB = mean(fibrillar_collagen_score),
              BM_minus_FIB = mean(bm_minus_fib),
              depth_peri = mean(log10_total_counts)), by = donor_id][n_peri >= opt$min_peri]
for (v in c("Peri_AT1R_response", "contractile", "inflammatory", "activated_migratory",
            "vascular_stabilizing", "BM", "FIB", "BM_minus_FIB"))
    peri[[paste0("depth__", v)]] <- peri$depth_peri
tabs[[length(tabs) + 1]] <- peri[, !c("n_peri", "depth_peri")]

meta_d <- unique(pb[ccc_group == "Pericytes", .(donor_id, study, dataset, disease_group)])
nd <- Reduce(function(a, b) merge(a, b, by = "donor_id", all = TRUE), tabs)
nd <- merge(meta_d, nd, by = "donor_id")          # pericyte-bearing donors only
fwrite(nd, "./ras_network_donor_nodes.tsv.gz", sep = "\t")
NODES <- setdiff(grep("^depth__", names(nd), value = TRUE, invert = TRUE),
                 c("donor_id", "study", "dataset", "disease_group"))
message(sprintf("donor node table: %d pericyte-bearing donors; donors per node:", nrow(nd)))
print(sapply(NODES, function(v) sum(is.finite(nd[[v]]))))

## ---- residualize -----------------------------------------------------------------------
residualize <- function(d, adj_disease = FALSE) {
    R <- d[, .(donor_id, study, dataset, disease_group)]
    for (v in NODES) {
        dv <- d[is.finite(get(v)) & is.finite(get(paste0("depth__", v)))]
        if (nrow(dv) < 15) {                 # e.g. mast cells in the Healthy-only arm
            message(sprintf("  node %s: %d donors, not residualized (left NA)", v, nrow(dv)))
            R[, (v) := NA_real_]; next
        }
        dv[, y := get(v)]; dv[, dep := get(paste0("depth__", v))]
        f <- if (adj_disease && uniqueN(dv$disease_group) > 1)
            y ~ dep + disease_group + (1 | study) else y ~ dep + (1 | study)
        fit <- try(suppressMessages(lmer(f, data = dv)), silent = TRUE)
        ## If the mixed model fails, keep the study guard as a fixed effect rather
        ## than dropping it.
        r <- if (inherits(fit, "try-error"))
            resid(lm(update(nobars(f), . ~ . + study), data = dv)) else resid(fit)
        R[dv, (v) := r / sd(r), on = "donor_id"]
    }
    R[, AGT_source_index := {
        M <- as.matrix(.SD); k <- rowSums(is.finite(M))
        fifelse(k >= 2, rowMeans(M, na.rm = TRUE), NA_real_)
    }, .SDcols = c("VSMC_AGT", "AlvFib_AGT", "AdvFib_AGT")]
    R[, processing_index := {
        M <- as.matrix(.SD); k <- rowSums(is.finite(M))
        fifelse(k >= 3, rowMeans(M, na.rm = TRUE), NA_real_)
    }, .SDcols = c("ECaero_ACE", "ECgcap_ACE", "AlvMac_ACE", "IntMac_ACE", "Mast_CMA1", "Mast_CTSG")]
    R
}

## child <- parents; the allowed-edge family
DAG <- list(
    Peri_AT1R_response = c("AGT_source_index", "processing_index", "Peri_AGTR1", "Peri_ACE2"),
    contractile = "Peri_AT1R_response", inflammatory = "Peri_AT1R_response",
    activated_migratory = "Peri_AT1R_response",
    BM_minus_FIB = "Peri_AT1R_response",
    EC_readout = c("BM_minus_FIB", "Peri_AT1R_response"))
ALLOWED <- rbindlist(lapply(names(DAG), function(ch)
    data.table(from = DAG[[ch]], to = ch,
               given = vapply(DAG[[ch]], function(p) paste(setdiff(DAG[[ch]], p), collapse = "+"), ""))))
ALLOWED[, family := "allowed"]
## implied conditional independencies -> should be null
FORBID <- data.table(
    from = c("AGT_source_index", "processing_index", "Mast_CMA1", "VSMC_AGT", "AGT_source_index"),
    to = c("BM_minus_FIB", "EC_readout", "BM_minus_FIB", "Peri_AGTR1", "Peri_AGTR1"),
    given = c("Peri_AT1R_response", "BM_minus_FIB+Peri_AT1R_response", "Peri_AT1R_response", "", ""),
    family = "forbidden")
## marginal, descriptive: each individual RAS node against the response
MARG <- data.table(from = c("VSMC_AGT", "AlvFib_AGT", "AdvFib_AGT", "ECaero_ACE", "ECgcap_ACE",
                            "AlvMac_ACE", "IntMac_ACE", "Mast_CMA1", "Mast_CTSG", "Peri_MAS1"),
                   to = "Peri_AT1R_response", given = "", family = "marginal_descriptive")
EDGES <- rbind(ALLOWED, FORBID, MARG)

pcor <- function(R, from, to, given) {
    g <- if (nzchar(given)) strsplit(given, "\\+")[[1]] else character(0)
    cols <- c(from, to, g)
    X <- R[, ..cols]; X <- X[complete.cases(X)]
    n <- nrow(X)
    if (n < 15) return(list(rho = NA_real_, p = NA_real_, n = n))
    rk <- as.data.table(lapply(X, rank))
    if (length(g)) {
        fx <- resid(lm(reformulate(g, from), data = rk))
        fy <- resid(lm(reformulate(g, to), data = rk))
    } else { fx <- rk[[from]]; fy <- rk[[to]] }
    ct <- suppressWarnings(cor.test(fx, fy))
    ## df correction for the conditioning set
    r <- unname(ct$estimate); df <- n - 2 - length(g)
    tt <- r * sqrt(df / max(1e-12, 1 - r^2))
    list(rho = r, p = 2 * pt(-abs(tt), df), n = n)
}
edge_table <- function(R) {
    out <- copy(EDGES)
    st <- lapply(seq_len(nrow(out)), function(i) pcor(R, out$from[i], out$to[i], out$given[i]))
    out[, `:=`(partial_rho = sapply(st, `[[`, "rho"), p_value = sapply(st, `[[`, "p"),
               n_donors = sapply(st, `[[`, "n"))]
    out[, p_BH := p.adjust(p_value, "BH"), by = family]
    out
}

arms <- list(primary = residualize(nd),
             `_disease_adj` = residualize(nd, adj_disease = TRUE),
             `_healthy` = residualize(nd[disease_group == "Healthy"]))
edges <- rbindlist(lapply(names(arms), function(a) edge_table(arms[[a]])[, arm := a]))

## ---- bootstrap + LODO on the primary residuals ------------------------------------------
R0 <- arms$primary
boot <- rbindlist(lapply(seq_len(opt$nboot), function(b) {
    Rb <- R0[sample.int(nrow(R0), replace = TRUE)]
    e <- edge_table(Rb)[, .(from, to, family, partial_rho)]
    e[, boot := b]
}))
fwrite(boot, file.path(opt$outdir, "ras_network_bootstrap.tsv"), sep = "\t")
bs <- boot[, .(boot_lo = quantile(partial_rho, 0.025, na.rm = TRUE),
               boot_hi = quantile(partial_rho, 0.975, na.rm = TRUE)), by = .(from, to, family)]
ds <- sort(unique(R0$dataset))
message(sprintf("LODO over %d datasets", length(ds)))
lodo <- rbindlist(lapply(ds, function(x) edge_table(R0[dataset != x])[
    , .(from, to, family, partial_rho, n_donors, dropped_dataset = x)]))
fwrite(lodo, file.path(opt$outdir, "ras_network_lodo.tsv"), sep = "\t")
ls <- merge(lodo, edges[arm == "primary", .(from, to, family, full = partial_rho)],
            by = c("from", "to", "family"))
ls <- ls[, .(lodo_min = min(partial_rho, na.rm = TRUE), lodo_max = max(partial_rho, na.rm = TRUE),
             lodo_sign_consistency = mean(sign(partial_rho) == sign(full[1]), na.rm = TRUE),
             n_lodo = sum(is.finite(partial_rho))), by = .(from, to, family)]
edges <- merge(edges, bs, by = c("from", "to", "family"), all.x = TRUE)
edges <- merge(edges, ls, by = c("from", "to", "family"), all.x = TRUE)
edges[arm != "primary", c("boot_lo", "boot_hi", "lodo_min", "lodo_max",
                          "lodo_sign_consistency", "n_lodo") := NA]
edges[, n_datasets_lodo := length(ds)]
write_tsv_safe(edges, file.path(opt$outdir, "ras_network_edges.tsv"))
print(edges[arm == "primary", .(family, from, to, partial_rho, p_value, p_BH, n_donors,
                                boot_lo, boot_hi, lodo_sign_consistency)])

## ---- sequential chain + d-separation ------------------------------------------------------
chain <- rbindlist(lapply(names(DAG), function(ch) {
    cols <- c(ch, DAG[[ch]]); X <- R0[, ..cols]; X <- X[complete.cases(X)]
    fit <- lm(reformulate(DAG[[ch]], ch), data = X)
    co <- summary(fit)$coefficients[-1, , drop = FALSE]
    data.table(child = ch, parent = rownames(co), estimate = co[, 1], SE = co[, 2],
               p_value = co[, 4], n_donors = nrow(X),
               model = paste0("lm(", ch, " ~ ", paste(DAG[[ch]], collapse = " + "),
                              ") on study/depth-residualized nodes"))
}))
write_tsv_safe(chain, file.path(opt$outdir, "ras_network_chain_models.tsv"))
if (requireNamespace("dagitty", quietly = TRUE)) {
    g <- dagitty::dagitty(paste0("dag {", paste(sprintf("%s -> %s", ALLOWED$from, ALLOWED$to),
                                                collapse = " ; "), "}"))
    vars <- unique(c(ALLOWED$from, ALLOWED$to))
    X <- as.data.frame(R0[, ..vars]); X <- X[complete.cases(X), ]
    lt <- try(dagitty::localTests(g, data = X, type = "cis"), silent = TRUE)
    if (!inherits(lt, "try-error")) {
        lt <- as.data.table(lt, keep.rownames = "implied_independence")
        lt[, n_donors := nrow(X)]
        write_tsv_safe(lt, file.path(opt$outdir, "ras_network_dsep.tsv"))
    }
} else message("dagitty not installed; d-separation tests skipped")

cat("\nReproducibility information:\n"); print(sessionInfo())
