## Figure 5F -- NicheNet with the pericyte as SENDER.
##
## cell_communication/_h/02.nichenet.R asks which niche ligands best explain a
## pericyte target programme. This inverts it: for each alveolar-capillary
## neighbour with a pre-specified RESPONSE programme (receiver_programs.py), which
## PERICYTE-expressed ligands best explain that programme in that neighbour?
##
## Setup mirrors 02.nichenet.R exactly (same priors, same 0.10 expression
## threshold, same receptor filter, same skip rule), with the sender restricted to
## Pericytes. Ligands are ordered by permutation z against random target sets of
## the same size, as 02b.nichenet_specificity.R does -- never by AUPR or by the
## floored empirical p (memory: nichenet-rank-by-z).
##
## Circularity guard: a programme gene that pericytes could themselves be scored
## for as a ligand is DROPPED from that programme and written to
## program_genes_dropped.tsv. The drop is visible, not silent.
##
## Also writes the gene list and design for the donor-level validation in
## 10.outgoing_donor_validation.R: the top K displayable ligands per target, and
## for each a pool of detection-matched pericyte-expressed NON-ligand genes.

suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(nichenetr); library(dplyr); library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default) {
    i <- which(args == flag); if (length(i)) args[i + 1] else default
}
PRIORS    <- parse_arg("--priors", "../../cell_communication/_m/nichenet_priors")
FRAC_FILE <- parse_arg("--frac-file",
                       "../../cell_communication/_m/expressed_fraction_per_group.tsv.gz")
PROGRAMS  <- parse_arg("--programs", "./receiver_programs.tsv")
TARGETS   <- parse_arg("--targets", "./outgoing_targets.tsv")
OUTDIR    <- parse_arg("--outdir", "./nichenet_outgoing")
EXPR_THR  <- as.numeric(parse_arg("--expr-thr", "0.10"))
N_PERM    <- as.integer(parse_arg("--n-perm", "1000"))
TOPK      <- as.integer(parse_arg("--top-k", "5"))
POOL_TOL  <- as.numeric(parse_arg("--pool-tol", "0.02"))
POOL_MAX  <- as.integer(parse_arg("--pool-max", "30"))
SEED      <- as.integer(parse_arg("--seed", "13"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

SENDER <- "Pericytes"
NON_LIGANDS <- c("COPA", "MMP14", "SIRPB2")   # figures/_h/_fig_common.R
SUBSTRATES  <- c("AGT")                       # renin substrate, not the AT1R ligand
RAS_GENES   <- c("AGT", "REN", "ACE", "ACE2", "CMA1", "CTSG", "CTSD", "ENPEP", "MME",
                 "AGTR1", "AGTR2", "LRP2", "MAS1")

ligand_target_matrix <- readRDS(file.path(PRIORS, "ligand_target_matrix_nsga2r_final.rds"))
lr_network <- readRDS(file.path(PRIORS, "lr_network_human_21122021.rds")) |>
    distinct(from, to)
all_ligands <- unique(lr_network$from); all_receptors <- unique(lr_network$to)

frac <- fread(FRAC_FILE, data.table = FALSE)
rownames(frac) <- frac[[1]]; frac[[1]] <- NULL
## AT2 is split upstream on AGTR2 detectability. For "expressed in AT2" take the
## larger of the two strata: a gene expressed in either stratum is expressed by AT2.
if (all(c("AT2_AGTR2det", "AT2_AGTR2undet") %in% colnames(frac)))
    frac$AT2 <- pmax(frac$AT2_AGTR2det, frac$AT2_AGTR2undet)
expressed_in <- function(g) if (g %in% colnames(frac)) rownames(frac)[frac[[g]] >= EXPR_THR] else character(0)

prog <- fread(PROGRAMS); tg <- fread(TARGETS)
## Collapse the two AT2 strata to one NicheNet target (programmes are per cell type).
tg[, nn_target := fifelse(grepl("^AT2_", target), "AT2", target)]
tg <- unique(tg[nzchar(program), .(nn_target, program)])
stopifnot(SENDER %in% colnames(frac))

sender_expr <- expressed_in(SENDER)
sender_ligands <- intersect(all_ligands, sender_expr)
message(sprintf("sender %s: %d expressed genes, %d expressed prior-network ligands",
                SENDER, length(sender_expr), length(sender_ligands)))

dropped <- list(); summ <- list(); design <- list()
for (i in seq_len(nrow(tg))) {
    tgt <- tg$nn_target[i]; pname <- tg$program[i]
    safe <- gsub("[^A-Za-z0-9]+", "_", tgt)
    message("=== ", SENDER, " -> ", tgt, " [", pname, "] ===")
    recv_expr <- expressed_in(tgt)
    if (length(recv_expr) < 50) { message("  too few expressed genes; skipping"); next }
    background <- intersect(recv_expr, rownames(ligand_target_matrix))
    potential <- lr_network |>
        filter(from %in% sender_ligands, to %in% intersect(all_receptors, recv_expr)) |>
        pull(from) |> unique() |> intersect(colnames(ligand_target_matrix))

    pg <- prog[program == pname, gene]
    in_bg <- intersect(pg, background)
    clash <- intersect(in_bg, c(sender_ligands, RAS_GENES))
    if (length(clash))
        dropped[[length(dropped) + 1]] <- data.table(target = tgt, program = pname,
                                                     gene = clash,
                                                     reason = "also a pericyte-expressed prior-network ligand or RAS gene")
    geneset <- setdiff(in_bg, clash)
    message(sprintf("  background %d | programme %d of %d expressed, %d dropped | potential ligands %d",
                    length(background), length(geneset), length(pg), length(clash),
                    length(potential)))
    if (length(potential) < 3 || length(geneset) < 5) {
        message("  insufficient ligands/targets; skipping")
        summ[[length(summ) + 1]] <- data.table(target = tgt, program = pname,
            n_background = length(background), n_geneset = length(geneset),
            n_potential = length(potential), status = "skipped (<3 ligands or <5 genes)")
        next
    }

    act <- predict_ligand_activities(
        geneset = geneset, background_expressed_genes = background,
        ligand_target_matrix = ligand_target_matrix, potential_ligands = potential) |>
        as.data.table()
    setorder(act, -aupr_corrected)
    act[, `:=`(rank_by_aupr = seq_len(.N), target = tgt, program = pname,
               not_a_ligand = test_ligand %in% NON_LIGANDS,
               substrate_not_ligand = test_ligand %in% SUBSTRATES)]
    fwrite(act, file.path(OUTDIR, paste0("ligand_activities_", safe, ".tsv")), sep = "\t")

    ## permutation null over the top ligands (AUPR per ligand does not depend on
    ## which other ligands are scored, so restricting is exact and faster)
    top <- head(act$test_ligand, 25)
    k <- length(geneset)
    nm <- matrix(NA_real_, length(top), N_PERM, dimnames = list(top, NULL))
    for (b in seq_len(N_PERM)) {
        av <- predict_ligand_activities(
            geneset = sample(background, k), background_expressed_genes = background,
            ligand_target_matrix = ligand_target_matrix, potential_ligands = top)
        nm[av$test_ligand, b] <- av$aupr_corrected
        if (b %% 250 == 0) message("  perm ", b, "/", N_PERM)
    }
    sp <- act[test_ligand %in% top]
    sp[, `:=`(null_mean = rowMeans(nm[test_ligand, , drop = FALSE], na.rm = TRUE),
              null_sd = apply(nm[test_ligand, , drop = FALSE], 1, sd, na.rm = TRUE))]
    sp[, z := (aupr_corrected - null_mean) / null_sd]
    sp[, p_emp := sapply(seq_len(.N), function(j)
        (1 + sum(nm[test_ligand[j], ] >= aupr_corrected[j], na.rm = TRUE)) / (N_PERM + 1))]
    sp[, `:=`(p_emp_at_floor = p_emp <= 1 / (N_PERM + 1), n_perm = N_PERM,
              rank_by_z = rank(-z, ties.method = "first"),
              p_emp_BH = p.adjust(p_emp, "BH"))]
    setorder(sp, rank_by_z)
    fwrite(sp, file.path(OUTDIR, paste0("specificity_", safe, ".tsv")), sep = "\t")

    disp <- sp[!not_a_ligand & !substrate_not_ligand]
    best <- head(disp$test_ligand, 20)
    links <- rbindlist(lapply(best, function(l)
        as.data.table(get_weighted_ligand_target_links(
            l, geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 200))),
        fill = TRUE)
    if (nrow(links)) links[, target_cell := tgt]
    fwrite(links, file.path(OUTDIR, paste0("links_", safe, ".tsv")), sep = "\t")

    summ[[length(summ) + 1]] <- data.table(target = tgt, program = pname,
        n_background = length(background), n_geneset = length(geneset),
        n_potential = length(potential), status = "run",
        top_ligands_by_z = paste(head(disp$test_ligand, 5), collapse = ";"),
        top_z = paste(sprintf("%.1f", head(disp$z, 5)), collapse = ";"))

    ## ---- validation design: top-K ligands + detection-matched non-ligand pools
    chosen <- head(disp$test_ligand, TOPK)
    pd <- frac[[SENDER]]; names(pd) <- rownames(frac)
    eligible <- setdiff(names(pd)[pd > 0], c(all_ligands, RAS_GENES, prog$gene))
    for (l in chosen) {
        d <- abs(pd[eligible] - pd[[l]])
        cand <- names(d)[d <= POOL_TOL]
        if (length(cand) > POOL_MAX) cand <- sample(cand, POOL_MAX)
        design[[length(design) + 1]] <- rbind(
            data.table(target = tgt, program = pname, role = "ligand", gene = l,
                       matched_to = l, pericyte_detect = pd[[l]]),
            if (length(cand)) data.table(target = tgt, program = pname, role = "null_pool",
                                         gene = cand, matched_to = l,
                                         pericyte_detect = unname(pd[cand])))
    }
    design[[length(design) + 1]] <- data.table(target = tgt, program = pname,
        role = "program", gene = geneset, matched_to = NA_character_,
        pericyte_detect = NA_real_)
}

fwrite(rbindlist(dropped, fill = TRUE), file.path(OUTDIR, "program_genes_dropped.tsv"), sep = "\t")
fwrite(rbindlist(summ, fill = TRUE), file.path(OUTDIR, "outgoing_nichenet_summary.tsv"), sep = "\t")
des <- rbindlist(design, fill = TRUE)
fwrite(des, file.path(OUTDIR, "outgoing_validation_design.tsv"), sep = "\t")
fwrite(data.table(gene = sort(unique(des$gene))), "./outgoing_pb_genes.tsv", sep = "\t")
message(sprintf("validation design: %d targets, %d unique genes for pseudobulk",
                uniqueN(des$target), uniqueN(des$gene)))
cat("\nReproducibility information:\n"); print(sessionInfo())
