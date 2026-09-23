## Mouse -> human, strict one-to-one, for the AngII signature (README step 7), and
## freeze the result into _h/signatures/ if no frozen copy exists yet.
suppressPackageStartupMessages({
    .libPaths(c("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung/.Rlib",
                .libPaths()))
    library(babelgene); library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
inp    <- if (length(args) >= 1) args[1] else "./angii_signature/mouse_signature.tsv"
frozen <- if (length(args) >= 2) args[2] else "../_h/signatures/angii_response_signature.tsv"

sig <- fread(inp)
orth <- as.data.table(orthologs(genes = sig$mouse_gene, species = "mouse", human = FALSE))
orth <- unique(orth[, .(mouse_gene = symbol, gene = human_symbol, support_n)])
## strict 1:1 -- each mouse gene to exactly one human gene, and vice versa
one <- orth[, .N, by = mouse_gene][N == 1, mouse_gene]
orth <- orth[mouse_gene %in% one]
back <- orth[, .N, by = gene][N == 1, gene]
orth <- orth[gene %in% back]
hs <- merge(sig, orth, by = "mouse_gene")
message(sprintf("orthologs: %d of %d mouse genes mapped 1:1", nrow(hs), nrow(sig)))
setcolorder(hs, c("gene", "direction", "log2FoldChange", "padj", "mouse_gene"))
fwrite(hs, sub("mouse_signature", "human_signature", inp), sep = "\t")

if (file.exists(frozen)) {
    message("frozen signature exists, NOT overwritten: ", frozen)
} else {
    hs[, frozen_on := format(Sys.Date())]
    fwrite(hs, frozen, sep = "\t")
    message("FROZE signature: ", frozen, " (", nrow(hs), " genes)")
}
