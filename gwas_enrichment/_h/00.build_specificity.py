#!/usr/bin/env python
"""
00.build_specificity.py — gene-level cell-type / pericyte-state specificity for
GWAS heritability enrichment (MAGMA gene-property + S-LDSC annotations).

Two specificity matrices (EWCE/MAGMA-Celltyping convention: spec_g = mean_g / sum_g mean):
  (1) across HLCA coarse lung cell types  (places the AGTR1+ pericyte among all lung cells)
  (2) within pericytes, across the 5 functional states (the headline: does disease
      heritability favour the injury-niche state over the vascular-stabilizing state?)

X-chromosome handling (the AGTR2 constraint): specificity is restricted to **autosomal
protein-coding genes** mapped to Entrez IDs via MAGMA's NCBI37.3.gene.loc (hg19). chrX/Y
genes — including AGTR2 (Xq23) — are dropped so the downstream SNP-heritability tests stay
on the autosomal LD reference. AGTR1 (chr3) is retained and flagged in the QC log.

Outputs (gwas_enrichment/_m/specificity/):
  specificity_celltype.tsv        Entrez + one column per coarse cell type
  specificity_pericyte_state.tsv  Entrez + one column per pericyte state
  specificity_all.tsv             union of covariate columns (MAGMA --gene-covar input)
  gene_coords_hg19.tsv            Entrez, symbol, chr, start, end (autosomal)
  beds/<annotation>.bed           top-decile genes per annotation ±100kb (hg19, for LDSC)
  qc_specificity.txt              gene counts, dropped X/Y, AGTR1/AGTR2 status
"""
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import session_info

WINDOW = 100_000   # ±100 kb gene window for LDSC annotations (LDSC-SEG convention)
TOP_FRAC = 0.10    # top decile of specificity defines each annotation's gene set

GENE_LOC = "/ocean/projects/bio250020p/shared/opt/magma-v1.10/NCBI37.3.gene.loc"


def load_gene_loc():
    """NCBI37.3.gene.loc: entrez, chr, start, end, strand, symbol (hg19)."""
    g = pd.read_csv(GENE_LOC, sep="\t", header=None,
                    names=["entrez", "chr", "start", "end", "strand", "symbol"],
                    dtype={"entrez": str, "chr": str})
    g = g[g["chr"].isin([str(i) for i in range(1, 23)])].copy()   # autosomes only
    g["chr"] = g["chr"].astype(int)
    # one Entrez per symbol (first occurrence) for the symbol->entrez join
    g = g.drop_duplicates(subset="symbol", keep="first")
    return g


def group_mean(adata, key):
    """Mean expression per group via one-hot matmul; returns DataFrame genes x groups."""
    labels = adata.obs[key].astype("category")
    groups = labels.cat.categories
    onehot = sp.csr_matrix(
        (np.ones(adata.n_obs),
         (np.arange(adata.n_obs), labels.cat.codes.values)),
        shape=(adata.n_obs, len(groups)))
    counts = np.asarray(onehot.sum(axis=0)).ravel()
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    sums = onehot.T @ X                      # groups x genes
    means = np.asarray(sums.todense()) / counts[:, None]
    return pd.DataFrame(means.T, index=adata.var_names, columns=list(groups))


def specificity(mean_df):
    """spec = group mean / row-sum across groups; drop all-zero genes."""
    rowsum = mean_df.sum(axis=1)
    keep = rowsum > 0
    spec = mean_df.loc[keep].div(rowsum[keep], axis=0)
    return spec


def to_symbol_index(adata):
    """Use feature_name (gene symbol) as var index, dedup."""
    if "feature_name" in adata.var.columns:
        adata.var["symbol"] = adata.var["feature_name"].astype(str)
    else:
        adata.var["symbol"] = adata.var_names.astype(str)
    adata.var_names = adata.var["symbol"].values
    adata.var_names_make_unique()
    return adata


def map_to_entrez(spec, gloc, qc):
    """Join symbol-indexed specificity to autosomal Entrez coordinates."""
    df = spec.copy()
    df.index.name = "symbol"
    df = df.reset_index()
    n_in = df.shape[0]
    merged = df.merge(gloc[["entrez", "symbol", "chr", "start", "end"]],
                      on="symbol", how="inner")
    qc.append(f"  {n_in} genes scored -> {merged.shape[0]} mapped to autosomal Entrez")
    return merged


def write_beds(spec_entrez, annot_cols, gloc, outdir, qc):
    beddir = os.path.join(outdir, "beds")
    os.makedirs(beddir, exist_ok=True)
    coords = gloc.set_index("entrez")[["chr", "start", "end"]]
    for col in annot_cols:
        s = spec_entrez[["entrez", col]].dropna()
        thr = s[col].quantile(1 - TOP_FRAC)
        top = s[s[col] >= thr]["entrez"].astype(str).tolist()
        c = coords.loc[[e for e in top if e in coords.index]].copy()
        c["bstart"] = np.maximum(0, c["start"] - WINDOW)
        c["bend"] = c["end"] + WINDOW
        bed = c[["chr", "bstart", "bend"]].copy()
        bed["chr"] = "chr" + bed["chr"].astype(str)
        safe = col.replace(" ", "_").replace("/", "_")
        bed.sort_values(["chr", "bstart"]).to_csv(
            os.path.join(beddir, f"{safe}_hg19.bed"), sep="\t", header=False, index=False)
        qc.append(f"  BED {safe}: {bed.shape[0]} genes (top {int(TOP_FRAC*100)}%)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hlca", default="../../inputs/hlca/_m/hlca_core.h5ad")
    ap.add_argument("--pericyte", default="../../pericyte_states/_m/pericyte_states.h5ad")
    ap.add_argument("--celltype-key", default="ann_coarse_for_GWAS_and_modeling")
    ap.add_argument("--state-key", default="pericyte_state")
    ap.add_argument("--outdir", default="../_m")
    args = ap.parse_args()

    outdir = os.path.join(args.outdir, "specificity")
    os.makedirs(outdir, exist_ok=True)
    qc = ["== build_specificity QC =="]
    gloc = load_gene_loc()
    qc.append(f"gene.loc autosomal genes: {gloc.shape[0]}")

    # (1) cross-cell-type specificity from HLCA core
    print("Loading HLCA core ...", flush=True)
    hl = sc.read_h5ad(args.hlca)
    hl = to_symbol_index(hl)
    key = args.celltype_key if args.celltype_key in hl.obs.columns else "cell_type"
    qc.append(f"cell-type key: {key} ({hl.obs[key].nunique()} groups)")
    ct_spec = specificity(group_mean(hl, key))
    ct_entrez = map_to_entrez(ct_spec, gloc, qc)
    ct_cols = [c for c in ct_spec.columns]
    del hl

    # (2) within-pericyte state specificity
    print("Loading pericyte states ...", flush=True)
    pc = sc.read_h5ad(args.pericyte)
    pc = to_symbol_index(pc)
    skey = args.state_key
    qc.append(f"state key: {skey} ({pc.obs[skey].nunique()} states)")
    st_spec = specificity(group_mean(pc, skey))
    st_entrez = map_to_entrez(st_spec, gloc, qc)
    st_cols = [c for c in st_spec.columns]

    # AGTR1 / AGTR2 status (the X-chromosome caveat, made explicit)
    for gene in ("AGTR1", "AGTR2"):
        row = gloc[gloc["symbol"] == gene]
        if row.empty:
            qc.append(f"  {gene}: NOT in autosomal gene.loc (expected for AGTR2/chrX)")
        else:
            qc.append(f"  {gene}: chr{row['chr'].iloc[0]} -> retained")

    # Write specificity tables (Entrez-keyed for MAGMA --gene-covar)
    ct_out = ct_entrez[["entrez"] + ct_cols].rename(columns={"entrez": "GENE"})
    st_out = st_entrez[["entrez"] + st_cols].rename(columns={"entrez": "GENE"})
    ct_out.to_csv(os.path.join(outdir, "specificity_celltype.tsv"), sep="\t", index=False)
    st_out.to_csv(os.path.join(outdir, "specificity_pericyte_state.tsv"), sep="\t", index=False)
    # Combined covariate file (outer join on Entrez GENE)
    allspec = ct_out.merge(st_out, on="GENE", how="outer")
    allspec.to_csv(os.path.join(outdir, "specificity_all.tsv"), sep="\t", index=False)

    gloc[["entrez", "symbol", "chr", "start", "end"]].to_csv(
        os.path.join(outdir, "gene_coords_hg19.tsv"), sep="\t", index=False)

    # BEDs for LDSC (cell types + pericyte states)
    write_beds(ct_entrez, ct_cols, gloc, outdir, qc)
    write_beds(st_entrez, st_cols, gloc, outdir, qc)

    with open(os.path.join(outdir, "qc_specificity.txt"), "w") as fh:
        fh.write("\n".join(qc) + "\n")
    print("\n".join(qc))
    print("\n[done] specificity written to", outdir)
    session_info.show()


if __name__ == "__main__":
    main()
