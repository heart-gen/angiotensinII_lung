"""
Derive the in-vivo AngII response signature from McLellan et al. 2020
(ArrayExpress E-MTAB-8810). The rule is pre-specified in _h/signatures/README.md;
this script implements it and nothing else.

Outputs (to --outdir):
  external/E-MTAB-8810/{full_count_matrix.tsv, E-MTAB-8810.sdrf.txt, PROVENANCE.txt}
  angii_signature/cluster_labels.tsv          cluster -> label, marker scores
  angii_signature/cells_per_sample_type.tsv   QC of the pseudobulk units
  angii_signature/de_<celltype>_<contrast>.tsv full DESeq2 tables
  angii_signature/mouse_signature.tsv          selected genes (mouse symbols)
"""
import argparse
import hashlib
import logging
import subprocess
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

BASE = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/810/E-MTAB-8810/Files/"
FILES = ["full_count_matrix.tsv", "E-MTAB-8810.sdrf.txt"]

MARKERS = {
    "pericyte": ["Kcnj8", "Abcc9", "Rgs5", "Vtn", "Higd1b", "Pdgfrb", "Cspg4", "Notch3"],
    "smooth_muscle": ["Myh11", "Acta2", "Cnn1", "Tagln", "Des"],
    "fibroblast": ["Col1a1", "Pdgfra", "Dcn", "Gsn", "Tcf21"],
    "endothelial": ["Pecam1", "Cdh5", "Kdr", "Fabp4"],
    "macrophage": ["Adgre1", "Cd68", "Lyz2", "C1qa"],
    "cardiomyocyte": ["Myh6", "Tnnt2", "Ttn", "Actc1"],
}
ORDER = [("pericyte", ["pericyte"]),
         ("mural", ["pericyte", "smooth_muscle"]),
         ("fibroblast", ["fibroblast"])]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        n = r.headers.get("Content-Length")
    return int(n) if n else None


def fetch(dest: Path, tries=6):
    """Resumable download. full_count_matrix.tsv is ~1 GB and a plain urlretrieve
    truncated it at 382 MB on a compute node (2026-09-22), so transfers resume from
    the partial file (curl -C -) and the final size is checked against the server."""
    dest.mkdir(parents=True, exist_ok=True)
    lines = []
    for f in FILES:
        out = dest / f
        url = BASE + f
        want = remote_size(url)
        for k in range(1, tries + 1):
            have = out.stat().st_size if out.exists() else 0
            if want is not None and have == want:
                break
            logging.info(f"{f}: {have}/{want} bytes, attempt {k}/{tries}")
            subprocess.run(["curl", "-fLsS", "--retry", "5", "--retry-delay", "10",
                            "-C", "-", "-o", str(out), url], check=False)
        have = out.stat().st_size if out.exists() else 0
        if want is not None and have != want:
            raise IOError(f"{f}: got {have} of {want} bytes after {tries} attempts")
        logging.info(f"{f}: {have} bytes")
        lines.append(f"{f}\t{url}\tbytes={have}\tmd5={md5(out)}")
    (dest / "PROVENANCE.txt").write_text("\n".join(lines) + "\n")


def read_matrix(path: Path, chunk=2000):
    """genes x cells dense TSV -> cells x genes CSR, chunked to bound memory."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    # The header has one fewer field than the rows (no gene-column name).
    blocks, genes = [], []
    for df in pd.read_csv(path, sep="\t", header=None, skiprows=1, index_col=0,
                          chunksize=chunk, dtype={0: str}):
        genes.extend(df.index.astype(str))
        blocks.append(sparse.csr_matrix(df.to_numpy(dtype=np.int32)))
        logging.info(f"  read {len(genes)} genes")
    M = sparse.vstack(blocks).T.tocsr()
    if M.shape[0] != len(header):
        raise AssertionError(f"{M.shape[0]} cell columns vs {len(header)} header names")
    a = ad.AnnData(M, obs=pd.DataFrame(index=list(header)),
                   var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    return a


def deseq(counts, meta, ref, tag, outdir):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    try:
        dds = DeseqDataSet(counts=counts, metadata=meta, design="~condition", quiet=True)
    except TypeError:
        dds = DeseqDataSet(counts=counts, metadata=meta, design_factors="condition",
                           quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "AngII", ref], quiet=True)
    st.summary()
    res = st.results_df.copy()
    res.index.name = "mouse_gene"
    res.to_csv(outdir / f"de_{tag}.tsv", sep="\t")
    return res


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--outdir", type=Path, default=Path("."))
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=0.5)
    p.add_argument("--max-per-dir", type=int, default=200)
    p.add_argument("--min-genes", type=int, default=15)
    p.add_argument("--seed", type=int, default=13)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    ext = a.outdir / "external" / "E-MTAB-8810"
    out = a.outdir / "angii_signature"
    out.mkdir(parents=True, exist_ok=True)
    fetch(ext)

    sdrf = pd.read_csv(ext / "E-MTAB-8810.sdrf.txt", sep="\t")
    comp = sdrf[["Source Name", "Factor Value[compound]"]].drop_duplicates()
    treat = dict(zip(comp["Source Name"], comp["Factor Value[compound]"].fillna("none")))
    logging.info(f"sdrf treatments: {treat}")

    adata = read_matrix(ext / "full_count_matrix.tsv")
    suf = adata.obs_names.str.rsplit("_", n=1).str[-1]
    if set(suf) != {str(i) for i in range(2, 10)}:
        raise AssertionError(f"unexpected barcode suffixes {sorted(set(suf))}")
    adata.obs["sample"] = "AP1800" + suf
    adata.obs["treatment"] = adata.obs["sample"].map(treat).astype(str)
    adata.obs["condition"] = np.where(adata.obs["treatment"].str.contains("angiotensin"),
                                      "AngII", "control")
    logging.info(adata.obs.groupby(["sample", "treatment"]).size().to_string())

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
    adata = adata[(adata.obs["n_genes_by_counts"] >= 200) &
                  (adata.obs["pct_counts_mt"] < 30)].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, n_comps=30, mask_var="highly_variable", random_state=a.seed)
    sc.pp.neighbors(adata, n_pcs=30, random_state=a.seed)
    try:
        sc.tl.leiden(adata, resolution=1.0, random_state=a.seed, key_added="leiden",
                     flavor="igraph", n_iterations=2, directed=False)
    except (ImportError, ValueError, TypeError) as e:
        logging.warning(f"igraph leiden unavailable ({e}); using leidenalg")
        sc.tl.leiden(adata, resolution=1.0, random_state=a.seed, key_added="leiden")

    for lab, genes in MARKERS.items():
        g = [x for x in genes if x in adata.var_names]
        sc.tl.score_genes(adata, g, score_name=f"score_{lab}", random_state=a.seed)
    sc_cols = [f"score_{k}" for k in MARKERS]
    cm = adata.obs.groupby("leiden", observed=True)[sc_cols].mean()
    cz = (cm - cm.mean()) / cm.std(ddof=0)
    lab = cz.idxmax(axis=1).str.replace("score_", "", regex=False)
    cm["label"] = lab
    cm["n_cells"] = adata.obs["leiden"].value_counts()
    cm.to_csv(out / "cluster_labels.tsv", sep="\t")
    adata.obs["celltype"] = adata.obs["leiden"].map(lab).astype(str)
    logging.info("\n" + cm[["label", "n_cells"]].to_string())

    cps = adata.obs.groupby(["celltype", "sample"]).size().unstack(fill_value=0)
    cps.to_csv(out / "cells_per_sample_type.tsv", sep="\t")
    logging.info("\n" + cps.to_string())

    samples = sorted(adata.obs["sample"].unique())
    meta_all = (adata.obs[["sample", "treatment", "condition"]]
                .drop_duplicates().set_index("sample").loc[samples])

    chosen = None
    for tag, types in ORDER:
        m = adata.obs["celltype"].isin(types).to_numpy()
        rows, keep = [], []
        for s in samples:
            mm = m & (adata.obs["sample"] == s).to_numpy()
            if mm.sum() >= a.min_cells:
                rows.append(np.asarray(adata.layers["counts"][mm].sum(axis=0)).ravel())
                keep.append(s)
        meta = meta_all.loc[keep].copy()
        n_ang, n_ctl = (meta["condition"] == "AngII").sum(), (meta["condition"] == "control").sum()
        logging.info(f"[{tag}] units: {len(keep)} (AngII {n_ang}, control {n_ctl})")
        if n_ang < 3 or n_ctl < 3:
            logging.warning(f"[{tag}] fewer than 3 units per arm; next cell type")
            continue
        counts = pd.DataFrame(np.vstack(rows), index=keep, columns=adata.var_names)
        counts = counts.loc[:, counts.sum(axis=0) >= 10]
        res = deseq(counts, meta[["condition"]], "control", f"{tag}_AngII_vs_control", out)
        # sensitivity: AngII vs saline only (n = 2 saline) -- sign concordance only
        sal = meta["treatment"].isin(["saline"]) | (meta["condition"] == "AngII")
        if (meta.loc[sal, "condition"] == "control").sum() >= 2:
            m2 = meta.loc[sal, ["condition"]].copy()
            deseq(counts.loc[m2.index], m2, "control", f"{tag}_AngII_vs_saline", out)
        sig = res[(res["padj"] < a.padj) & (res["log2FoldChange"].abs() >= a.lfc)].copy()
        sig["direction"] = np.where(sig["log2FoldChange"] > 0, "up", "down")
        sig = (sig.sort_values("padj").groupby("direction", group_keys=False)
                  .head(a.max_per_dir))
        logging.info(f"[{tag}] passing genes: {len(sig)} "
                     f"({(sig['direction'] == 'up').sum()} up)")
        # The >= min_genes floor is applied AFTER ortholog mapping and disjointness
        # pruning (README step 6); a generous pre-check avoids a doomed cell type.
        if len(sig) >= a.min_genes * 2 or chosen is None:
            chosen = (tag, sig)
        if len(sig) >= a.min_genes * 2:
            break
    if chosen is None:
        raise RuntimeError("no cell type had >= 3 units per arm")
    tag, sig = chosen
    sig = sig.reset_index()[["mouse_gene", "direction", "log2FoldChange", "padj", "baseMean"]]
    sig["cell_type"] = tag
    sig["source_id"] = "E-MTAB-8810 (McLellan 2020, doi:10.1161/CIRCULATIONAHA.119.045115)"
    sig["contrast"] = "AngII (n=4) vs saline+untreated (n=4), pseudobulk DESeq2"
    sig.to_csv(out / "mouse_signature.tsv", sep="\t", index=False)
    logging.info(f"signature cell type: {tag}; {len(sig)} mouse genes")


if __name__ == "__main__":
    main()
