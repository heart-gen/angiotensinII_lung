"""
Prepare the CoGAPS input matrix from the HLCA pericyte subset.

CoGAPS (NMF-family Bayesian matrix factorization) is run here as a DATA-DRIVEN
validation / refinement of the curated 5-program, NVU-pattern pericyte-state
model (see pericyte_states/_h/00.state_discovery.py). The patterns CoGAPS learns
are unsupervised; downstream (02) we ask whether they corroborate the curated
programs and where the data diverge from them.

Scope: HLCA pericyte subset ONLY (~11.7k cells), not the full niche -- the full
niche would be too broad / expensive for a first pass.

Input is the already-state-annotated pericyte object so cell barcodes, stable
states, and curated scores line up with the validation step.

Exports a genes x cells matrix (CoGAPS factorizes genes x samples) as Matrix
Market sparse plus row/column names:
  - cogaps_input.mtx          (genes x cells, log-normalized, non-negative)
  - cogaps_genes.tsv          (gene symbols, row order)
  - cogaps_barcodes.tsv       (cell barcodes, column order)
  - cogaps_hvg_info.tsv.gz    (per-gene HVG / panel-membership flags)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse, io as sio

# Curated panels (mirror pericyte_states/_h/00.state_discovery.py) -- forced into
# the matrix even if not flagged HVG, so the validation step can map patterns to
# programs on a common gene space.
STATE_PANELS = {
    "vascular_stabilizing": ["RGS5", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "HIGD1B",
                             "COX4I2", "NDUFA4L2", "GJA4", "CSPG4", "PLXDC1"],
    "inflammatory": ["IL6", "CCL2", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CXCL6",
                     "CXCL8", "CXCL10", "IL1A", "IL1B", "MIF", "ICAM1", "VCAM1",
                     "SELE", "NFKBIA"],
    "synthetic_contractile": ["ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "DES",
                              "VIM", "PDGFRB"],
    "activated_migratory": ["ADAMTS1", "THBS1", "TIMP1", "MMP2", "MMP3", "MMP9",
                            "SERPINE1", "POSTN"],
    "fibroblast_like": ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "LUM",
                        "DCN", "PDGFA", "FBLN1"],
}
EXTRA_GENES = ["AGTR1", "AGTR2"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="State-annotated pericyte h5ad (pericyte_states.h5ad).")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--n-hvg", type=int, default=2000,
                   help="Number of highly variable genes to retain.")
    p.add_argument("--min-cells", type=int, default=10,
                   help="Drop genes detected in fewer than this many cells.")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def load_anndata(path: Path):
    adata = sc.read_h5ad(path)
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        if not np.array_equal(symbols.to_numpy(), adata.var_names.to_numpy()):
            adata.var["ensembl_id"] = adata.var_names
            adata.var_names = adata.var.index = symbols
    adata.raw = None
    adata.var_names_make_unique()
    return adata


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    np.random.seed(args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)
    adata.X = adata.layers["logcounts"]
    logging.info(f"loaded {adata.n_obs} pericytes x {adata.n_vars} genes")

    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    # HVGs on the log-normalized data; union with curated panels + receptors.
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=args.n_hvg)
    hvg = set(adata.var_names[adata.var["highly_variable"]].tolist())
    panel_genes = {g for genes in STATE_PANELS.values() for g in genes}
    panel_genes |= set(EXTRA_GENES)
    panel_present = {g for g in panel_genes if g in adata.var_names}
    missing = sorted(panel_genes - panel_present)
    if missing:
        logging.warning(f"panel genes absent from object: {missing}")

    keep = sorted(hvg | panel_present)
    sub = adata[:, keep].copy()
    logging.info(f"matrix: {len(keep)} genes ({len(hvg)} HVG, "
                 f"+{len(panel_present - hvg)} forced panel) x {sub.n_obs} cells")

    # genes x cells, non-negative log-normalized values
    X = sub.X
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)
    if X.min() < 0:
        raise ValueError("input has negative values; CoGAPS requires non-negative")
    gxc = sparse.csc_matrix(X.T)  # genes x cells

    sio.mmwrite(str(args.outdir / "cogaps_input.mtx"), gxc, field="real")
    pd.Series(sub.var_names).to_csv(args.outdir / "cogaps_genes.tsv",
                                    index=False, header=False)
    pd.Series(sub.obs_names).to_csv(args.outdir / "cogaps_barcodes.tsv",
                                    index=False, header=False)

    gene2panel = {}
    for prog, genes in STATE_PANELS.items():
        for g in genes:
            gene2panel.setdefault(g, []).append(prog)
    info = pd.DataFrame({"gene": sub.var_names})
    info["highly_variable"] = info["gene"].isin(hvg)
    info["panel"] = info["gene"].map(lambda g: ";".join(gene2panel.get(g, [])))
    info.to_csv(args.outdir / "cogaps_hvg_info.tsv.gz", sep="\t", index=False)

    logging.info(f"wrote genes x cells = {gxc.shape[0]} x {gxc.shape[1]}")
    session_info.show()


if __name__ == "__main__":
    main()
