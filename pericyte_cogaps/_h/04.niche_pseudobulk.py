"""
Pseudobulk the CCC niche onto the CoGAPS gene space for projectR transfer.

The pericyte CoGAPS patterns (pericyte_cogaps/) were learned on pericytes ONLY.
To ask whether those latent programs -- in particular the vascular-stabilizing
and the injury (inflammatory / fibroblast-like) patterns -- are MIRRORED in the
neighbouring niche cell types, we project the learned gene loadings onto the rest
of the niche with projectR (05.project_niche.R). projectR is a linear projection,
so it is run on a clean, donor-aware granularity rather than 423k raw cells:
cell-type x donor PSEUDOBULK (mean log-normalized expression). This both bounds
memory and gives the projected pattern usage a natural unit (donor) for the
downstream mixed-model test.

Scope note: this is the deliberate, targeted "expand to the niche" step. We do NOT
re-run CoGAPS on the whole niche (that would mostly relearn cell-type identity);
we TRANSFER the validated pericyte patterns outward, which is projectR's intended
use and directly tests the coordinated-niche-program hypothesis.

Outputs (to ../_m):
  - niche_pseudobulk_logmean.tsv.gz   genes x (cell_type||donor) samples
  - niche_pseudobulk_samples.tsv      sample -> cell_type, donor_id, disease_group, n_cells
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--niche", type=Path,
                   default=Path("../../cell_communication/_m/ccc_niche.h5ad"))
    p.add_argument("--genes", type=Path, default=Path("../_m/cogaps_genes.tsv"),
                   help="CoGAPS gene order (one symbol per line)")
    p.add_argument("--outdir", type=Path, default=Path("../_m"))
    p.add_argument("--min-cells", type=int, default=10,
                   help="Min cells per (cell_type x donor) sample to keep")
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(args.genes, header=None)[0].astype(str).tolist()
    logging.info(f"CoGAPS gene space: {len(genes)} genes")

    adata = sc.read_h5ad(args.niche)
    if "logcounts" in adata.layers:
        adata.X = adata.layers["logcounts"]
    adata.var_names_make_unique()

    shared = [g for g in genes if g in set(adata.var_names)]
    logging.info(f"Shared genes niche∩CoGAPS: {len(shared)}/{len(genes)}")
    adata = adata[:, shared].copy()

    obs = adata.obs
    sample = (obs["cell_type"].astype(str) + "||" + obs["donor_id"].astype(str))
    sample = sample.to_numpy()
    codes, uniq = pd.factorize(sample)

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    # mean log-normalized expression per sample = (group-sum) / (group-n)
    n_samp = len(uniq)
    onehot = sparse.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)),
                               shape=(len(codes), n_samp))
    sums = onehot.T @ X                      # samples x genes
    counts = np.asarray(onehot.sum(axis=0)).ravel()  # cells per sample
    means = np.asarray(sums.todense()) / counts[:, None]

    pb = pd.DataFrame(means, index=uniq, columns=shared)  # samples x genes
    meta = pd.DataFrame({"sample": uniq})
    meta[["cell_type", "donor_id"]] = meta["sample"].str.split(r"\|\|", n=1, expand=True)
    # disease_group is donor-level; take first per sample
    dmap = (obs.assign(_s=sample)[["_s", "disease_group"]]
            .drop_duplicates("_s").set_index("_s")["disease_group"])
    meta["disease_group"] = meta["sample"].map(dmap).astype(str)
    meta["n_cells"] = counts.astype(int)

    keep = meta["n_cells"] >= args.min_cells
    logging.info(f"Samples: {n_samp} total, {keep.sum()} with >= {args.min_cells} cells")
    pb, meta = pb.loc[keep.to_numpy()], meta.loc[keep.to_numpy()]

    # write genes x samples (projectR factorizes genes x samples, like CoGAPS input)
    pb.T.to_csv(args.outdir / "niche_pseudobulk_logmean.tsv.gz", sep="\t")
    meta.to_csv(args.outdir / "niche_pseudobulk_samples.tsv", sep="\t", index=False)
    logging.info(f"Wrote pseudobulk {pb.T.shape} (genes x samples) and {len(meta)} sample rows")
    session_info.show()


if __name__ == "__main__":
    main()
