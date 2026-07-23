"""
Integrate the mouse lung niche across CELLxGENE datasets with scVI.

The downloaded mouse niche pools 4 datasets of very unequal size; without
integration, global marker z-scoring is dominated by the largest dataset and
state assignment is confounded with batch. This step removes dataset batch
effects with scVI (matching the human-side scVI/SCANVI workflow used in this
repo) so that state assignment can be made on an integrated representation.

Pipeline (mirrors disease_association/pericyte_analysis/_h/02.train_model.py and
the nvu-neuroimmune-reference 02_integration scVI utilities):
  raw counts -> batch-aware seurat_v3 HVG (force-include state markers) ->
  scVI(batch_key=dataset_id, n_latent=30, n_layers=2) -> latent + Leiden + UMAP +
  batch-corrected expression layer (transform_batch = largest dataset).

Output: mouse_integrated.h5ad (HVG genes; layers counts/lognorm/scvi_corrected;
obsm X_scVI/X_umap; obs leiden_scvi).
"""
import numpy as np
import scanpy as sc
import scvi
import session_info
import logging, argparse
from pathlib import Path

# State marker programs (mouse orthologs) + receptors -- force-kept in the HVG
# set so state scoring and receptor readout use integrated genes.
STATE_PANELS_MOUSE = {
    "vascular_stabilizing": ["Rgs5", "Pdgfrb", "Notch3", "Kcnj8", "Abcc9",
                             "Higd1b", "Cox4i2", "Ndufa4l2", "Gja4", "Cspg4"],
    "inflammatory": ["Il6", "Ccl2", "Ccl20", "Cxcl1", "Cxcl2", "Cxcl3",
                     "Il1a", "Il1b", "Mif", "Icam1", "Vcam1", "Sele", "Nfkbia"],
    "synthetic_contractile": ["Acta2", "Myh11", "Tagln", "Cnn1", "Myl9", "Des", "Vim"],
    "activated_migratory": ["Adamts1", "Thbs1", "Timp1", "Mmp2", "Mmp3", "Serpine1", "Postn"],
    "fibroblast_like": ["Col1a1", "Col1a2", "Col3a1", "Col4a1", "Fn1", "Lum", "Dcn", "Pdgfa"],
    # Basement-membrane deposition, mirroring the human panel in
    # basement_membrane/_h/bm_panels.py (all 13 orthologs are 1:1 and are present
    # in mouse_lung.h5ad). Kept separate from fibroblast_like because the two
    # matrices are biologically distinct and, in human, near-orthogonal once the
    # shared Col4a1 is removed.
    "basement_membrane": ["Col4a1", "Col4a2", "Col18a1", "Lama3", "Lama4",
                          "Lama5", "Lamb1", "Lamb2", "Lamc1", "Nid1", "Nid2",
                          "Hspg2", "Agrn"],
}
RECEPTORS = ["Agtr1a", "Agtr1b", "Agtr2"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="mouse_lung.h5ad (raw counts)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--batch-key", default="dataset_id")
    p.add_argument("--n-top-genes", type=int, default=3000)
    p.add_argument("--n-latent", type=int, default=30)
    p.add_argument("--max-epochs", type=int, default=400)
    p.add_argument("--leiden-resolution", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    scvi.settings.seed = args.seed
    sc.settings.verbosity = 1

    adata = sc.read_h5ad(args.adata)
    # Census X is raw counts
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()

    # Batch-aware HVG (skip batch if any dataset too small), then force-include markers
    bkey = args.batch_key
    counts_ok = adata.obs[bkey].value_counts().min() >= 100
    sc.pp.highly_variable_genes(
        adata, n_top_genes=args.n_top_genes, flavor="seurat_v3", layer="counts",
        batch_key=bkey if counts_ok else None, subset=False)
    force = set(g for panel in STATE_PANELS_MOUSE.values() for g in panel) | set(RECEPTORS)
    in_obj = [g for g in force if g in adata.var_names]
    adata.var.loc[adata.var_names.isin(in_obj), "highly_variable"] = True
    logging.info(f"forced {len(in_obj)} marker/receptor genes into HVG set")

    ref = adata[:, adata.var["highly_variable"]].copy()
    logging.info(f"scVI input: {ref.shape}; batches={ref.obs[bkey].nunique()}")

    scvi.model.SCVI.setup_anndata(ref, layer="counts", batch_key=bkey)
    vae = scvi.model.SCVI(ref, n_latent=args.n_latent, n_layers=2)
    vae.train(max_epochs=args.max_epochs, early_stopping=True,
              early_stopping_patience=15, accelerator="auto")

    ref.obsm["X_scVI"] = vae.get_latent_representation()
    # Batch-corrected expression (transform to the largest dataset as reference)
    ref_batch = ref.obs[bkey].value_counts().idxmax()
    corr = vae.get_normalized_expression(library_size=1e4, transform_batch=ref_batch,
                                         return_numpy=True)
    ref.layers["scvi_corrected"] = np.log1p(corr)

    # Integrated neighbors / clustering / embedding
    sc.pp.neighbors(ref, use_rep="X_scVI", random_state=args.seed)
    sc.tl.leiden(ref, resolution=args.leiden_resolution, key_added="leiden_scvi",
                 random_state=args.seed)
    sc.tl.umap(ref, random_state=args.seed)

    vae.save(str(args.outdir / "scvi_model"), overwrite=True, save_anndata=False)
    if ref.var.index.name in ref.var.columns:
        ref.var.index.name = None
    ref.write(args.outdir / "mouse_integrated.h5ad")
    logging.info(f"wrote mouse_integrated.h5ad: {ref.shape}; "
                 f"leiden clusters={ref.obs['leiden_scvi'].nunique()}")
    session_info.show()


if __name__ == "__main__":
    main()
