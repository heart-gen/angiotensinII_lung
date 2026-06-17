"""
Download a mouse lung stromal/vascular niche from CELLxGENE Census for
cross-species conservation of pericyte states + angiotensin receptors.

NOTE: annotated mouse lung "pericyte" cells are sparse in the Census (~99),
so we download the broader mural/vascular/stromal niche (pericytes + vascular
SMC + fibroblasts + endothelium + AT1/AT2) and focus conservation on the
mural/pericyte-like compartment. Run on a login/data node (Census needs
internet; PSC compute nodes do not have outbound network). Cell count is capped
to keep memory modest.
"""
import numpy as np
import pandas as pd
import logging, argparse
from pathlib import Path

NICHE_MOUSE = [
    "pericyte", "vascular associated smooth muscle cell",
    "aortic smooth muscle cell", "bronchial smooth muscle cell",
    "smooth muscle cell of the pulmonary artery",
    "adventitial fibroblast", "alveolar adventitial fibroblast",
    "fibroblast of lung", "alveolar type 1 fibroblast cell",
    "pulmonary interstitial fibroblast", "mesenchymal cell",
    "capillary endothelial cell", "alveolar capillary type 1 endothelial cell",
    "alveolar capillary type 2 endothelial cell", "endothelial cell of artery",
    "vein endothelial cell",
    "pulmonary alveolar type 1 cell", "pulmonary alveolar type 2 cell",
]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--census-version", default="stable")
    p.add_argument("--max-cells", type=int, default=60000)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    import cellxgene_census

    types_str = ", ".join(f"'{t}'" for t in NICHE_MOUSE)
    vf = (f"tissue_general=='lung' and is_primary_data==True "
          f"and cell_type in [{types_str}]")

    with cellxgene_census.open_soma(census_version=args.census_version) as census:
        obs = cellxgene_census.get_obs(
            census, "Mus musculus", value_filter=vf,
            column_names=["soma_joinid", "cell_type", "dataset_id", "donor_id",
                          "sex", "disease", "assay"])
        logging.info(f"Niche cells available: {len(obs)}")
        obs.to_csv(args.outdir / "mouse_lung_niche_obs_inventory.tsv.gz",
                   sep="\t", index=False)

        # Cap by stratified subsample per cell_type (observed groups only)
        if len(obs) > args.max_cells:
            obs = obs.copy()
            obs["cell_type"] = obs["cell_type"].astype(str)
            frac = args.max_cells / len(obs)
            obs = (obs.groupby("cell_type", group_keys=False, observed=True)
                      .apply(lambda d: d.sample(max(1, int(round(len(d) * frac))),
                                                random_state=args.seed)))
            logging.info(f"Subsampled to {len(obs)} cells")
        joinids = obs["soma_joinid"].to_numpy().tolist()

        adata = cellxgene_census.get_anndata(
            census, organism="Mus musculus",
            obs_value_filter=vf, obs_coords=joinids)

    # var_names -> mouse symbols
    if "feature_name" in adata.var.columns:
        adata.var["ensembl_id"] = adata.var_names
        adata.var_names = adata.var.index = adata.var["feature_name"].astype(str)
    adata.var_names_make_unique()
    # Avoid h5ad write error when var index name collides with a var column
    adata.var.index.name = None
    adata.write(args.outdir / "mouse_lung.h5ad")
    logging.info(f"Wrote mouse_lung.h5ad: {adata.shape}")


if __name__ == "__main__":
    main()
