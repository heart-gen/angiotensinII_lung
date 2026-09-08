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
    # PIN THE RELEASE (P2-13, 2026-09-08). This was "stable", a moving alias, and
    # no download log was ever written -- so the release that produced the shipped
    # h5ad is not recorded anywhere and a rerun could silently mix a Census
    # version change into any result. The alias currently resolves to 2025-11-08
    # and has not moved since, so pinning it here reproduces the existing data
    # rather than replacing it. Change this deliberately, never by drift.
    p.add_argument("--census-version", default="2025-11-08")
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

        # Cap by stratified subsample per cell_type (observed groups only).
        #
        # THE MURAL EXEMPTION (P2-13, added 2026-09-08). `max(1, ...)` guards a
        # cell type against reaching zero. It does NOT guard the rarest and most
        # analytically central population against losing 59% of itself, which is
        # what happened: at frac 0.4168 the 99 Census pericytes became 41, and
        # this module's docstring one screen above names pericyte scarcity as its
        # central problem. Every "n = 41, so this can only be a compartment-level
        # claim" sentence downstream rested on a sampling decision, not on the
        # data available.
        #
        # Proportional sampling is unbiased, so this was never a validity error.
        # It is a POWER error, and power is what decides what this module may
        # claim. The three VASCULAR mural types are therefore exempt and taken
        # whole. Cost: 1,601 extra cells, 2.7% of the 60,000 budget, in exchange
        # for 99 pericytes from 26 donors (2.4x the cells, 1.4x the donors) and
        # the full 209 PA-SMC.
        #
        # Note which types are NOT in the exemption: aortic (4,631, one donor)
        # and bronchial (4,456) smooth muscle are airway/systemic, not the lung
        # vascular mural compartment this module compares against pericytes.
        # Exempting those too would cost 6,901 cells for no analytical gain.
        #
        # The non-mural cells stay subsampled: scVI needs them for a well-behaved
        # latent space, and none of them is ever the unit of a claim here.
        MURAL_EXEMPT = ("pericyte",
                        "smooth muscle cell of the pulmonary artery",
                        "vascular associated smooth muscle cell")
        if len(obs) > args.max_cells:
            obs = obs.copy()
            obs["cell_type"] = obs["cell_type"].astype(str)
            exempt = obs[obs.cell_type.isin(MURAL_EXEMPT)]
            rest = obs[~obs.cell_type.isin(MURAL_EXEMPT)]
            # Budget is spent on the non-exempt cells so the cap still means
            # something; the exemption is additive and its size is logged.
            frac = max(0.0, (args.max_cells - len(exempt))) / max(1, len(rest))
            rest = (rest.groupby("cell_type", group_keys=False, observed=True)
                        .apply(lambda d: d.sample(max(1, int(round(len(d) * frac))),
                                                  random_state=args.seed)))
            obs = pd.concat([exempt, rest], ignore_index=True)
            logging.info(
                f"Subsampled to {len(obs)} cells "
                f"(mural exempt and taken whole: {len(exempt)} cells across "
                f"{exempt.cell_type.nunique()} types; frac applied to the rest: "
                f"{frac:.6f})")
            for ctname in MURAL_EXEMPT:
                n = int((obs.cell_type == ctname).sum())
                nd = obs.loc[obs.cell_type == ctname, "donor_id"].nunique()
                logging.info(f"  exempt {ctname}: {n} cells / {nd} donors")
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
