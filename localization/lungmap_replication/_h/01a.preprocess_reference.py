"""Preprocess reference data for localization workflow."""

import pandas as pd
import scanpy as sc
import session_info
from pandas.api.types import CategoricalDtype
from pathlib import Path
from pyhere import here


def _sanitize_var_for_h5ad(adata):
    if isinstance(adata.var.index.dtype, CategoricalDtype):
        adata.var.index = pd.Index(adata.var.index.astype(str))
    else:
        adata.var.index = pd.Index(adata.var.index.astype(str))

    adata.var.index.name = None
    adata.var_names_make_unique()
    return adata


def load_reference():
    input_path = Path(here("inputs/hlca/_m/hlca_core.h5ad"))
    adata = sc.read_h5ad(input_path)
    if "counts" not in adata.layers:
        if "soupX" in adata.layers:
            adata.layers["counts"] = adata.layers["soupX"]
        elif adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
        else:
            raise ValueError("No suitable count layer found (expected 'counts', 'soupX', or .raw).")

    required_cats = ["Vascular smooth muscle", "Mesothelium", "Myofibroblasts"]
    if not isinstance(adata.obs["ann_level_4"].dtype, CategoricalDtype):
        adata.obs["ann_level_4"] = adata.obs["ann_level_4"].astype("category")

    for cat in required_cats:
        if cat not in adata.obs["ann_level_4"].cat.categories:
            adata.obs["ann_level_4"] = adata.obs["ann_level_4"].cat.add_categories([cat])

    vsm_clusters = {"Smooth muscle", "Smooth muscle FAM83D+", "SM activated stress response"}
    fixes = {
        "Vascular smooth muscle": vsm_clusters,
        "Mesothelium": {"Mesothelium"},
        "Myofibroblasts": {"Myofibroblasts"},
    }

    for label, fine_clusters in fixes.items():
        cond = (
            adata.obs["ann_finest_level"].isin(fine_clusters)
            & (adata.obs["ann_level_4"].isna() | (adata.obs["ann_level_4"].isin(["None", ""])))
        )
        adata.obs.loc[cond, "ann_level_4"] = label

    obs_map = {
        "subclusters": "ann_finest_level",
        "cell_type": "ann_level_4",
        "clusters": "ann_level_4",
        "compartment": "ann_level_1",
        "patient": "donor_id",
    }
    for new, old in obs_map.items():
        if old in adata.obs:
            adata.obs[new] = adata.obs[old]
        else:
            print(f"Warning: '{old}' not found in obs; skipping {new} mapping.")

    if "study" in adata.obs:
        study_counts = adata.obs["study"].value_counts()
        valid_studies = study_counts[study_counts >= 20].index
        adata = adata[adata.obs["study"].isin(valid_studies)].copy()
    else:
        print("Warning: Variable 'study' not found in observation metadata; skipping study filter.")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="study")
    return adata


def main():
    ref_adata = load_reference()
    ref_adata = _sanitize_var_for_h5ad(ref_adata)
    ref_adata.write_h5ad("ref_preprocessed.h5ad", compression="gzip")
    session_info.show()


if __name__ == "__main__":
    main()
