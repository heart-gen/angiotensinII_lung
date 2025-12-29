"""Preprocess reference data for localization workflow."""
import pandas as pd
import scanpy as sc
import session_info
import harmonypy as hm
from pyhere import here
from pathlib import Path
from pandas.api.types import CategoricalDtype

def subset_data():
    input_path = Path(here("inputs/hlca/_m/hlca_full.h5ad"))
    adata = sc.read_h5ad(input_path)
    # Ensure count layer
    if "counts" not in adata.layers:
        if "soupX" in adata.layers:
            adata.layers["counts"] = adata.layers["soupX"]
        elif adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
        else:
            raise ValueError("No suitable count layer found (expected 'counts', 'soupX', or .raw).")

    required_cats = [
        "Vascular smooth muscle", "Mesothelium", "Myofibroblasts",
        "AT1", "AT2", "Hematopoietic stem cells", "Lymphatic EC",
        'Mast cells', "Monocyte-derived Mph"
    ]
    if not isinstance(adata.obs["ann_level_4"].dtype, CategoricalDtype):
        adata.obs["ann_level_4"] = adata.obs["ann_level_4"].astype("category")

    for cat in required_cats:
        if cat not in adata.obs["ann_level_4"].cat.categories:
            adata.obs["ann_level_4"] = adata.obs["ann_level_4"].cat.add_categories([cat])

    vsm_clusters = {"Smooth muscle", "Smooth muscle FAM83D+",
                    "SM activated stress response"}
    lec_clusters = {"Lymphatic EC differentiating",
                    "Lymphatic EC mature",
                    "Lymphatic EC proliferating"}
    fixes = {
        "Vascular smooth muscle": vsm_clusters,
        "Mesothelium": {"Mesothelium"},
        "Myofibroblasts": {"Myofibroblasts"},
        "AT1": {"AT1"}, "AT2": {"AT2"},
        "Hematopoietic stem cells": {"Hematopoietic stem cells"},
        "Lymphatic EC": lec_clusters,
        "Mast cells": {"Mast cells"},
        "Monocyte-derived Mph": {"Monocyte-derived Mph"},
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

    mask = adata.obs["compartment"].eq("Stroma")
    return adata, adata[mask].copy()


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""
    batch_vars = [
        "donor_id", "data", "assay", "tissue_sampling_method", "sequencing_platform",
        "development_stage", "tissue", "subject_type", "study",
        "lung_condition", "sex", "self_reported_ethnicity", "age_or_mean_of_age_range"
    ]

    present_vars = [v for v in batch_vars if v in adata.obs.columns]
    missing_vars = [v for v in batch_vars if v not in adata.obs.columns]

    if missing_vars:
        print(f"Warning: Missing batch variables: {', '.join(missing_vars)}")

    # Preprocess data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, batch_key='study'
    )
    sc.tl.pca(
        adata, n_comps=50, mask_var="highly_variable",
        svd_solver="arpack"
    )

    # Run harmony
    if present_vars:
        meta = adata.obs[present_vars].copy()
        
        for col in present_vars:
            meta[col] = meta[col].astype("category").astype(str)

        harmony_embed = hm.run_harmony(
            adata.obsm["X_pca"], meta, vars_use=present_vars,
            max_iter_harmony=max_iter, epsilon_harmony=1e-5, 
            random_state=seed
        )
        adata.obsm["X_pca_harmony"] = harmony_embed.Z_corr.T
    else:
        print("Warning: No batch variables found for Harmony correction.")

    return adata


def main():
    fadata, adata = subset_data()

    # Process all data
    sc.pp.normalize_total(fadata, target_sum=1e4)
    sc.pp.log1p(fadata)
    sc.pp.highly_variable_genes(
        fadata, n_top_genes=2000, batch_key='study'
    )
    fadata.write("hlca_full.dataset.h5ad")
    del fadata

    # Process stromal only
    adata = preprocess_data(adata, max_iter=75)
    adata.write("stroma.hlca_full.dataset.h5ad")
    session_info.show()


if __name__ == "__main__":
    main()
