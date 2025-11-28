## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.

import argparse
import session_info
import scanpy as sc
import harmonypy as hm
from pyhere import here
from pathlib import Path

def subset_data(input_file, COMPARTMENT=False):
    """Subset lung single-cell data, harmonizing annotations and filtering studies."""
    # Load AnnData
    adata = sc.read_h5ad(here(input_file))

    # Ensure count layer
    if "counts" not in adata.layers:
        if "soupX" in adata.layers:
            adata.layers["counts"] = adata.layers["soupX"]
        elif adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
        else:
            raise ValueError("No suitable count layer found (expected 'counts', 'soupX', or .raw).")

    # Define categories to add
    required_cats = ["Vascular smooth muscle", "Mesothelium", "Myofibroblasts"]
    for cat in required_cats:
        if cat not in adata.obs["ann_level_4"].cat.categories:
            adata.obs["ann_level_4"] = adata.obs["ann_level_4"].cat.add_categories([cat])

    # Handle annotation issues
    vsm_clusters = {"Smooth muscle", "Smooth muscle FAM83D+",
                    "SM activated stress response"}
    fixes = {
        "Vascular smooth muscle": vsm_clusters,
        "Mesothelium": {"Mesothelium"},
        "Myofibroblasts": {"Myofibroblasts"},
    }

    for label, fine_clusters in fixes.items():
        cond = (
            adata.obs["ann_finest_level"].isin(fine_clusters)
            & (adata.obs["ann_level_4"].isna() | (adata.obs["ann_level_4"] == "None"))
        )
        adata.obs.loc[cond, "ann_level_4"] = label

    # Update annotation columns
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

    # Filter studies (>= 20 cells)
    if "study" in adata.obs:
        study_counts = adata.obs["study"].value_counts()
        valid_studies = study_counts[study_counts >= 20].index
        adata = adata[adata.obs["study"].isin(valid_studies)].copy()
    else:
        print("Warning: Variable 'study' not found in observation metadata; skipping study filter.")

    # Subset data
    if COMPARTMENT:
        mask = adata.obs["compartment"].eq("Stroma")
    else:
        mask = adata.obs["subclusters"].eq(["Pericytes", 'EC general capillary', 'EC aerocyte capillary'])

    adata = adata[mask].copy()

    return adata


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
    parser = argparse.ArgumentParser(description="Subset and preprocess lung single-cell data.")
    parser.add_argument("--model", type=str, default="core",
                        help="Model name (string). Default: 'core'")
    parser.add_argument("--m_iter", type=int, default=30,
                        help="Harmony max iterations. Default: 30.")
    args = parser.parse_args()

    model = args.model
    print("Model selected:", model)

    for COMPARTMENT in False: # Just runt he pericyte+ subclustering
        label = "stroma" if COMPARTMENT else "pericyte"
        out_file = f"{label}.hlca_{model}.dataset.h5ad"
        in_file = Path("inputs/hlca/_m") / f"hlca_{model}.h5ad"

        # Run processing
        adata = subset_data(in_file, COMPARTMENT)
        adata = preprocess_data(adata, max_iter=args.m_iter)

        # Save
        adata.write(out_file)
        print(f"Saved: {out_file}")

    # Reproducibility info
    import session_info
    session_info.show()


if __name__ == "__main__":
    main()

    
