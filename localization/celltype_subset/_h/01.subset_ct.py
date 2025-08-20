import os
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
import harmonypy as hm
from pathlib import Path
from warnings import warn

def subset_data(input_file, COMPARTMENT=False):
    """Subset lung single-cell data."""
    # Load AnnData
    adata = sc.read_h5ad(input_file)

    # Ensure count layer
    if "soupX" in adata.layers:
        adata.layers["counts"] = adata.layers["counts"] if "counts" in adata.layers else adata.layers["soupX"]
    else:
        if "counts" not in adata.layers and adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()

    # Remove missing or unknown annotations (level 4)
    mask = (
        adata.obs["ann_level_4"].notna() &
        ~adata.obs["ann_level_4"].isin(["None", "Unknown"])
    )
    adata = adata[mask].copy()

    # Update annotation columns
    adata.obs["subclusters"] = adata.obs["ann_finest_level"]
    adata.obs["cell_type"] = adata.obs["ann_level_4"]
    adata.obs["clusters"] = adata.obs["cell_type"]
    adata.obs["compartment"] = adata.obs["ann_level_1"]
    adata.obs["patient"] = adata.obs["donor_id"]

    # Remove studies with fewer than 20 pericytes
    if "study" in adata.obs:
        study_counts = adata.obs["study"].value_counts()
        valid_studies = study_counts[study_counts >= 20].index
        adata = adata[adata.obs["study"].isin(valid_studies)].copy()
    else:
        warn("Variable 'study' not found in observation metadata; skipping study filter.")

    # Subset data
    if COMPARTMENT:
        adata = adata[adata.obs["compartment"].isin(["Stroma"])].copy()
    else:
        adata = adata[adata.obs["subclusters"].isin(["Pericytes"])].copy()
    return adata


def preprocess_data(adata):
    """Preprocess and batch-correct data with Harmony."""
    batch_vars = [
        "donor_id", "data", "assay", "tissue_sampling_method", "sequencing_platform",
        "development_stage", "tissue", "subject_type", "study",
        "lung_condition", "sex", "self_reported_ethnicity", "age_or_mean_of_age_range"
    ]
    present_vars = [v for v in batch_vars if v in adata.obs.columns]
    missing_vars = [v for v in batch_vars if v not in adata.obs.columns]
    if missing_vars:
        warn(f"Missing batch variables: {', '.join(missing_vars)}")

    # Normalize total counts and log1p transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Run PCA
    sc.tl.pca(adata, svd_solver="arpack")

    # Run Harmony if batch vars exist
    if present_vars:
        meta = adata.obs[present_vars]
        harmony_embed = hm.run_harmony(adata.obsm["X_pca"], meta, vars_use=present_vars)
        adata.obsm["X_pca_harmony"] = harmony_embed.Z_corr.T
    else:
        warn("No batch variables found for Harmony correction.")
    return adata


def main():
    parser = argparse.ArgumentParser(description="Subset and preprocess lung single-cell data.")
    parser.add_argument("model", type=str, help="Model name (string)")
    args = parser.parse_args()

    model = args.model
    print("Model selected:", model)

    for COMPARTMENT in [False, True]:
        label = "stroma" if COMPARTMENT else "pericyte"
        out_file = f"{label}.hlca_{model}.dataset.h5ad"
        in_file = Path("inputs/hlca/_m") / f"hlca_{model}.h5ad"

        # Run processing
        adata = subset_data(in_file, COMPARTMENT)
        adata = preprocess_data(adata)

        # Save
        adata.write(out_file)
        print(f"Saved: {out_file}")

    # Reproducibility info
    import sessioninfo
    sessioninfo.show()


if __name__ == "__main__":
    main()

    
