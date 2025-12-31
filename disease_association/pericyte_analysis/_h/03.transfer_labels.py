"""Single-cell label transfer and QC."""
import scvi
import logging
import numpy as np
import session_info
import pandas as pd
import scanpy as sc
import harmonypy as hm
from pathlib import Path
from anndata import AnnData
import matplotlib.pyplot as plt

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def set_seed(SEED):
    np.random.seed(13)
    sc.settings.seed = 13


def load_anndata(path: Path, label: str):
    logging.info("Loading %s AnnData from %s", label, path)
    adata = sc.read_h5ad(path)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"{label} AnnData is empty: {path}")
    return adata


def subset_data(input_path: Path, subset_key: str = "cell_type",
                subset_value: str = "Pericyte"):
    logging.info("Subsetting IPF data for %s", subset_value)
    adata = sc.read_h5ad(input_path)
    
    # Ensure count layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if subset_key not in adata.obs:
        raise KeyError(f"{subset_key} not found in adata.obs")
    mask = adata.obs[subset_key].eq(subset_value)

    return adata[mask].copy()


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    logging.info("Preprocessing the IPF query data")
    # Preprocess and dimensionality reduction
    sc.pp.normalize_total(adata, layer="counts")
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Run harmony
    batch_vars = ["patient"]
    meta = adata.obs[batch_vars].copy()
        
    for col in batch_vars:
        meta[col] = meta[col].astype("category").astype(str)

    harmony_embed = hm.run_harmony(
        adata.obsm["X_pca"], meta, vars_use=batch_vars,
        max_iter_harmony=max_iter, epsilon_harmony=1e-5, 
        random_state=seed
    )
    adata.obsm["X_pca_harmony"] = harmony_embed.Z_corr.T

    return adata


def process_query_data(path: Path):
    adata = subset_data(path)
    adata = preprocess_data(adata)
    return adata


def align_query_objects(query_hvg: AnnData, query_full: AnnData):
    missing = query_hvg.obs_names.difference(query_full.obs_names)
    if len(missing) > 0:
        raise ValueError("Query HVG object contains cells absent from full query data")
    aligned = query_full[query_hvg.obs_names].copy()
    if not np.array_equal(aligned.obs_names.values, query_hvg.obs_names.values):
        raise AssertionError("Failed to align query HVG and full objects")
    return aligned


def load_model(model_dir: Path, ref_hvg: AnnData) -> scvi.model.SCANVI:
    logging.info("Loading SCANVI reference model from %s", model_dir)
    ref_hvg.obs["batch"] = ref_hvg.obs["dataset"]
    model = scvi.model.SCANVI.load(model_dir, adata=ref_hvg)
    return model


def add_latent_representations(
    scanvi_ref: scvi.model.SCANVI,  ref_hvg: AnnData, query_hvg: AnnData,
    query_full: AnnData,
):
    logging.info("Computing latent representations")
    ref_latent   = scanvi_ref.get_latent_representation(ref_hvg)
    query_latent = scanvi_ref.get_latent_representation(query_hvg)
    ref_hvg.obsm["X_scANVI"]    = ref_latent
    query_hvg.obsm["X_scANVI"]  = query_latent
    query_full.obsm["X_scANVI"] = query_latent.copy()


def predict_labels(scanvi_ref: scvi.model.SCANVI, query_hvg: AnnData, query_full: AnnData):
    logging.info("Predicting soft labels for query data")
    prob_df = scanvi_ref.predict(query_hvg, soft=True)
    if not isinstance(prob_df, pd.DataFrame):
        raise TypeError("Expected SCANVI.predict(soft=True) to return a DataFrame")

    pred = prob_df.idxmax(axis=1).astype(str)
    logging.info("Predicting soft probs and labels for query data")
    query_full.obsm["pred_probs"] = prob_df
    query_full.obs["predicted_labels"] = pred.astype("category")
    query_hvg.obs["predicted_labels"]  = pred.astype("category")
    return prob_df


def compute_uncertainty_metrics(prob_df: pd.DataFrame):
    probs = prob_df.to_numpy()
    eps   = np.finfo(np.float64).eps
    clipped      = np.clip(probs, eps, 1.0)
    confidences  = clipped.max(axis=1)
    sorted_probs = np.sort(clipped, axis=1)
    margins   = sorted_probs[:, -1] - sorted_probs[:, -2] if clipped.shape[1] > 1 else np.zeros(clipped.shape[0])
    entropies = -np.sum(clipped * np.log(clipped), axis=1)
    # Probability histogram
    plt.hist(confidences, bins=30)
    plt.xlabel("Prediction confidence")
    plt.ylabel("Cell count")
    plt.title("Prediction confidence distribution")
    plt.savefig("prediction_confidence_hist.png"); plt.close()
    return pd.DataFrame({"prediction_confidence": confidences,
                         "prediction_entropy": entropies,
                         "margin": margins}, index=prob_df.index)


def generate_filtering_summary(adata: AnnData, threshold: float, outdir: Path):
    if "predicted_labels" not in adata.obs or "prediction_confidence" not in adata.obs:
        logging.warning("Missing columns for filtering summary")
        return

    df = adata.obs[["predicted_labels", "prediction_confidence"]].copy()
    df["retained"] = df["prediction_confidence"] >= threshold
    summary = (
        df.groupby("predicted_labels")
        .agg(total_cells=("predicted_labels", "size"), retained_cells=("retained", "sum"))
        .assign(percent_retained=lambda x: 100 * x["retained_cells"] / x["total_cells"].clip(lower=1))
    )
    path = outdir / "annotation_filtering_summary.tsv"
    logging.info("Writing table to %s", path)
    summary.to_csv(path, sep="\t")


def build_graph_and_cluster(adata: AnnData, neighbors: int, resolution: float, seed: int):
    logging.info("Constructing neighborhood graph (n_neighbors=%d)", neighbors)
    sc.pp.neighbors(adata, n_neighbors=neighbors, use_rep="X_scANVI", random_state=seed)
    logging.info("Computing UMAP embedding")
    sc.tl.umap(adata, random_state=seed)
    logging.info("Running Leiden clustering (resolution=%.3f)", resolution)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                 random_state=seed)


def visualize_clusters(adata: AnnData, outdir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='disease', ax=axes[1],show=False, title='Disease')
    sc.pl.umap(adata, color='patient', ax=axes[0], show=False, title='Donors')
    plt.tight_layout()
    plt.savefig(outdir / 'harmony_batch_correction.clustering.png',
                dpi=300, bbox_inches="tight")
    plt.savefig(outdir / 'harmony_batch_correction.clustering.pdf',
                bbox_inches="tight")
    plt.close(fig)


def additional_viz(adata: AnnData, outdir: Path):
    # Clusters - UMAP
    sc.pl.umap(adata, color=['leiden', 'predicted_labels'],
               show=False, legend_loc="on data")
    plt.savefig(outdir / 'umap.pred_celltypes.png', dpi=300, bbox_inches='tight')
    plt.savefig(outdir / 'umap.pred_celltypes.pdf', bbox_inches='tight')
    plt.close()


def _prepare_pred_probs_for_h5(adata: AnnData, key: str = "pred_probs") -> None:
    """Make obsm[key] HDF5-safe by storing as matrix and saving columns in .uns."""
    pp = adata.obsm.get(key)
    if pp is None:
        return
    if isinstance(pp, pd.DataFrame):
        # keep original column names safely in .uns
        adata.uns[f"{key}_columns"] = [str(c) for c in pp.columns]
        # store data as a dense float32 matrix for compactness/compatibility
        adata.obsm[key] = pp.to_numpy(dtype=np.float32, copy=False)


def save_adata(adata: AnnData, path: Path) -> None:
    logging.info("Saving AnnData to %s", path)
    _prepare_pred_probs_for_h5(adata, "pred_probs")
    adata.write(path, compression="gzip")


def main():
    # Configure environment
    configure_logging()
    set_seed(13)

    # Make directory
    outdir = Path("results")
    outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    input_path = Path(here('disease_association/ipf_analysis',
                           '_m/ipf_dataset.h5ad'))
    ref_hvg = load_anndata(Path("ref_hvg.h5ad"), "reference HVG")
    query_hvg = load_anndata(Path("query_hvg.h5ad"), "query HVG")
    query_full = process_query_data(input_path)
    query_full.X = query_full.layers["counts"]
    query_full = align_query_objects(query_hvg, query_full)

    # Transfer model
    scanvi_ref = load_model(Path("scanvi_model/"), ref_hvg)
    add_latent_representations(scanvi_ref, ref_hvg, query_hvg, query_full)
    prob_df = predict_labels(scanvi_ref, query_hvg, query_full)
    uncertainty = compute_uncertainty_metrics(prob_df)
    for column in uncertainty.columns:
        query_full.obs[column] = uncertainty[column]
        query_hvg.obs[column]  = uncertainty[column]

    # Cluster data and visualization
    build_graph_and_cluster(query_full, 50, 0.45, 13)
    visualize_clusters(query_full, outdir)
    additional_viz(query_full, outdir)

    # Filter data summary
    generate_filtering_summary(query_full, 0.75, outdir)

    # Save data
    clustered_path = outdir / "clustered_data.h5ad"
    save_adata(query_full, clustered_path)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
