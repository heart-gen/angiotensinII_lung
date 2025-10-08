## This is a converted script from our original R version.
## Issues with loading the full model with `zellkonverter`.
import os
import torch, scvi
import numpy as np
import session_info
import scanpy as sc
import harmonypy as hm
from pyhere import here
from pathlib import Path
import matplotlib.pyplot as plt

torch.set_float32_matmul_precision("high")

def subset_data(subset_key: str = "cell_type",
                subset_value: str = "Pericyte"):
    """Subset lung single-cell data for label transfer."""
    # Load AnnData
    input_path = Path('../../_m/ipf_dataset.h5ad')
    adata = sc.read_h5ad(Path(input_path))
    
    # Ensure count layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    mask = adata.obs[subset_key].eq(subset_value)

    return adata[mask].copy()


def check_data(adata, outdir="qc_plots"):
    # Preprocess and dimensionality reduction
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    # Plot & Save
    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)

    # Define the filename base
    umap_file = os.path.join(outdir, "umap_overview")

    # Create and save the plot
    sc.pl.umap(adata, color=["Library_Identity", "patient", "disease"],
               wspace=0.4, ncols=1, save=None, show=False)

    # Save manually for full control
    plt.tight_layout()
    plt.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{umap_file}.pdf", bbox_inches="tight")
    plt.close()
    return adata


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""    
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


def process_query_data():
    # Preprocessing query data
    adata = subset_data()
    adata = check_data(adata, outdir="qc_plots")
    adata = preprocess_data(adata)
    return adata


def load_reference():
    input_path = Path(here("localization/celltype_subset/_m",
                           "pericyte.hlca_full.subclustered.analysis.h5ad"))
    adata = sc.read_h5ad(input_path)
    adata.obs["celltype"] = adata.obs["leiden"]
    mask = adata.var["feature_name"].notna() & (adata.var["feature_name"] != "")
    adata.var_names = np.where(mask, adata.var["feature_name"], adata.var_names)
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    return adata


def prepare_data(query_adata, ref_adata):
    # Ensure common genes
    common_genes = ref_adata.var_names.intersection(query_adata.var_names)
    ref_adata = ref_adata[:, common_genes].copy()
    query_adata = query_adata[:, common_genes].copy()

    # Preprocessing
    hvgs = ref_adata.var.highly_variable
    ref_hvg = ref_adata[:, hvgs].copy()
    query_hvg = query_adata[:, hvgs].copy()
    return ref_hvg, query_hvg


def train_model(ref_hvg):
    # Setup and train model SCANVI model
    scvi.model.SCVI.setup_anndata(ref_hvg, layer="counts")
    vae = scvi.model.SCVI(ref_hvg, n_latent=30, n_layers=2)
    vae.train()

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="unknown",
        labels_key="celltype"
    )
    scanvi_model.train(
        max_epochs=500,
        early_stopping=True,
        early_stopping_patience=10,
        early_stopping_monitor="elbo_validation",
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=5
    )

    # Save model for downstream
    scanvi_model.save("scanvi_model/", overwrite=True)


def transfer_labels(ref_hvg, query_hvg, query_adata, confidence_threshold=0.7):
    # Load model with reference data
    scanvi_model = scvi.model.SCANVI.load("scanvi_model/", adata=ref_hvg)

    # Add latent variables
    ref_hvg.obsm["X_scANVI"] = scanvi_model.get_latent_representation(ref_hvg)
    query_adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(query_hvg)
    query_hvg.obsm["X_scANVI"] = scanvi_model.get_latent_representation(query_hvg)

    # Labels transfer
    query_labels = scanvi_model.predict(query_hvg)
    query_adata.obs["predicted_labels"] = query_labels

    # Get prediction probs
    prediction_probs = scanvi_model.predict(query_hvg, soft=True)
    query_adata.obs["prediction_confidence"] = np.max(prediction_probs, axis=1)
    query_adata.obs["high_confidence"] = \
        query_adata.obs["prediction_confidence"] > confidence_threshold

    # Filter based on confidence
    high_confidence_mask = query_adata.obs["prediction_confidence"] >= confidence_threshold
    filtered_adata = query_adata[high_confidence_mask, :].copy()
    print(f"95th percentile: {np.percentile(query_adata.obs['prediction_confidence'], 95):.2f}")
    print(f"Median confidence: {np.median(query_adata.obs['prediction_confidence']):.2f}")

    # Probability histogram
    plt.hist(query_adata.obs["prediction_confidence"], bins=30)
    plt.xlabel("Prediction confidence")
    plt.ylabel("Cell count")
    plt.title("Prediction confidence distribution")
    plt.savefig("prediction_confidence_hist.png"); plt.close()

    # Reliability curve
    confidences = np.max(prediction_probs, axis=1)
    predictions = np.argmax(prediction_probs, axis=1)
    # Treat predicted label correctness as "true" (no ground truth)
    frac_pos = np.ones_like(predictions)  
    prob_true, prob_pred = calibration_curve(frac_pos, confidences,
                                             n_bins=10, strategy="uniform")
    plt.plot(prob_pred, prob_true, marker="o")
    plt.plot([0, 1], [0, 1], "--", color="gray")
    plt.xlabel("Predicted confidence")
    plt.ylabel("Fraction positive (proxy)")
    plt.title("Reliability curve (approx, no GT)")
    plt.savefig("reliability_curve.png"); plt.close()

    # Cell type composition summary
    original_counts = query_adata.obs["predicted_labels"].value_counts()
    filtered_counts = filtered_adata.obs["predicted_labels"].value_counts()
    df = pd.DataFrame({"Original": original_counts,
                       "Filtered": filtered_counts})\
           .reset_index()
    df.columns = ['Cell_Type', 'Original', 'Filtered']
    df['Difference'] = df['Original'] - df['Filtered']
    df['Percent Retained'] = (df['Filtered'] / df['Original'] * 100).round(2)
    df.sort_values('Difference', ascending=False)\
      .to_csv("annotated_clusters_filtering_summary.tsv", sep="\t",
              index=False)

    return filtered_adata


def cluster_data(adata):
    one_hot = pd.get_dummies(adata.obs["predicted_labels"])
    weighted_one_hot = one_hot.multiply(adata.obs["prediction_confidence"], axis=0)
    weighted_one_hot_norm = (weighted_one_hot - weighted_one_hot.mean()) / weighted_one_hot.std()
    combined_rep = np.column_stack([adata.obsm["X_pca_harmony"],
                                    weighted_one_hot_norm.values])
    adata.obsm["X_combined"] = combined_rep
    sc.pp.neighbors(adata, use_rep='X_combined', n_neighbors=15, n_pcs=30)
    return adata


def visualize_clusters(adata, resolution):
    sc.tl.leiden(adata, resolution=resolution)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='disease', ax=axes[1],show=False,title='Disease')
    plt.tight_layout()
    plt.savefig('harmony_batch_correction.clustering.png')
    plt.savefig('harmony_batch_correction.clustering.pdf')
    plt.close(fig)
    return adata


def additional_viz(adata):
    # Clusters - UMAP
    sc.pl.umap(adata, color=['leiden', 'predicted_labels'],
               show=False, legend_loc="on data")
    plt.savefig('umap.pred_celltypes.png', dpi=300, bbox_inches='tight')
    plt.savefig('umap.pred_celltypes.pdf', bbox_inches='tight')

    # Clusters t-SNE
    sc.tl.tsne(adata)
    sc.pl.tsne(adata, color=['leiden', 'predicted_labels'],
               show=False, legend_loc="on data")
    plt.savefig('tsne.pred_celltypes.png', dpi=300, bbox_inches='tight')
    plt.savefig('tsne.pred_celltypes.pdf', bbox_inches='tight')


def neighborhood_purity(adata, label_key="predicted_labels", k=30):
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=k)
    conn = adata.obsp["connectivities"].tocoo()
    purity_scores = []
    for i in range(adata.n_obs):
        neigh = conn.col[conn.row == i]
        same = np.mean(adata.obs[label_key].iloc[neigh] == adata.obs[label_key].iloc[i])
        purity_scores.append(same)
    mean_purity = np.mean(purity_scores)
    print(f"Neighborhood purity (k={k}): {mean_purity:.2f}")


def mapping_distance(ref_adata, query_adata, k=10):
    ref_latent = ref_adata.obsm["X_scANVI"]
    query_latent = query_adata.obsm["X_scANVI"]
    tree = cKDTree(ref_latent)
    dists, _ = tree.query(query_latent, k=k)
    mean_dist = np.mean(dists)
    print(f"Mean {k}-NN mapping distance (query -> reference): {mean_dist:.3f}")


def composition_parity(ref_adata, query_adata):
    ref_counts = ref_adata.obs["celltype"].value_counts(normalize=True)
    query_counts = query_adata.obs["predicted_labels"].value_counts(normalize=True)
    common = ref_counts.index.intersection(query_counts.index)
    chi2, p = chisquare(query_counts.loc[common], f_exp=ref_counts.loc[common] * query_counts.sum())
    print(f"Composition parity chi-square={chi2:.2f}, p={p:.3e}")


def main():
    # Load data
    query_adata = process_query_data()
    ref_adata = load_reference()
    if "counts" not in ref_adata.layers:
        ref_adata.layers["counts"] = ref_adata.

    # Train model
    ref_hvg, query_hvg = prepare_data(query_adata, ref_adata)
    train_model(ref_hvg)

    # Transfer labels
    adata = transfer_labels(ref_hvg, query_hvg, query_adata, 0.9)
    
    # Save preprocessed objects for GPU step
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    # Reproducibility info
    session_info.show()


if __name__ == "__main__":
    main()

    
