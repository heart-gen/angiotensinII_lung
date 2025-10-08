"""
Single-cell label transfer.
"""
import scvi, os
import numpy as np
import session_info
import pandas as pd
import scanpy as sc
import seaborn as sns
import harmonypy as hm
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import chisquare
from sklearn.calibration import calibration_curve

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


def preprocess_data(adata, max_iter: int = 30, seed: int = 13):
    """Preprocess and batch-correct data with Harmony."""    
    # Preprocess and dimensionality reduction
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    sc.pp.scale(adata, zero_center=False)
    sc.tl.pca(adata)
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


def process_query_data():
    # Preprocessing query data
    adata = subset_data()
    adata = preprocess_data(adata)
    return adata


def transfer_labels(confidence_threshold=0.7, outdir="qc_plots"):
    # Load data
    ref_hvg = sc.read_h5ad("ref_hvg.h5ad")
    query_adata = process_query_data()
    query_adata.X = query_adata.layers["counts"]
    query_hvg = sc.read_h5ad("query_hvg.h5ad")

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
    hist_file = os.path.join(outdir, "prediction_confidence_hist")
    plt.hist(query_adata.obs["prediction_confidence"], bins=30)
    plt.xlabel("Prediction confidence")
    plt.ylabel("Cell count")
    plt.title("Prediction confidence distribution")
    plt.savefig(f"{hist_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{hist_file}.pdf", bbox_inches="tight")
    plt.close()

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

    return filtered_adata, query_hvg, ref_hvg


def cluster_data(adata):
    one_hot = pd.get_dummies(adata.obs["predicted_labels"])
    weighted_one_hot = one_hot.multiply(adata.obs["prediction_confidence"], axis=0)
    weighted_one_hot_norm = (weighted_one_hot - weighted_one_hot.mean()) / weighted_one_hot.std()
    combined_rep = np.column_stack([adata.obsm["X_pca_harmony"],
                                    weighted_one_hot_norm.values])
    adata.obsm["X_combined"] = combined_rep
    sc.pp.neighbors(adata, use_rep='X_combined', n_neighbors=15, n_pcs=30)
    return adata


def visualize_clusters(adata, resolution, outdir):
    sc.tl.leiden(adata, resolution=resolution)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='disease', ax=axes[1],show=False, title='Disease')
    sc.pl.umap(adata, color='patient', ax=axes[0], show=False,
               title='Donors')
    plt.tight_layout()
    plt.savefig(f'{outdir}/harmony_batch_correction.clustering.png',
                dpi=300, bbox_inches="tight")
    plt.savefig(f'{outdir}/harmony_batch_correction.clustering.pdf',
                bbox_inches="tight")
    plt.close(fig)
    return adata


def additional_viz(adata, outdir):
    # Clusters - UMAP
    sc.pl.umap(adata, color=['leiden', 'predicted_labels'],
               show=False, legend_loc="on data")
    plt.savefig(f'{outdir}/umap.pred_celltypes.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outdir}/umap.pred_celltypes.pdf', bbox_inches='tight')

    # Clusters t-SNE
    sc.tl.tsne(adata)
    sc.pl.tsne(adata, color=['leiden', 'predicted_labels'],
               show=False, legend_loc="on data")
    plt.savefig(f'{outdir}/tsne.pred_celltypes.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{outdir}/tsne.pred_celltypes.pdf', bbox_inches='tight')
    plt.close()


def neighborhood_purity(adata, label_key="predicted_labels", k=30, outdir="qc_plots"):
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=k)
    conn = adata.obsp["connectivities"].tocoo()

    purity_scores = np.zeros(adata.n_obs)
    for i in range(adata.n_obs):
        neigh = conn.col[conn.row == i]
        same = np.mean(adata.obs[label_key].iloc[neigh] == adata.obs[label_key].iloc[i])
        purity_scores.append(same)

    adata.obs[f"purity_k{k}"] = purity_scores
    mean_purity = purity_scores.mean()
    print(f"Neighborhood purity (k={k}): {mean_purity:.2f}")

    # Visualization
    purity_file = os.path.join(outdir, f"neighborhood_purity_k{k}")
    plt.figure(figsize=(6, 4))
    sns.histplot(purity_scores, bins=30, kde=True, color="steelblue")
    plt.axvline(mean_purity, color="red", linestyle="--", label=f"Mean = {mean_purity:.2f}")
    plt.xlabel("Neighborhood Purity")
    plt.ylabel("Cell Count")
    plt.title(f"Neighborhood Purity Distribution (k={k})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{purity_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{purity_file}.pdf", bbox_inches="tight")
    plt.close()


def mapping_distance(ref_adata, query_adata, k=10, outdir="qc_plots"):
    ref_latent = ref_adata.obsm["X_scANVI"]
    query_latent = query_adata.obsm["X_scANVI"]
    tree = cKDTree(ref_latent)
    dists, _ = tree.query(query_latent, k=k)
    mean_dist = dists.mean()
    query_adata.obs[f"mapping_dist_k{k}"] = dists.mean(axis=1)

    print(f"Mean {k}-NN mapping distance (query -> reference): {mean_dist:.3f}")

    # Visualization
    dists_file = os.path.join(outdir, f"mapping_distance_k{k}")
    plt.figure(figsize=(6, 4))
    sns.histplot(query_adata.obs[f"mapping_dist_k{k}"], bins=30, kde=True, color="seagreen")
    plt.axvline(mean_dist, color="red", linestyle="--", label=f"Mean = {mean_dist:.3f}")
    plt.xlabel("Mean k-NN Distance to Reference")
    plt.ylabel("Query Cells")
    plt.title(f"Mapping Distance Distribution (k={k})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{dists_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{dists_file}.pdf", bbox_inches="tight")
    plt.close()

    if "X_umap" in query_adata.obsm:
        umap_file = os.path.join(outdir, f"umap_diffuse_{k}")
        sc.pl.umap(query_adata, color=f"mapping_dist_k{k}", cmap="viridis",
                   title=f"Mapping Distance (k={k})", save=None, show=False)
        # Save manually for full control
        plt.tight_layout()
        plt.savefig(f"{umap_file}.png", dpi=300, bbox_inches="tight")
        plt.savefig(f"{umap_file}.pdf", bbox_inches="tight")
        plt.close()


def composition_parity(ref_adata, query_adata, outdir="qc_plots"):
    ref_counts = ref_adata.obs["celltype"].value_counts(normalize=True)
    query_counts = query_adata.obs["predicted_labels"].value_counts(normalize=True)
    common = ref_counts.index.intersection(query_counts.index)

    chi2, p = chisquare(query_counts.loc[common],
                        f_exp=ref_counts.loc[common] * query_counts.sum())
    print(f"Composition parity chi-square={chi2:.2f}, p={p:.3e}")

    # Visualization
    comp_file = os.path.join(outdir, "composition_parity")
    comp_df = (
        pd.DataFrame({
            "Reference": ref_counts.loc[common],
            "Query": query_counts.loc[common]
        })
        .reset_index()
        .melt(id_vars="index", var_name="Dataset", value_name="Proportion")
        .rename(columns={"index": "CellType"})
    )

    plt.figure(figsize=(8, 5))
    sns.barplot(data=comp_df, x="CellType", y="Proportion",
                hue="Dataset", palette="Set2")
    plt.xticks(rotation=45, ha="right")
    plt.title()
    plt.tight_layout()
    plt.savefig(f"{comp_file}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{comp_file}.pdf", bbox_inches="tight")
    plt.close()


def main():
    # Make directory
    outdir = "qc_plots"
    os.makedirs(outdir, exist_ok=True)

    # Transfer labels    
    adata, query_hvg, ref_hvg = transfer_labels(0.9, outdir=outdir)

    # Cluster data and visualization
    adata = cluster_data(adata)
    adata = visualize_clusters(adata, 0.45, outdir)
    additional_viz(adata, outdir)

    # Quality control analysis
    neighborhood_purity(adata, outdir=outdir)
    mapping_distance(ref_hvg, query_hvg, outdir=outdir)
    composition_parity(ref_hvg, adata, outdir=outdir)
    
    # Save data
    adata.write_h5ad("clustered_data.h5ad", compression="gzip")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
