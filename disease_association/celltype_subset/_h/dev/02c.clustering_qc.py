"""
Single-Cell Data Clustering and Visualization

This script clusters the data and extracts marker genes for cell type annotation.

Author: Kynon J Benjamin
Date: August 23, 2024
Revision: September 27, 2025
"""
import scvi
import numpy as np
import session_info
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import chisquare
from sklearn.calibration import calibration_curve

def load_query_data():
    adata = sc.read_h5ad("../../_m/integrated_data.h5ad")
    adata.obs["celltype"] = "unknown"
    return adata


def transfer_labels(confidence_threshold=0.7):
    # Load data
    ref_hvg = sc.read_h5ad("ref_hvg.h5ad")
    query_adata = load_query_data()
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


def visualize_clusters(adata, resolution):
    sc.tl.leiden(adata, resolution=resolution)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='disease', ax=axes[1],show=False,title='Disease')
    sc.pl.umap(adata, color='genome', ax=axes[0], show=False,
               title='Genomes')
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


def identify_marker_genes(adata, label):
    marker_genes = {
        'Radial glia': ['PAX6', 'SOX2', 'VIM'],
        'Intermediate progenitors': ['EOMES', 'NEUROG1'],
        'Excitatory neurons': ['NEUROD6', 'SLC17A7', 'SATB2'],
        'Inhibitory neurons': ['DLX1', 'GAD1', 'GAD2'],
        'Astrocytes': ['GFAP', 'AQP4', 'S100B']
    }

    # Report and filter missing genes
    all_markers = set(sum(marker_genes.values(), []))
    missing = all_markers - set(adata.var_names)
    if missing:
        print(f"[identify_marker_genes] Skipping {len(missing)} markers not found in dataset: {', '.join(sorted(missing))}")

    marker_genes_filtered = {
        k: [g for g in v if g in adata.var_names]
        for k, v in marker_genes.items()
    }

    # Remove empty groups
    marker_genes_filtered = {
        k: v for k, v in marker_genes_filtered.items() if v
    }

    # Differential expression analysis
    sc.tl.rank_genes_groups(adata, label, method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25,
                            save=f"rank_genes_groups-{label}.png")

    # Canonical marker dotplot
    if marker_genes_filtered:
        sc.pl.dotplot(adata, marker_genes_filtered, groupby=label,
                      standard_scale="var", swap_axes=True,
                      save=f"marker_genes_dotplot-{label}.png")
    else:
        print(f"[identify_marker_genes] No canonical markers found for {label}, skipping dotplot.")

    # Save ranked DEGs table
    result_df = sc.get.rank_genes_groups_df(adata, group=None)
    result_df = result_df.sort_values(['group', 'scores'],
                                      ascending=[True, False])
    result_df.to_csv(f"marker_genes_per_cluster-{label}.tsv",
                     sep='\t', index=False)


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
    # Transfer labels
    adata, query_hvg, ref_hvg = transfer_labels(0.9)

    # Cluster data and visualization
    adata = cluster_data(adata)
    adata = visualize_clusters(adata, 0.45)
    additional_viz(adata)
    identify_marker_genes(adata, "leiden")
    identify_marker_genes(adata, "predicted_labels")

    # Quality control analysis
    neighborhood_purity(adata)
    mapping_distance(ref_hvg, query_hvg)
    composition_parity(ref_hvg, adata)
    
    # Save data
    adata.write_h5ad("clustered_data.h5ad", compression="gzip")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
