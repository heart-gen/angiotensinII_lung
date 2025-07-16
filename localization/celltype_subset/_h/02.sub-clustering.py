## This script subclusters pericytes

import phate
import numpy as np
import scanpy as sc
import session_info
import matplotlib.pyplot as plt

def preprocess_adata(adata):
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)
    return adata


def add_phate_embedding(adata, knn=5, decay=40):
    phate_op = phate.PHATE(knn=knn, decay=decay, n_jobs=-1)
    phate_emb = phate_op.fit_transform(adata.obsm['X_pca'])
    adata.obsm['X_phate'] = phate_emb
    return adata


def perform_clustering(adata, resolution=1.0):
    sc.tl.leiden(adata, resolution=resolution)
    return adata


def analyze_marker_genes(adata, groupby='leiden', method='wilcoxon'):
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
    return adata


def normalize_marker_expression(adata, markers, ref_adata=None,
                                ref_label='Alveolar fibroblasts'):
    ensembl_markers = adata.var[adata.var.feature_name.isin(markers)].index.values
    pericyte_expr = adata[:, ensembl_markers].X.mean(axis=1)
    
    if ref_adata is not None and 'cell_type' in ref_adata.obs:
        fibroblast_cells = ref_adata.obs['cell_type'] == ref_label
        fibroblast_count = fibroblast_cells.sum()
    else:
        fibroblast_cells = adata.obs['cell_type'] == ref_label
        fibroblast_count = fibroblast_cells.sum()

    if fibroblast_count > 0:
        normalized_expr = pericyte_expr / fibroblast_count
    else:
        print("Warning: No fibroblast cells found; skipping normalization")
        normalized_expr = pericyte_expr

    adata.obs['normalized_pericyte_marker_expr'] = normalized_expr
    return adata


def normalize_by_fibroblast_count(adata, ref_adata, ref_label='Alveolar fibroblasts'):
    if ref_adata is None or 'cell_type' not in ref_adata.obs:
        raise ValueError("Reference AnnData (`ref_adata`) with 'cell_type' annotation required.")
        
    fibroblast_count = (ref_adata.obs['cell_type'] == ref_label).sum()
    if fibroblast_count == 0:
        print(f"Warning: No fibroblast cells labeled '{ref_label}' found. No normalization performed.")
        return adata  # or adata unchanged

    normalized_total = adata.X.sum(axis=1) / fibroblast_count
    adata.obs['normalized_total_expression_by_fibroblast'] = normalized_total
    return adata


def plot_clusters_and_markers(adata, marker_genes, cluster_key='leiden'):
    # Plot Leiden Clusters in a separate figure
    plt.figure(figsize=(7, 6))
    sc.pl.umap(adata, color=cluster_key, title='Leiden Clusters', show=True)
    plt.savefig('leiden_clusters.png')
    plt.savefig('leiden_clusters.pdf')
    plt.close()
    
    plt.figure(figsize=(7, 6))
    sc.pl.embedding(adata, basis='X_phate', color=cluster_key, title='PHATE Clusters', show=True)
    plt.savefig('phate_clusters.png')
    plt.savefig('phate_clusters.pdf')
    plt.close()

    # Plot gene expression
    for gene in marker_genes:
        ensembl_id = adata.var[adata.var.feature_name == gene].index.values[0]
        if ensembl_id in adata.var_names:
            plt.figure(figsize=(7, 6))
            sc.pl.umap(adata, color=ensembl_id, title=f'UMAP: {gene} expression', show=True)
            plt.savefig(gene.lower() + '_umap.png')
            plt.savefig(gene.lower() + '_umap.pdf')
            plt.close()

            plt.figure(figsize=(7, 6))
            sc.pl.embedding(adata, basis='X_phate', color=ensembl_id, title=f'PHATE: {gene} expression', show=True)
            plt.savefig(gene.lower() + '_phate.png')
            plt.savefig(gene.lower() + '_phate.pdf')
            plt.close()
        else:
            print(f"Warning: Gene {gene} not found in dataset")

    if 'normalized_pericyte_marker_expr' in adata.obs:
        plt.figure(figsize=(7, 6))
        sc.pl.umap(adata, color='normalized_pericyte_marker_expr',
                   title='Normalized Pericyte Marker Expression', show=True)
        plt.savefig('pericyte_marker_expression.png')
        plt.savefig('pericyte_marker_expression.pdf')
        plt.close()
        
        plt.figure(figsize=(7, 6))
        sc.pl.embedding(adata, basis='X_phate', color='normalized_pericyte_marker_expr',
                        title='Normalized Pericyte Marker Expression', show=True)
        plt.savefig('pericyte_marker_expression_phate.png')
        plt.savefig('pericyte_marker_expression_phate.pdf')
        plt.close()

    if 'normalized_total_expression_by_fibroblast' in adata.obs:
        plt.figure(figsize=(7, 6))
        sc.pl.umap(adata, color='normalized_total_expression_by_fibroblast',
                   title='Normalized Pericyte Expression', show=True)
        plt.savefig('pericyte_fibroblast_expression.png')
        plt.savefig('pericyte_fibroblast_expression.pdf')
        plt.close()

        plt.figure(figsize=(7, 6))
        sc.pl.embedding(adata, basis='X_phate',
                        color='normalized_total_expression_by_fibroblast',
                        title='Normalized Pericyte Expression', show=True)
        plt.savefig('pericyte_fibroblast_expression_phate.png')
        plt.savefig('pericyte_fibroblast_expression_phate.pdf')
        plt.close()


def subcluster_pericytes(
        adata, ref_adata=None, fibroblast_label='Alveolar fibroblasts',
        pericyte_markers=['HIGD1B', 'PDGFRB', 'CSPG4'],
        marker_genes=['AGTR1', 'ACTA2'], phate_knn=5, phate_decay=40,
        leiden_resolution=0.5):
    adata = preprocess_adata(adata)
    adata = add_phate_embedding(adata, knn=phate_knn, decay=phate_decay)
    adata = perform_clustering(adata, resolution=leiden_resolution)
    adata = analyze_marker_genes(adata)
    adata = normalize_marker_expression(adata, pericyte_markers, ref_adata, fibroblast_label)
    adata = normalize_by_fibroblast_count(adata, ref_adata, fibroblast_label)
    plot_clusters_and_markers(adata, marker_genes)
    return adata


def main():
    # Load data
    adata = sc.read_h5ad('../_m/pericyte.hlca_core.dataset.h5ad')
    ref_adata = sc.read_h5ad("../_m/stroma.hlca_core.dataset.h5ad")
    # Subcluster
    adata = subcluster_pericytes(adata, ref_adata, leiden_resolution=0.25)
    # Save the subclusters
    adata.write('pericyte.hlca_core.subclustered.h5ad')
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
