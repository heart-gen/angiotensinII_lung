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


def normalize_marker_expression(adata, markers, ref_adata=None, ref_label='Fibroblasts'):
    pericyte_expr = adata[:, markers].X.mean(axis=1)
    
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


def plot_clusters_and_markers(adata, marker_genes, cluster_key='leiden'):
    sc.pl.umap(adata, color=cluster_key, title='Leiden Clusters')
    sc.pl.embedding(adata, basis='X_phate', color=cluster_key, title='PHATE Clusters')

    for gene in marker_genes:
        if gene in adata.var_names:
            sc.pl.umap(adata, color=gene, title=f'UMAP: {gene} expression')
            sc.pl.embedding(adata, basis='X_phate', color=gene, title=f'PHATE: {gene} expression')
        else:
            print(f"Warning: Gene {gene} not found in dataset")

    if 'normalized_pericyte_marker_expr' in adata.obs:
        sc.pl.umap(adata, color='normalized_pericyte_marker_expr',
                   title='Normalized Pericyte Marker Expression')


def subcluster_pericytes(
        adata, ref_adata=None, fibroblast_label='Fibroblasts',
        pericyte_markers=['HIGD1B', 'PDGFRB', 'CSPG4'],
        marker_genes=['AGTR1', 'ACTA2'], phate_knn=5, phate_decay=40,
        leiden_resolution=1.0):
    preprocess_adata(adata)
    add_phate_embedding(adata, knn=phate_knn, decay=phate_decay)
    perform_clustering(adata, resolution=leiden_resolution)
    analyze_marker_genes(adata)
    normalize_marker_expression(adata, pericyte_markers, ref_adata, fibroblast_label)
    plot_clusters_and_markers(adata, marker_genes)
    return adata


def main():
    adata = sc.read_h5ad('../_m/pericyte.hlca_core.dataset.h5ad')
    adata = subcluster_pericytes(adata)
    adata.write('pericyte.hlca_core.subclustered.h5ad')


if __name__ == "__main__":
    main()
