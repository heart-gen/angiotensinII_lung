## This script subclusters pericytes

import phate
import numpy as np
import scanpy as sc
import session_info
import matplotlib.pyplot as plt

def subcluster_and_analyze(adata, fibroblast_label='Fibroblasts',
                           pericyte_markers=['HIGD1B', 'PDGFRB', 'CSPG4'],
                           marker_genes=['AGTR1', 'ACTA2'], phate_knn=5,
                           phate_decay=40, leiden_resolution=1.0):
    # Step 4: Dimensionality reduction
    if 'X_pca' not in adata.obsm:
        sc.tl.pca(adata, svd_solver='arpack')
    
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)
    
    # PHATE embedding for trajectory analysis
    phate_op = phate.PHATE(knn=phate_knn, decay=phate_decay, n_jobs=-1)
    phate_emb = phate_op.fit_transform(adata.obsm['X_pca'])
    adata.obsm['X_phate'] = phate_emb

    # Step 5: Clustering
    sc.tl.leiden(adata, resolution=leiden_resolution)

    # Step 6: Marker gene analysis
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

    # Step 8: Visualization

    # Plot UMAP colored by Leiden clusters
    sc.pl.umap(adata, color='leiden', title='Leiden Clusters')

    # Plot PHATE colored by Leiden clusters
    sc.pl.embedding(adata, basis='X_phate', color='leiden', title='PHATE Clusters')

    # Plot expression of key marker genes on UMAP and PHATE
    for gene in marker_genes:
        if gene in adata.var_names:
            sc.pl.umap(adata, color=gene, title=f'UMAP: {gene} expression')
            sc.pl.embedding(adata, basis='X_phate', color=gene, title=f'PHATE: {gene} expression')
        else:
            print(f"Warning: Gene {gene} not found in dataset")

    # Normalize pericyte marker expression to fibroblast count
    # Calculate mean expression of pericyte markers per cell
    pericyte_expr = adata[:, pericyte_markers].X.mean(axis=1)

    # Identify fibroblast cells
    fibroblast_cells = adata.obs['cell_type'] == fibroblast_label
    fibroblast_count = fibroblast_cells.sum()

    if fibroblast_count > 0:
        normalized_expr = pericyte_expr / fibroblast_count
    else:
        print("Warning: No fibroblast cells found; skipping normalization")
        normalized_expr = pericyte_expr

    # Add normalized expression to adata.obs for visualization
    adata.obs['normalized_pericyte_marker_expr'] = normalized_expr

    # Plot normalized pericyte marker expression on UMAP
    sc.pl.umap(adata, color='normalized_pericyte_marker_expr', title='Normalized Pericyte Marker Expression')
    return adata


def main():
    adata = sc.read_h5ad('../_m/pericyte.hlca_core.dataset.h5ad')
    adata = subcluster_and_analyze(adata)
    adata.write('pericyte.hlca_core.subclustered.h5ad')
