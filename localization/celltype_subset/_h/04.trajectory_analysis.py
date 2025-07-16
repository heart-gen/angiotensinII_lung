import scanpy as sc
import phate
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF, FactorAnalysis
import gseapy as gp

# Optional: import TFvelo (ensure it's installed)
try:
    import tfvelo
except ImportError:
    tfvelo = None
    print("TFvelo not installed; RNA velocity step will be skipped.")

def load_adata(filepath):
    '''Load AnnData object from .h5ad file'''
    return sc.read_h5ad(filepath)

def compute_pseudotime(adata, root_cell=None, phate_key='X_phate'):
    '''Compute diffusion pseudotime using PHATE embedding and set root'''
    sc.pp.neighbors(adata, use_rep=phate_key)
    if root_cell is None:
        root_cell = np.argmin(adata.obsm[phate_key][:, 0])
    sc.tl.dpt(adata, root=root_cell)
    print("Pseudotime added to adata.obs['dpt_pseudotime']")
    return adata

def run_tfvelo(adata, phate_key='X_phate', cluster_key='leiden'):
    '''Estimate RNA velocity with TFvelo and embed results on PHATE space'''
    if tfvelo is None:
        print("TFvelo not installed.")
        return
    tfvelo.tl.velocity(adata)
    tfvelo.tl.velocity_graph(adata)
    tfvelo.pl.velocity_embedding_stream(adata, basis=phate_key, color=cluster_key)
    print("RNA velocity analysis complete.")

def apply_dimred_gene_programs(adata, n_components=10, method='nmf'):
    '''Apply NMF or factor analysis for gene program extraction'''
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    if method == 'nmf':
        model = NMF(n_components=n_components, init='nndsvda', random_state=42)
    elif method == 'fa':
        model = FactorAnalysis(n_components=n_components, random_state=42)
    else:
        raise ValueError("Specify method as 'nmf' or 'fa'")
    W = model.fit_transform(X)
    H = model.components_
    adata.obsm[f'{method}_cell_usage'] = W
    adata.varm[f'{method}_gene_loadings'] = H.T
    print(f"{method.upper()} complete: cell usage in adata.obsm, gene loadings in adata.varm")
    return W, H

def perform_gsea(adata, groupby='leiden', top_n=100, gene_sets='GO_Biological_Process_2021'):
    '''Run GSEA for each Leiden cluster using marker genes'''
    sc.tl.rank_genes_groups(adata, groupby=groupby, method='wilcoxon')
    enrichment_results = {}
    for cluster in adata.obs[groupby].unique():
        genes = sc.get.rank_genes_groups_df(adata, group=cluster)['names'].head(top_n).tolist()
        enr = gp.enrichr(gene_list=genes, gene_sets=gene_sets, organism='Human', outdir=None, no_plot=True)
        enrichment_results[cluster] = enr.results
        print(f"GSEA complete for cluster {cluster}")
    return enrichment_results

def main():
    # 1. Load data
    adata = load_adata('pericyte.hlca_core.subclustered.h5ad')

    # 2. Trajectory inference & pseudotime
    adata = compute_pseudotime(adata)

    # 3. RNA velocity (TFvelo)
    run_tfvelo(adata)

    # 4. Gene program extraction
    nmf_W, nmf_H = apply_dimred_gene_programs(adata, n_components=10, method='nmf')
    fa_W, fa_H = apply_dimred_gene_programs(adata, n_components=10, method='fa')

    # 5. GSEA per Leiden cluster
    gsea_results = perform_gsea(adata)

    # 6. Save processed data
    adata.write('pericyte.hlca_core.subclustered.analysis.h5ad')

if __name__ == '__main__':
    main()
