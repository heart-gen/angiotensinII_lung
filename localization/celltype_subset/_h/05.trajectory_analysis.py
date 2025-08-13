## This script performs trajectory analysis
import phate
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from sklearn.decomposition import NMF, FactorAnalysis

def load_adata(filepath):
    '''Load AnnData object from .h5ad file'''
    return sc.read_h5ad(filepath)


def compute_pseudotime(adata, root_cell=None, phate_key='X_phate'):
    '''Compute diffusion pseudotime using PHATE embedding and set root'''
    sc.pp.neighbors(adata, use_rep=phate_key, random_state=13)
    # PAGA
    sc.tl.draw_graph(adata, random_state=13)
    sc.tl.paga(adata, groups="leiden")
    if root_cell is None:
        root_cell = np.argmin(adata.obsm[phate_key][:, 0])
        adata.uns['iroot'] = root_cell
    sc.tl.dpt(adata)
    print("Pseudotime added to adata.obs['dpt_pseudotime']")
    return adata


def apply_dimred_gene_programs(adata, n_components=10, method='nmf'):
    '''Apply NMF or factor analysis for gene program extraction'''
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    if method == 'nmf':
        model = NMF(n_components=n_components, init='nndsvda', random_state=13)
    elif method == 'fa':
        model = FactorAnalysis(n_components=n_components, random_state=13)
    else:
        raise ValueError("Specify method as 'nmf' or 'fa'")
    W = model.fit_transform(X)
    H = model.components_
    adata.obsm[f'{method}_cell_usage'] = W
    adata.varm[f'{method}_gene_loadings'] = H.T
    print(f"{method.upper()} complete: cell usage in adata.obsm, gene loadings in adata.varm")
    return W, H


def main():
    adata = load_adata('../_m/pericyte.hlca_core.subclustered.h5ad')
    # Trajectory inference & pseudotime
    adata = compute_pseudotime(adata)
    # Gene program extraction
    nmf_W, nmf_H = apply_dimred_gene_programs(adata, n_components=10, method='nmf')
    fa_W, fa_H = apply_dimred_gene_programs(adata, n_components=10, method='fa')
    # Save processed data
    adata.write('pericyte.hlca_core.subclustered.analysis.h5ad')
    # Session information
    session_info.show()


if __name__ == '__main__':
    main()
