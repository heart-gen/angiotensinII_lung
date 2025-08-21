import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from os import makedirs, path
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF, FactorAnalysis

def load_adata(filepath):
    '''Load AnnData object from .h5ad file'''
    return sc.read_h5ad(filepath)


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
    # Gene program extraction
    nmf_W, nmf_H = apply_dimred_gene_programs(adata, n_components=10,
                                              method='nmf')
    fa_W, fa_H = apply_dimred_gene_programs(adata, n_components=10,
                                            method='fa')
    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
