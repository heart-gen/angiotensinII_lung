## Analyze Lung SmartSeq2 scRNA-seq data for angiotensin II.

import torch
import numpy as np
import pandas as pd
import scanpy as sc
import scvi, session_info
import matplotlib.pyplot as plt
from functools import lru_cache
from scipy.sparse import csr_matrix

@lru_cache()
def load_data():
    fn = "../_m/lungmap_dataset.h5ad"
    return sc.read_h5ad(fn)


@lru_cache()
def preprocessing_data():
    adata = load_data()
    # Select genes
    scvi.data.poisson_gene_selection(adata)
    adata = adata[:, adata.var["highly_variable"]]
    adata.layers["counts"] = csr_matrix(adata.X).copy()
    return adata


def train_model(adata):
    scvi.model.SCVI.setup_anndata(adata, layer="counts",
                                  batch_key="patient")
    torch.set_float32_matmul_precision('high')
    # define and train model
    model = scvi.model.SCVI(adata, gene_likelihood="nb")
    model.train(check_val_every_n_epoch=1,
                max_epochs=400, early_stopping=True,
                early_stopping_patience=20,
                early_stopping_monitor="elbo_validation",)
    # ensure convergence
    train_test_results = model.history["elbo_train"]
    train_test_results["elbo_validation"] = model.history["elbo_validation"]
    train_test_results.iloc[10:].plot(logy=True)
    plt.savefig("convergence_elbo.pdf")
    # latent space and UMAP
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.pl.umap(adata, color="cell_type",
               save=".annotated_cluster.batch_corrected.pdf")
    # save model
    model.save("scVI_model/", overwrite=True)


def diff_expr(adata):
    # load model
    model = scvi.model.SCVI.load("scVI_model/", adata=adata,
                                 use_gpu=True)
    cats = adata.obs.cell_type.cat.categories
    # get latent variables
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    # de analysis
    de_df = model.differential_expression(groupby="cell_type")
    # extract top markers for each cluster
    markers1 = {}; markers2 = {};
    for i, c in enumerate(cats):
        cid = f"{c} vs Rest"
        cell_type_df = de_df.loc[de_df.comparison == cid]
        cell_type_df = cell_type_df[cell_type_df.lfc_mean > 0]
        cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 3]
        cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]
        markers1[c] = cell_type_df.index.tolist()[:3]
        markers2[c] = cell_type_df.index.tolist()[:10] # top 10
    pd.DataFrame.from_records(markers2)\
                .to_csv("celltype_markers_top10.csv", index=False)
    sc.tl.dendrogram(adata, groupby="cell_type", use_rep="X_scVI")
    # plot dotplot
    sc.pl.dotplot(adata, markers1, groupby="cell_type",
                  dendrogram=True, color_map="Blues", swap_axes=True,
                  save="celltype.markers.pdf")


def main():
    ## Main
    seed = 20230706; np.random.seed(seed)
    adata = preprocessing_data()
    train_model(adata)
    diff_expr(adata)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()
