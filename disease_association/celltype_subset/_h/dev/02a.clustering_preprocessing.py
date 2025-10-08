"""
Single-Cell Data Clustering and Visualization

Preprocessing step. Broken up to use with Bridges-2.
"""
import session_info
import scanpy as sc

def load_query_data():
    adata = sc.read_h5ad("../../_m/integrated_data.h5ad")
    adata.obs["celltype"] = "unknown"
    return adata


def load_reference():
    adata = sc.read_h5ad("organoid_data.h5ad")
    adata.X = adata.X.tocsr()
    adata.obs["celltype"] = adata.obs["FinalName"]
    return adata


def main():
    # Load data
    query_adata = load_query_data()
    ref_adata = load_reference()
    ref_adata.layers["counts"] = ref_adata.X

    # Ensure common genes
    common_genes = ref_adata.var_names.intersection(query_adata.var_names)
    ref_adata = ref_adata[:, common_genes].copy()
    query_adata = query_adata[:, common_genes].copy()

    # Preprocessing on CPU
    sc.pp.normalize_total(ref_adata, target_sum=1e4)
    sc.pp.log1p(ref_adata)
    sc.pp.highly_variable_genes(ref_adata, n_top_genes=2000)

    hvgs = ref_adata.var.highly_variable
    ref_hvg = ref_adata[:, hvgs].copy()
    query_hvg = query_adata[:, hvgs].copy()

    # Save preprocessed objects for GPU step
    ref_hvg.write_h5ad("ref_hvg.h5ad", compression="gzip")
    query_hvg.write_h5ad("query_hvg.h5ad", compression="gzip")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
