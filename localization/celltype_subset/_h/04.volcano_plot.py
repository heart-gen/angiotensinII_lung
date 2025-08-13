## This script generates volcano plots
import re
import numpy as np
import session_info
import scanpy as sc
import pandas as pd
from os import makedirs, path
import matplotlib.pyplot as plt
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests

def load_ensembl_to_geneid_map(adata_file_path):
    """
    Load a mapping table with columns:
        ensembl_id   gene_id
    Returns a dictionary {ensembl_id: gene_id}.
    """
    adata = sc.read_h5ad(adata_file_path)
    mapping_df = pd.DataFrame(adata.var)\
                   .reset_index().rename(columns={"index":"ensembl_id"})
    return dict(zip(mapping_df['ensembl_id'], mapping_df['feature_name']))


def load_all_clusters(filepath, mapping_file=None):
    """
    Load a precomputed rank_genes_groups export (wide format) into
    a dict of DataFrames, keyed by cluster ID.
    If mapping_file is provided, convert Ensembl IDs to gene symbols.
    """
    df = pd.read_csv(filepath, sep="\t", compression="infer", index_col=0)

    # Map gene ids, optional
    ensembl_to_gene = {}
    if mapping_file:
        ensembl_to_gene = load_ensembl_to_geneid_map(mapping_file)

    # Extract unique cluster IDs from column names by regex
    cluster_ids = sorted({
        re.match(r"^(\d+)_names$", col).group(1)
        for col in df.columns if re.match(r"^\d+_names$", col)
    })

    clusters_data = {}
    for cid in cluster_ids:
        name_col = f"{cid}_names"
        pval_col = f"{cid}_pvals"
        logfc_col = f"{cid}_logfoldchanges"

        cluster_df = pd.DataFrame({
            'gene_id': df[name_col],
            'pvals': df[pval_col],
            'logfoldchanges': df[logfc_col]
        })

        # Convert Ensembl IDs → gene symbols if mapping provided
        if ensembl_to_gene:
            cluster_df['gene_name'] = cluster_df['gene_id'].map(
                lambda x: ensembl_to_gene.get(x, x)
            )

        # Calculate adjusted p-values (FDR Benjamini-Hochberg)
        cluster_df['pvals_adj'] = multipletests(
            cluster_df['pvals'], method='fdr_bh'
        )[1]
        cluster_df['-log10_fdr'] = -np.log10(cluster_df['pvals_adj'] + 1e-300)
        clusters_data[cid] = cluster_df

    return clusters_data


def plot_volcano(cluster_df, cluster_id, pval_cutoff=0.05, logfc_cutoff=1.0,
                 top_n_labels=10, outdir="volcano_plots"):
    """
    Make a volcano plot for one cluster and save to file.
    """
    makedirs(outdir, exist_ok=True)
    sig_mask = (cluster_df['pvals_adj'] < pval_cutoff) & \
               (np.abs(cluster_df['logfoldchanges']) > logfc_cutoff)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        cluster_df['logfoldchanges'], cluster_df['-log10_fdr'],
        c=sig_mask, cmap='coolwarm', s=10, alpha=0.7, edgecolor='none'
    )
    plt.axhline(-np.log10(pval_cutoff), color='grey', ls='--', lw=1)
    plt.axvline(logfc_cutoff, color='grey', ls='--', lw=1)
    plt.axvline(-logfc_cutoff, color='grey', ls='--', lw=1)
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10(FDR)')
    plt.title(f'Volcano Plot — Cluster {cluster_id}')

    # Label top significant genes
    top_hits = cluster_df[sig_mask].nlargest(top_n_labels, '-log10_fdr')
    texts = []
    for _, row in top_hits.iterrows():
        texts.append(
            plt.text(row['logfoldchanges'], row['-log10_fdr'], row['gene_name'],
                 fontsize=7, ha='right')
        )

    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

    outfile = path.join(outdir, f"volcano_cluster_{cluster_id}.png")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)

    pdf_outfile = outfile.replace(".png", ".pdf")
    plt.savefig(pdf_outfile)

    plt.close()
    print(f"[OK] Saved {outfile} and {pdf_outfile}")


def main():
    # File paths
    filepath = "rank_genes_groups_results.tsv.gz"
    adata_file_path = "pericyte.hlca_core.subclustered.h5ad"

    # Load and process cluster data
    clusters_data = load_all_clusters(filepath, mapping_file=adata_file_path)

    # Generate volcano plots
    for cid, cdata in clusters_data.items():
        plot_volcano(cdata, cluster_id=cid, pval_cutoff=0.05, logfc_cutoff=1.0)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
