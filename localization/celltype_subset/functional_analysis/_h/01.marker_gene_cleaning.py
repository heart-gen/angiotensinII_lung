import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import gseapy as gp
import seaborn as sns
import matplotlib.pyplot as plt
from os import makedirs, path, listdir

def load_adata(filepath):
    """Load AnnData object from .h5ad file."""
    return sc.read_h5ad(filepath)


def split_rank_results(res_file, adata, outdir="cluster_markers", uniq_thresh=None):
    """
    Split wide-form rank_genes_groups results into per-cluster files.
    Produces two outputs per cluster: all DEGs and unique DEGs.
    """
    # Load wide rank_genes_groups results
    df = pd.read_csv(res_file, sep="\t", index_col=0)
    makedirs(outdir, exist_ok=True)

    # Build ensembl_id to gene name
    ensg_to_gene = dict(zip(adata.var.index, adata.var["feature_name"]))

    # Identify clusters
    clusters = sorted({col.split("_")[0] for col in df.columns})

    # Build combined long-form dataframe
    results = []
    for cl in clusters:
        sub = pd.DataFrame({
            "cluster": cl,
            "gene_id": df[f"{cl}_names"],
            "pval": df[f"{cl}_pvals"],
            "logFC": df[f"{cl}_logfoldchanges"],
            })
        results.append(sub)
    combined = pd.concat(results, ignore_index=True)

    # Map ensembl_id to gene name
    combined["gene_name"] = combined["gene_id"].map(ensg_to_gene).fillna(combined["gene_id"])

    # Filter out ribosomal and mitochondrial
    mask_ribo = combined["gene_name"].str.match(r"^RP[SL]")
    mask_mt = combined["gene_name"].str.match(r"^MT-")
    combined = combined[~(mask_ribo | mask_mt)]

    # Determine uniqueness threshold
    if uniq_thresh is None:
        uniq_thresh = len(clusters) // 2

    sig = combined[combined["pval"] < 0.05].copy()
    counts = sig.groupby("gene_id")["cluster"].nunique()
    mcluster_genes = counts[counts > uniq_thresh].index

    # Export per cluster
    for cl in clusters:
        sub = combined[combined["clusterr"] == cl].sort_values("pval")
        # All DEGs
        sub.to_csv(path.join(outdir, f"cluster_{cl}.all.tsv"),
                   sep="\t", index=False)
        # Unique DEGs (remove ubiquitous genes)
        sub_unique = sub[~sub["gene_id"].isin(mcluster_genes)].copy()
        sub_unique.to_csv(path.join(outdir, f"cluster_{cl}.unique.tsv"),
                          sep="\t", index=False)
        print(f"Saved per-cluster 'all' and 'unique' marker files to {outdir}")


def pathway_heatmap(marker_dir, outdir="figures", prefix="pathway_heatmap"):
    """
    Run pathway enrichment per cluster and generate a heatmap of -log10(FDR).
    """
    # Gene sets
    gene_sets=['GO_Biological_Process_2021',
               'GO_Molecular_Function_2021',
               'KEGG_2021_Human']
    filenames = [path.join(marker_dir, f) for f in listdir(marker_dir) 
                 if f.startswith("cluster_") and f.endswith(".tsv")]
    filenames = sorted(filenames)

    enrich_mat = {}
    for fn in filenames:
        # Load data
        cl = path.basename(fn).split(".")[0].split("_")[1]
        df = pd.read_csv(fn, sep="\t")

        # Select significant genes
        genes = df.loc[df["pval"] < 0.05, "gene_name"].dropna().astype(str).tolist()
        if not genes:
            continue

        # Run enrichment
        enr = gp.enrichr(gene_list=genes, gene_sets=gene_sets, organism="Human")
        if enr.results is not None and not enr.results.empty:
            enr.results["-log10FDR"] = -np.log10(enr.results["Adjusted P-value"].clip(lower=1e-300))
            enr_unique = enr.results.drop_duplicates(subset=["Term"]).set_index("Term")
            enrich_mat[cl] = enr_unique["-log10FDR"]

    if not enrich_mat:
        print("No enrichment results found.")
        return

    # Combine into one matrix
    enrich_df = pd.DataFrame(enrich_mat).fillna(0)

    # Make output directory
    makedirs(outdir, exist_ok=True)
    
    # Plot heatmap
    g = sns.clustermap(enrich_df, cmap="viridis", figsize=(10,12))
    plt.title("Pathway Enrichment Heatmap")

    # Save plot
    png_file = path.join(outdir, f"{prefix}.png")
    pdf_file = path.join(outdir, f"{prefix}.pdf")
    g.savefig(png_file, dpi=300, bbox_inches="tight")
    g.savefig(pdf_file, bbox_inches="tight")
    plt.close()
    print(f"Saved heatmap to {png_file} and {pdf_file}")


def main(model):
    # General outputs
    marker_dir = "cluster_markers"
    outdir = path.join(marker_dir, model)
    adata_file = f"../../_m/pericyte.hlca_{model}.subclustered.analysis.h5ad"
    res_file = f"../../_m/figures/{model}/rank_genes_groups_results.txt.gz"

    # Load data
    adata = load_adata(adata_file)

    # Split cluster
    split_rank_results(res_file, adata, outdir)

    # Pathway heatmap
    pathway_heatmap(marker_dir, outdir=path.join("figures", model),
                    "pathway_heatmap")

    # Session information
    session_info.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster-specific gene program extraction and analysis with NMF and Factor Analysis"
    )
    parser.add_argument("--model", type=str, default="core",
                        help="Model type: 'core' or 'full'. Default: core")
    args = parser.parse_args()

    # Load data
    fname = f'../../_m/pericyte.hlca_{args.model}.subclustered.analysis.h5ad'
    adata = sc.read_h5ad(fname)
    outdir = path.join(args.output_dir, args.model)
    main(adata, args.cluster_key, outdir)
