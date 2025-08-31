import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import gseapy as gp
from os import makedirs, path
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import NMF, FactorAnalysis
from scipy.stats import spearmanr, pointbiserialr, kruskal

RANDOM_SEED = 13

def load_adata(filepath):
    """Load AnnData object from .h5ad file."""
    return sc.read_h5ad(filepath)


def subset_adata_by_cluster(adata, cluster_key, cluster_value):
    """Subset AnnData object for specific cluster."""
    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster key '{cluster_key}' not found in adata.obs")
    subset = adata[adata.obs[cluster_key] == cluster_value].copy()
    if subset.n_obs == 0:
        raise ValueError(f"No cells found for cluster '{cluster_value}' with key '{cluster_key}'")
    print(f"[INFO] Subset to cluster '{cluster_value}' with {subset.n_obs} cells")
    return subset


def apply_dimred_gene_programs(adata, n_components=10, method='nmf', tune=False,
                               candidate_range=range(5, 26), random_state=13, max_iter=500):
    """
    Apply NMF or factor analysis for gene program extraction, with optional tuning of
    n_components.
    """
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    best_score = None
    best_n = n_components
    best_model = None

    if tune:
        print(f"Tuning n_components for {method.upper()}...")
        for k in candidate_range:
            if method == 'nmf':
                model = NMF(n_components=k, init='nndsvda',
                            random_state=random_state, max_iter=max_iter)
                W = model.fit_transform(X)
                H = model.components_
                score = model.reconstruction_err_  # lower is better
            elif method == 'fa':
                model = FactorAnalysis(n_components=k, random_state=random_state)
                model.fit(X)
                # Use log-likelihood or explained variance proxy
                score = model.score(X)  # higher is better
            else:
                raise ValueError("Specify method as 'nmf' or 'fa'")

            if best_score is None:
                best_score, best_n, best_model = score, k, model
            else:
                if method == 'nmf' and score < best_score:  # minimize error
                    best_score, best_n, best_model = score, k, model
                elif method == 'fa' and score > best_score:  # maximize likelihood
                    best_score, best_n, best_model = score, k, model
        print(f"Best n_components for {method.upper()}: {best_n} (score={best_score:.4f})")
    else:
        if method == 'nmf':
            best_model = NMF(n_components=n_components, init='nndsvda',
                             random_state=random_state, max_iter=max_iter)
        elif method == 'fa':
            best_model = FactorAnalysis(n_components=n_components,
                                        random_state=random_state)
        else:
            raise ValueError("Specify method as 'nmf' or 'fa'")

        best_model.fit(X)

    # Final decomposition using best model
    W = best_model.fit_transform(X)
    H = best_model.components_
    adata.obsm[f'{method}_cell_usage'] = W
    adata.varm[f'{method}_gene_loadings'] = H.T
    
    return W, H, best_n


def compare_factors(nmf_H, fa_H, cluster_value):
    """Compare NMF and FA gene loadings via Spearman correlation to find robust programs."""
    results = []
    for i in range(nmf_H.shape[0]):
        for j in range(fa_H.shape[0]):
            corr, pval = spearmanr(nmf_H[i, :], fa_H[j, :])
            results.append({'cluster': cluster_value, 'NMF_factor': i, 'FA_factor': j,
                            'Spearman_corr': corr, 'Spearman_pval': pval})

    df = pd.DataFrame(results).sort_values(by='Spearman_corr', ascending=False)
    n_feat = df[(np.abs(df['Spearman_corr']) >= 0.5)].shape[0]
    print(f"[INFO] Cluster {cluster_value}: Found {n_feat} factor pairs with abs correlation >= 0.5")
    return df


def run_go_enrichment(genes, background_genes=None, outdir=None):
    """Run gene ontology enrichment on a list of genes using gseapy."""
    gene_sets = ['GO_Biological_Process_2021', 'GO_Molecular_Function_2021',
                 'KEGG_2021_Human']
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        background=background_genes,
        cutoff=0.05,
        outdir=outdir,
        no_plot=True
    )

    if enr is None or enr.results.empty:
        print("[WARN] No GO enrichment results found.")
        return pd.DataFrame()

    print(f"[INFO] GO enrichment found {len(enr.results)} significant terms")
    return enr.results


def correlate_program_activity_with_phenotypes(
        adata, method, cluster_value, phenotypes):
    """Correlate gene program activities with phenotype variables."""
    usage = adata.obsm[f'{method}_cell_usage']
    records = []

    for i in range(usage.shape[1]):
        program_scores = usage[:, i]

        for phenotype, values in phenotypes.items():
            # Clean data
            values = pd.Series(values)

            valid_mask = values.notnull() & ~values.astype(str).str.lower().isin(['unknown', 'na', 'nan', ''])
            filtered_scores = program_scores[valid_mask.values]
            filtered_values = values.loc[valid_mask]

            if filtered_values.empty:
                continue

            # Determine data type for phenotype
            unique_vals = filtered_values.unique()

            if hasattr(filtered_values.dtype, "categories"):
                try:
                    filtered_values = filtered_values.cat.codes
                    mask = filtered_values != -1
                    filtered_scores = filtered_scores[mask]
                    filtered_values = filtered_values[mask]
                    unique_vals = filtered_values.unique()
                except Exception as e:
                    print(f"[WARN] Failed to convert categorical codes for phenotype {phenotype} in cluster {cluster_value}, program {i}: {e}")
                    continue
            try:
                if np.issubdtype(filtered_values.dtype, np.number):
                    # Numeric phenotype
                    if len(unique_vals) == 1:
                        continue  # no variability
                    if len(unique_vals) == 2:
                        # Binary numeric, use point biserial
                        corr, pval = pointbiserialr(filtered_values, filtered_scores)
                    elif len(unique_vals) > 2:
                        # Continuous or multiple numeric values: Spearman
                        corr, pval = spearmanr(filtered_scores, filtered_values)
                    else:
                        continue
                else:
                    # For non-numeric (should not happen after conversion), skip
                    continue
            except Exception as e:
                print(f"[WARN] Correlation failed for phenotype {phenotype} in cluster {cluster_value}, program {i}: {e}")
                continue
            
            records.append({
                'cluster': cluster_value,
                'program': i,
                'method': method,
                'phenotype': phenotype,
                'corr': corr,
                'pval': pval
            })
    df = pd.DataFrame(records)
    df['pval_adj'] = np.nan
    for pheno, group in df.groupby(["phenotype"]):
        pvals = group["pval"].values
        _, fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
        df.loc[group.index, 'pval_adj'] = fdr
    print(f"[INFO] Cluster {cluster_value}: Correlation results calculated")
    return df


def get_cluster_values(adata, cluster_key):
    """Extract the unique cluster values to test."""
    return adata.obs[cluster_key].unique()


def main(adata, cluster_key, outdir, model):
    cluster_values = get_cluster_values(adata, cluster_key)

    all_factor_comp = []
    all_corrs = []
    all_go_res = []
    
    for c_value in cluster_values:
        adata_sub = subset_adata_by_cluster(adata, cluster_key, c_value)
        gene_names = adata_sub.var.feature_name.to_list()

        # Gene program extraction
        nmf_w, nmf_h, nmf_n = apply_dimred_gene_programs(
            adata_sub, method='nmf', tune=True, max_iter=2_500
        )
        fa_w, fa_h, _ = apply_dimred_gene_programs(
            adata_sub, n_components=nmf_n, method='fa'
        )

        # Factor comparison store results
        factor_comparison = compare_factors(nmf_h, fa_h, c_value)
        all_factor_comp.append(factor_comparison)
        
        # Phenotype correlations
        phenotypes = {}
        pheno_cols = ['sex', 'smoking_status', 'self_reported_ethnicity', 'age_or_mean_of_age_range']
        if model == 'full':
            pheno_cols.append('disease')

        for pheno in pheno_cols:
            phenotypes[pheno] = adata_sub.obs[pheno].values

        nmf_corr = correlate_program_activity_with_phenotypes(
            adata_sub, 'nmf', c_value, phenotypes)
        fa_corr = correlate_program_activity_with_phenotypes(
            adata_sub, 'fa', c_value, phenotypes)
        corr_df = pd.concat([nmf_corr, fa_corr], axis=0)
        all_corrs.append(corr_df)

        # GO enrichment
        temp_go_res = []
        for i in range(nmf_n)):
            top_genes = [gene_names[idx] for idx in np.argsort(nmf_h[i, :])[::-1][:100]]
            go_results = run_go_enrichment(top_genes, background_genes=gene_names,
                                           outdir=None)
            go_results["NMF_Module"] = i+1
            go_results["cluster"] = c_value
            temp_go_res.append(go_results)
        all_go_res.append(pd.concat(temp_go_res, ignore_index=True))

    # Combine data
    factor_comp_df = pd.concat(all_factor_comp, ignore_index=True)
    corrs_df = pd.concat(all_corrs, ignore_index=True)
    go_res_df = pd.concat(all_go_res, ignore_index=True)

    # Save data
    if outdir and not path.exists(outdir):
        makedirs(outdir, exist_ok=True)

    factor_comp_df.to_csv(path.join(outdir, "factor_compare.all_clusters.tsv"),
                          sep="\t", index=False)
    nmf_corrs_df.to_csv(path.join(outdir, "nmf_pheno_corr.all_clusters.tsv"),
                        sep="\t", index=False)
    fa_corrs_df.to_csv(path.join(outdir, "fa_pheno_corr.all_clusters.tsv"),
                       sep="\t", index=False)
    go_res_df.to_csv(path.join(outdir, "go_results.all_clusters.tsv"),
                     sep="\t", index=False)

    # Session information
    session_info.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster-specific gene program extraction and analysis with NMF and Factor Analysis"
    )
    parser.add_argument("--model", type=str, default="core",
                        help="Model type: 'core' or 'full'. Default: core")
    parser.add_argument('--cluster_key', type=str, default='leiden',
                        help="Key in adata.obs for cluster labels")
    parser.add_argument('--output_dir', type=str, default='results',
                        help="Directory to save GO enrichment results")
    args = parser.parse_args()

    # Load data
    fname = f'../../_m/pericyte.hlca_{args.model}.subclustered.analysis.h5ad'
    adata = sc.read_h5ad(fname)
    outdir = path.join(args.output_dir, args.model)
    main(adata, args.cluster_key, outdir)
