"""
Pericyte state classification using airspace proximity and marker signatures.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from pathlib import Path
from anndata import AnnData
import logging, argparse, yaml
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from typing import Sequence, Dict, List
from statsmodels.stats.anova import anova_lm
from scipy.sparse import issparse, csr_matrix
from sklearn.neighbors import NearestNeighbors
from scipy.stats import f_oneway, kruskal, rankdata, chi2_contingency
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sns.set_context("talk")
sns.set_style("whitegrid")

def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata", type=Path, required=True,
                        help="Input AnnData with airspace score")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Output directory for QC and updated AnnData")
    parser.add_argument("--state-yaml", required=True, type=Path,
                        help="YAML with pericyte states genes")
    parser.add_argument("--use-rep", type=str, default="X_pca_harmony",
                        help="obsm key for latent representation (default: X_pca_harmony)")
    parser.add_argument("--cluster-key", type=str, default="subclusters",
                        help="obs column with fine cell types (default: subclusters)")
    parser.add_argument("--pericyte-label", type=str, default="Pericytes",
                        help="Value in cluster-key used for pericytes")
    parser.add_argument("--ec-labels", nargs="+",
                        default=["EC aerocyte capillary", "EC general capillary"],
                        help="Labels treated as alveolar capillary EC")
    parser.add_argument("--neighbors", type=int, default=50,
                        help="n_neighbors for graph construction if neighbors not present")
    parser.add_argument("--seed", type=int, default=13, help="Random seed for neighbors")
    return parser.parse_args()


def load_anndata(path: Path) -> AnnData:
    logging.info("Loading AnnData from %s", path)
    adata = sc.read_h5ad(path)
    return adata


def ensure_neighbors(adata: AnnData, use_rep: str, n_neighbors: int, seed: int) -> None:
    """Ensure that a neighbor graph + distances exist for the chosen representation."""
    if "neighbors" in adata.uns and adata.uns["neighbors"].get("use_rep", None) == use_rep:
        logging.info("Neighbor graph already present for %s; reusing.", use_rep)
        return

    logging.info(
        "Computing neighbors for %s (n_neighbors=%d, random_state=%d)",
        use_rep, n_neighbors, seed
    )
    sc.pp.neighbors(
        adata, use_rep=use_rep, n_neighbors=n_neighbors,
        metric="cosine", random_state=seed,
    )


def pericyte_to_nearest_ec_latent(
    adata: AnnData, use_rep: str="X_pca_harmony", cluster_key: str="subclusters",
    pericyte_label: str="Pericytes", n_neighbors: int=1, metric: str="euclidean",
    ec_labels: Sequence[str]=("EC aerocyte capillary", "EC general capillary"),
) -> np.ndarray:
    if use_rep not in adata.obsm_keys():
        raise KeyError(f"{use_rep!r} not found in adata.obsm")

    Z = adata.obsm[use_rep]
    labs = adata.obs[cluster_key].astype(str)
    mask_peri = labs.eq(pericyte_label).values
    mask_ec   = labs.isin(ec_labels).values

    n_cells = adata.n_obs
    n_per   = int(mask_peri.sum())
    n_ec    = int(mask_ec.sum())

    if n_ec == 0:
        raise ValueError(f"No EC cells found for labels={ec_labels}")

    logging.info(
        "[NN-latent] use_rep=%s | metric=%s | n_cells=%d | pericytes=%d | EC=%d | k=%d",
        use_rep, metric, n_cells, n_per, n_ec, n_neighbors
    )
    
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric=metric).fit(Z[mask_ec])
    dist, _ = nbrs.kneighbors(Z[mask_peri], n_neighbors=n_neighbors,
                              return_distance=True)
    dmin = dist[:, 0].astype(np.float32)
    pericyte_to_ec = np.full(adata.n_obs, np.nan, dtype=np.float32)
    pericyte_to_ec[mask_peri] = dmin
    return pericyte_to_ec


def load_state_genes(path: Path) -> Dict[str, List[str]]:
    """Load pericyte-state marker panels from YAML."""
    data = yaml.safe_load(path.read_text())
    if "pericyte_state" not in data:
        raise ValueError("YAML must contain top-level key 'pericyte_state'")
    return data["pericyte_state"]


def rank_module_score_chunked(
    adata: AnnData, genes: Sequence[str], score_name: str,
    layer: str = "logcounts", chunk_size: int = 1000,
):
    """
    AUCell-ish: per cell, rank genes; score = mean percentile of target genes.
    Memory-efficient gene module scoring using ranks and chunked processing.
    """
    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    var_names = np.asarray(adata.var_names)

    idx = np.where(np.isin(var_names, genes))[0]
    if len(idx) == 0:
        adata.obs[score_name] = np.nan
        return adata

    n_cells, n_genes = X.shape
    scores = np.full(n_cells, np.nan, dtype=np.float32)

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        X_chunk = X[start:end]
        X_chunk = X_chunk.toarray() if issparse(X_chunk) else np.asarray(X_chunk)

        r1    = np.argsort(X_chunk, axis=1).astype(np.int32)
        ranks = np.argsort(r1, axis=1).astype(np.float32)
        percentiles = ranks / (n_genes - 1 + 1e-9)
        scores[start:end] = percentiles[:, idx].mean(axis=1)

    adata.obs[score_name] = scores
    return adata


def score_signatures(
    adata: AnnData, states: Dict[str, Dict[str, str]],
    layer: str | None = "logcounts"):
    """
    Compute simple module scores for capillary vs arteriolar signatures.
    """
    cap = [g for g in states["capillary"] if g in adata.var_names]
    art = [g for g in states["arteriolar"] if g in adata.var_names]

    logging.info("Capillary signature genes present: %s", cap)
    logging.info("Arteriolar signature genes present: %s", art)

    rank_module_score_chunked(adata, cap, "capillary_sig", layer=layer)
    rank_module_score_chunked(adata, art, "arteriolar_sig", layer=layer)
    return adata


def per_donor_z(adata, cols, donor_key="donor_id"):
    """Z-score columns within donor to remove donor effects."""
    df = adata.obs
    for col in cols:
        x = df[col]
        df[col + "_z"] = (
            x.groupby(df[donor_key], observed=False)
             .transform(lambda v: (v - v.mean()) / (v.std(ddof=0) + 1e-8))
        )
    return adata


def assign_pericyte_states(
    adata: AnnData, cluster_key: str, pericyte_label: str,
    donor_key: str = "donor_id", margin=0.4,
) -> Dict[str, int]:
    """
    Combine z-scored distance + signature scores into a pericyte_state.
    """
    peri_mask = adata.obs[cluster_key].astype(str) == pericyte_label

    for col in ["capillary_sig", "arteriolar_sig", "airspace_dist"]:
        if col in adata.obs.columns:
            per_donor_z(adata, [col], donor_key=donor_key)

    ad = adata.obs
    ad["sig_diff_z"] = ad["capillary_sig_z"] - ad["arteriolar_sig_z"]
    # proximity feature: smaller distance -> more capillary-like
    prox = -ad["airspace_dist_z"]
    prox = prox.fillna(prox.median())
    sigd = ad["sig_diff_z"].fillna(0)

    # Simple weight
    w_sigd, w_prox = 1.0, 1.0
    combined = (w_sigd * sigd) + (w_prox * prox)
    ad["pericyte_combined_score"] = combined

    # Assign states only within pericytes
    state = pd.Series(index=ad.index, dtype=object)
    state[:] = np.nan
    state[(peri_mask) & (combined >=  margin)] = "capillary_like"
    state[(peri_mask) & (combined <= -margin)] = "arteriolar_like"
    state[(peri_mask) & state.isna()] = "ambiguous"
    adata.obs["pericyte_state"] = state

    return {
        s: ad.loc[(state == s) & peri_mask, donor_key].dropna().nunique()
        for s in ["capillary_like", "arteriolar_like", "ambiguous"]
    }


def plot_violin_by_state(
    adata: AnnData, cluster_key: str, pericyte_label: str, outdir: Path,
) -> None:
    peri = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()

    def _violin(col: str, fname: str, ylabel: str):
        plt.figure(figsize=(5, 5))
        sns.violinplot(
            data=peri, x="pericyte_state", y=col,
            order=["capillary_like", "arteriolar_like"], cut=0,
        )
        sns.stripplot(
            data=peri, x="pericyte_state", y=col,
            order=["capillary_like", "arteriolar_like"],
            color="black", size=2, alpha=0.4,
        )
        plt.ylabel(ylabel)
        plt.xlabel("Pericyte state")
        plt.tight_layout()
        plt.savefig(outdir / f"{fname}.png", dpi=300)
        plt.savefig(outdir / f"{fname}.pdf")
        plt.close()

    _violin("airspace_dist_z", "violin_airspace_dist_z", "Airspace distance (z)")
    _violin("sig_diff_z", "violin_sig_diff_z", "Capillary − arteriolar signature (z)")


def plot_batch_composition(
    adata: AnnData, outdir: Path, state_key: str = "pericyte_state",
    batch_keys: Sequence[str] = ("donor_id", "data", "assay", "tissue_sampling_method"),
    cluster_key: str = "subclusters", pericyte_label: str = "Pericytes",
) -> Dict[str, pd.DataFrame]:
    """Stacked bar plots of pericyte_state composition across batch/dataset."""
    peri = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    results: Dict[str, pd.DataFrame] = {}

    for key in batch_keys:
        if key not in peri.columns:
            continue
        tab = pd.crosstab(peri[key], peri[state_key], normalize="index")
        results[key] = tab

        plt.figure(figsize=(max(5, 1.2 * tab.shape[0]), 4))
        tab.plot(kind="bar", stacked=True, ax=plt.gca(), colormap="tab20")
        plt.ylabel("Proportion")
        plt.xlabel(key)
        plt.legend(title="Pericyte state", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(outdir / f"state_by_{key}.png", dpi=300)
        plt.savefig(outdir / f"state_by_{key}.pdf")
        plt.close()

    return results


def run_leiden_on_pericytes(
    adata: AnnData, use_rep: str, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", resolution: float = 0.5,
    n_neighbors: int = 50, leiden_key: str = "leiden_pericytes",
    seed: int = 13
) -> AnnData:
    """
    Run Leiden clustering on pericytes only, using given representation.
    """
    peri = adata[adata.obs[cluster_key] == pericyte_label].copy()
    sc.pp.neighbors(peri, use_rep=use_rep, n_neighbors=n_neighbors,
                    metric="cosine", random_state=seed)
    sc.tl.umap(peri, min_dist=0.9, spread=1.5, random_state=seed)
    sc.tl.leiden(peri, key_added=leiden_key, resolution=resolution)

    # Map cluster labels back to full AnnData
    adata.obs[leiden_key] = np.nan
    adata.obs.loc[peri.obs_names, leiden_key] = peri.obs[leiden_key].astype(str).values

    return adata


def analyze_pericytes_subclusters(
    adata: AnnData, outdir: Path, *, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", leiden_key: str = "leiden_pericytes",
    state_key: str = "pericyte_state", score_key: str = "airspace_dist_z",
    donor_key: str = "donor_id", min_cells_per_donor: int = 20,
) -> None:
    pc_dir = outdir / "per_cell"; pc_dir.mkdir(parents=True, exist_ok=True)
    pd_dir = outdir / "per_donor"; pd_dir.mkdir(parents=True, exist_ok=True)

    cols_needed = [cluster_key, leiden_key, state_key, score_key, donor_key,
                   "sex", "disease", "self_reported_ethnicity", "age_or_mean_of_age_range"]
    cols_existing = [col for col in cols_needed if col in adata.obs.columns]
    obs = adata.obs[cols_existing].copy()
    peri = obs[obs[cluster_key].astype(str) == pericyte_label].copy()
    if state_key in peri.columns and leiden_key in peri.columns:
        df = peri.loc[peri[state_key].notna() & peri[leiden_key].notna(), [state_key, leiden_key]].copy()
        if not df.empty and df[leiden_key].nunique() > 1 and df[state_key].nunique() > 1:
            # Contingency (row-normalized by Leiden)
            tab = pd.crosstab(df[leiden_key], df[state_key], normalize="index")
            tab.to_csv(pc_dir / "leiden_vs_state.csv")

            # ARI/NMI
            ari = adjusted_rand_score(df[state_key], df[leiden_key])
            nmi = normalized_mutual_info_score(df[state_key], df[leiden_key])
            with open(pc_dir / "leiden_vs_state_summary.txt", "w") as fh:
                fh.write(f"Adjusted Rand Index (ARI): {ari:.3f}\n")
                fh.write(f"Normalized Mutual Information (NMI): {nmi:.3f}\n\n")
                fh.write("Row-normalized contingency table (leiden x pericyte_state):\n")
                fh.write(tab.to_string())
                fh.write("\n")

            # Plot
            plt.figure(figsize=(max(6, 1.2 * tab.shape[0]), 4))
            tab.plot(kind="bar", stacked=True, ax=plt.gca())
            plt.ylabel("Proportion within Leiden cluster")
            plt.xlabel("Leiden cluster")
            plt.legend(title="Pericyte state", bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(pc_dir / "leiden_vs_state.png", dpi=300)
            plt.savefig(pc_dir / "leiden_vs_state.pdf")
            plt.close()

            # Chi-square + Cramer's V
            raw_tab = pd.crosstab(df[state_key], df[leiden_key])
            chi2, p, dof, _ = chi2_contingency(raw_tab)
            n = raw_tab.to_numpy().sum()
            r, k = raw_tab.shape
            cramer_v = np.sqrt(chi2 / (n * (min(r, k) - 1)))
            with open(pc_dir / "leiden_vs_state_chisq.txt", "w") as fh:
                fh.write(raw_tab.to_string() + "\n\n")
                fh.write(f"Chi^2 = {chi2:.4f}, df = {dof}, p = {p:.4e}\n")
                fh.write(f"Cramer's V = {cramer_v:.4f}\n")

    # Per cell: score ~ leiden
    if score_key in peri.columns and leiden_key in peri.columns:
        df = peri.loc[peri[score_key].notna() & peri[leiden_key].notna(), [score_key, leiden_key]].copy()
        if not df.empty and df[leiden_key].nunique() > 1:
            # Summary by Leiden
            stats = (
                df.groupby(leiden_key)[score_key]
                  .agg(["mean", "std", "count"])
                  .rename(columns={"count": "n_cells"})
                  .reset_index()
            )
            stats.to_csv(pc_dir / f"{score_key}_by_leiden.csv", index=False)

            # Violin
            plt.figure(figsize=(max(6, 1.2 * stats.shape[0]), 4))
            sns.violinplot(data=df, x=leiden_key, y=score_key, cut=0)
            plt.ylabel(score_key); plt.xlabel("Leiden cluster")
            plt.tight_layout()
            plt.savefig(pc_dir / f"{score_key}_by_leiden_violin.png", dpi=300)
            plt.savefig(pc_dir / f"{score_key}_by_leiden_violin.pdf")
            plt.close()

            # ANOVA + Kruskal
            groups = [vals[score_key].dropna().values for _, vals in df.groupby(leiden_key) if len(vals) > 1]
            anova_F, anova_p = f_oneway(*groups)
            kw_H, kw_p = kruskal(*groups)
            
            with open(pc_dir / f"{score_key}_leiden_stats.txt", "w") as fh:
                fh.write(f"ANOVA F = {anova_F:.4f}, p = {anova_p:.4e}\n")
                fh.write(f"Kruskal-Wallis H = {kw_H:.4f}, p = {kw_p:.4e}\n")

    # Per donor summaries
    if donor_key in peri.columns:
        # Filter tiny donors
        per_donor_counts = peri.groupby(donor_key, observed=False).size()
        keep_donors = per_donor_counts[per_donor_counts >= min_cells_per_donor].index
        dperi = peri.loc[peri[donor_key].isin(keep_donors)].copy()

        # (a) Donor-level state composition
        if state_key in dperi.columns:
            df = dperi.loc[dperi[state_key].notna(), [donor_key, state_key]]
            if not df.empty and df[state_key].nunique() > 1:
                tab = pd.crosstab(df[donor_key], df[state_key], normalize="index")
                tab.to_csv(pd_dir / "state_by_donor.csv")

        # (b) ARI/NMI per donor (state vs leiden)
        if all(k in dperi.columns for k in (state_key, leiden_key)):
            rows = []
            for donor, sub in dperi.groupby(donor_key, observed=False):
                sub = sub.dropna(subset=[state_key, leiden_key])
                if sub[state_key].nunique() > 1 and sub[leiden_key].nunique() > 1:
                    rows.append({
                        donor_key: donor,
                        "n_cells": int(sub.shape[0]),
                        "ARI": adjusted_rand_score(sub[state_key], sub[leiden_key]),
                        "NMI": normalized_mutual_info_score(sub[state_key], sub[leiden_key]),
                    })
            if rows:
                pd.DataFrame(rows).sort_values(donor_key).to_csv(
                    pd_dir / "leiden_state_metrics_by_donor.csv", index=False
                )

        # (c) Continuous score by donor X leiden + donor-aware OLS
        if score_key in dperi.columns and leiden_key in dperi.columns:
            df = dperi.dropna(subset=[score_key, leiden_key]).copy()
            if not df.empty and df[leiden_key].nunique() > 1:
                grp = (
                    df.groupby([donor_key, leiden_key], observed=False)[score_key]
                      .agg(["mean", "std", "count"])
                      .rename(columns={"count": "n_cells"})
                      .reset_index()
                )
                grp.to_csv(pd_dir / f"{score_key}_by_donor_leiden.csv", index=False)

                valid = df.groupby(donor_key, observed=False)[leiden_key].nunique() >= 2
                df2 = df[df[donor_key].isin(valid[valid].index)].copy()
                if not df2.empty and df2[leiden_key].nunique() > 1 and df2[donor_key].nunique() > 1:
                    formula = f"{score_key} ~ C({leiden_key}) + C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range"
                    fit = smf.ols(formula, data=df2).fit(cov_type="HC3") # unequal variances
                    anova = fit.compare_f_test(
                        smf.ols(f"{score_key} ~ C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range", data=df2).fit()
                    )
                    anova_table = anova_lm(fit, typ=2)
                    anova_table.to_csv(pd_dir / f"anova_{score_key}_ols.csv")
                    with open(pd_dir / f"donor_aware_ols_{score_key}.txt", "w") as fh:
                        fh.write(fit.summary().as_text())
                        fh.write("\n\nDonor-adjusted effect of Leiden:\n")
                        fh.write(f"F = {anova[0]:.4f}, p = {anova[1]:.4e}, df_diff = {anova[2]}\n")
                        
                    pairwise = fit.t_test_pairwise(f"C({leiden_key})")
                    summary = pairwise.result_frame
                    summary.to_csv(pd_dir / f"pairwise_leiden_{score_key}_ols.csv", index=False)


def analyze_pericytes_subclusters(
    df: pd.DataFrame, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", leiden_key: str = "leiden_pericytes",
    score_key: str = "airspace_score", imputed_key: str = "AGTR1_scvi",
    donor_key: str = "donor_id", min_cells_per_donor: int = 20,
) -> None:

    cols_needed = [cluster_key, leiden_key, score_key, imputed_key,
                   donor_key, "sex", "disease", "self_reported_ethnicity",
                   "age_or_mean_of_age_range"]
    cols_existing = [col for col in cols_needed if col in df.columns]
    peri = df[df[cluster_key].astype(str) == pericyte_label].copy()
    # Per donor summaries
    per_donor_counts = peri.groupby(donor_key, observed=False).size()
    keep_donors = per_donor_counts[per_donor_counts >= min_cells_per_donor].index
    dperi = peri.loc[peri[donor_key].isin(keep_donors)].copy()
    # Airspace vs imputed AGTR1
    df = dperi.dropna(subset=[score_key, leiden_key, imputed_key]).copy()
    formula = f"{score_key} ~ C({imputed_key}) + C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range"
    fit = smf.ols(formula, data=df).fit(cov_type="HC3") # unequal variances
    anova_table = anova_lm(fit, typ=2)
    anova_table.to_csv(outdir / f"anova_ols.{score_key}.{imputed_key}.tsv", sep="\t")
        
    # Airspace vs leident clusters
    valid = df.groupby(donor_key, observed=False)[leiden_key].nunique() >= 2
    df2 = df[df[donor_key].isin(valid[valid].index)].copy()
    if not df2.empty and df2[leiden_key].nunique() > 1 and df2[donor_key].nunique() > 1:
        formula = f"{score_key} ~ C({leiden_key}) + C({donor_key}) + C(sex) + C(disease) + C(self_reported_ethnicity) + age_or_mean_of_age_range"
        fit = smf.ols(formula, data=df2).fit(cov_type="HC3") # unequal variances
        anova_table = anova_lm(fit, typ=2)
        anova_table.to_csv(outdir / f"anova_{score_key}_ols.tsv", sep="\t")
        with open(outdir / f"donor_aware_ols_{score_key}.txt", "w") as fh:
            fh.write(fit.summary().as_text())
            fh.write("\n\nDonor-adjusted effect of Leiden:\n")
        pairwise = fit.t_test_pairwise(f"C({leiden_key})")
        summary  = pairwise.result_frame
        summary.to_csv(outdir / f"pairwise_leiden_{score_key}_ols.csv", index=False)


def main():
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)

    # Neighbor graph
    ensure_neighbors(adata, use_rep=args.use_rep,
                     n_neighbors=args.neighbors, seed=args.seed)

    # Graph distance to EC
    nn_dist = pericyte_to_nearest_ec_latent(
        adata, use_rep=args.use_rep, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, n_neighbors=1,
        ec_labels=tuple(args.ec_labels)
    )
    adata.obs["airspace_dist"] = nn_dist

    # Signature scores
    layer = "logcounts" if "logcounts" in adata.layers else None
    states = load_state_genes(args.state_yaml)
    score_signatures(adata, states, layer=layer)

    # Assign states
    donor_counts = assign_pericyte_states(
        adata, cluster_key=args.cluster_key, pericyte_label=args.pericyte_label,
        donor_key="donor_id", margin=0.4
    )

    # QC plots
    qc_dir = outdir / "qc_stats"
    qc_dir.mkdir(exist_ok=True, parents=True)

    plot_violin_by_state(
        adata, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, outdir=qc_dir,
    )
    batch_vars = ("donor_id", "data", "assay", "tissue_sampling_method",
                  "sequencing_platform", "development_stage", "tissue",
                  "subject_type", "study", "lung_condition", "sex",
                  "self_reported_ethnicity", "age_or_mean_of_age_range")
    batch_tables = plot_batch_composition(
        adata, outdir=qc_dir, state_key="pericyte_state",
        batch_keys=batch_vars, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label,
    )

    # De novo clustering on pericytes
    adata = run_leiden_on_pericytes(
        adata, use_rep=args.use_rep, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, resolution=0.5,
        leiden_key="leiden_pericytes", n_neighbors=args.neighbors, seed=args.seed
    )

    # Per cell and per donor analysis
    summary = analyze_pericytes_subclusters(
        adata, outdir=outdir / "pericyte_leiden_analysis",
    )
    print(summary)

    # Imputed denoising per donor
    impute_dir = outdir / "airspace"
    imputed_path = impute_dir / "pericytes_airspace_denoising.tsv"
    imputed_df = pd.read_csv(imputed_path, sep="\t", index_col=0)
    obs = adata.obs.copy()
    df  = obs.merge(imputed_df, left_index=True, right_index=True)
    analyze_pericytes_subclusters(df, impute_dir)
    
    # Save AnnData with pericyte_state annotations
    adata.write(outdir / "airspace_pericyte_states.h5ad")
    logging.info("Pericyte state analysis complete.")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
