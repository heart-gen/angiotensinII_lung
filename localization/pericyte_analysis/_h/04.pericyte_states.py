"""
Pericyte state classification using airspace proximity and marker signatures.

Description:
- Compute graph distance from each pericyte to alveolar capillary EC clusters
  ("EC aerocyte capillary", "EC general capillary") in latent space.
- Compute signature scores for:
    * capillary-like: KCNJ8, ABCC9
    * arteriolar-like: NOTCH3, THY1
- Combine distances + signatures into a single score to define:
    adata.obs['pericyte_state'] ∈ {'capillary_like', 'arteriolar_like'}

Outputs:
- adata.obs['pericyte_state']
- adata.obs['airspace_dist_z'] (z-scored graph distance for pericytes)
- adata.obs['capillary_sig_z'], adata.obs['arteriolar_sig_z'],
  adata.obs['sig_diff_z'] (z-scored signature scores & difference)
- QC markdown: pericyte_state_qc.md
"""

import argparse
import logging
from pathlib import Path
from typing import Sequence, Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
from scipy import sparse
from anndata import AnnData
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import chi2_contingency
from scipy.stats import f_oneway, kruskal
from scipy.sparse.csgraph import shortest_path
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
    parser.add_argument(
        "--adata", type=Path, required=True,
        help="Input AnnData (e.g., airspace.hlca_full.dataset.h5ad)"
    )
    parser.add_argument(
        "--outdir", type=Path, required=True,
        help="Output directory for QC and updated AnnData"
    )
    parser.add_argument(
        "--use-rep", type=str, default="X_pca_harmony",
        help="obsm key for latent representation (default: X_pca_harmony)"
    )
    parser.add_argument(
        "--cluster-key", type=str, default="subclusters",
        help="obs column with fine cell types (default: subclusters)"
    )
    parser.add_argument(
        "--pericyte-label", type=str, default="Pericytes",
        help="Value in cluster-key used for pericytes"
    )
    parser.add_argument(
        "--ec-labels", nargs="+",
        default=["EC aerocyte capillary", "EC general capillary"],
        help="Labels treated as alveolar capillary EC"
    )
    parser.add_argument(
        "--neighbors", type=int, default=30,
        help="n_neighbors for graph construction if neighbors not present"
    )
    parser.add_argument(
        "--seed", type=int, default=13, help="Random seed for neighbors"
    )
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
        metric="euclidean", random_state=seed,
    )


def compute_graph_distance_to_ec(
    adata: AnnData, cluster_key: str, pericyte_label: str, ec_labels: Sequence[str],
) -> np.ndarray:
    """
    Compute shortest-path distance from each pericyte to the nearest EC
    (from specified ec_labels) using the neighbors graph.
    """
    if "distances" not in adata.obsp:
        raise ValueError("adata.obsp['distances'] missing; run sc.pp.neighbors first.")

    dist_graph = adata.obsp["distances"]
    if not issparse(dist_graph):
        dist_graph = sparse.csr_matrix(dist_graph)

    dist_graph = dist_graph.maximum(dist_graph.T)

    n_cells = adata.n_obs
    cluster = adata.obs[cluster_key].astype(str)

    pericyte_mask = cluster == pericyte_label
    ec_mask = cluster.isin(ec_labels)

    if ec_mask.sum() == 0:
        raise ValueError(f"No EC cells found for labels={ec_labels}")

    pericyte_idx = np.where(pericyte_mask)[0]
    ec_idx = np.where(ec_mask)[0]

    logging.info(
        "Computing shortest paths from %d EC cells to all %d cells.",
        ec_idx.size, n_cells,
    )

    # distances from each EC cell to all cells
    dist_from_ec = shortest_path(
        dist_graph, directed=False, unweighted=False,
        indices=ec_idx, overwrite=False,
    )

    # min distance from any EC to each cell
    min_to_ec = np.min(dist_from_ec, axis=0)
    pericyte_to_ec = np.full(n_cells, np.nan, dtype=float)
    pericyte_to_ec[pericyte_idx] = min_to_ec[pericyte_idx]
    return pericyte_to_ec


def score_signatures(
    adata: AnnData, capillary_genes: Sequence[str] = ("KCNJ8", "ABCC9"),
    arteriolar_genes: Sequence[str] = ("NOTCH3", "RGS5", "ACTA2"),
    layer: str | None = "logcounts",
) -> None:
    """
    Compute simple module scores for capillary vs arteriolar signatures.
    Uses scanpy.tl.score_genes.
    """
    # Filter gene lists to those present
    cap_genes_present = [g for g in capillary_genes if g in adata.var_names]
    art_genes_present = [g for g in arteriolar_genes if g in adata.var_names]

    logging.info("Capillary signature genes present: %s", cap_genes_present)
    logging.info("Arteriolar signature genes present: %s", art_genes_present)

    if not cap_genes_present:
        logging.warning("No capillary genes found in var_names!")
    if not art_genes_present:
        logging.warning("No arteriolar genes found in var_names!")

    sc.tl.score_genes(
        adata, gene_list=cap_genes_present, score_name="capillary_sig",
        layer=layer, use_raw=False,
    )
    sc.tl.score_genes(
        adata, gene_list=art_genes_present, score_name="arteriolar_sig",
        layer=layer, use_raw=False,
    )


def zscore_series(x: pd.Series) -> pd.Series:
    """Simple z-score with NaN handling."""
    m = x.mean(skipna=True)
    s = x.std(skipna=True)
    if s == 0 or np.isnan(s):
        return pd.Series(np.zeros_like(x), index=x.index)
    return (x - m) / s


def assign_pericyte_states(
    adata: AnnData, cluster_key: str, pericyte_label: str,
    min_donors_per_state: int = 5, donor_key: str = "donor_id",
) -> Dict[str, int]:
    """
    Combine z-scored distance + signature scores into a pericyte_state:
        state = 'capillary_like' if combined_score >= 0 else 'arteriolar_like'
    combined_score = -airspace_dist_z + sig_diff_z

    Returns a dict with donor counts per state.
    """
    obs = adata.obs
    pericyte_mask = obs[cluster_key].astype(str) == pericyte_label

    # Z-score distance & signatures within pericytes only
    dist_z = pd.Series(np.nan, index=obs.index)
    sig_cap_z = pd.Series(np.nan, index=obs.index)
    sig_art_z = pd.Series(np.nan, index=obs.index)

    dist_z.loc[pericyte_mask] = zscore_series(obs.loc[pericyte_mask, "airspace_dist"])
    sig_cap_z.loc[pericyte_mask] = zscore_series(obs.loc[pericyte_mask, "capillary_sig"])
    sig_art_z.loc[pericyte_mask] = zscore_series(obs.loc[pericyte_mask, "arteriolar_sig"])

    obs["airspace_dist_z"] = dist_z
    obs["capillary_sig_z"] = sig_cap_z
    obs["arteriolar_sig_z"] = sig_art_z
    obs["sig_diff_z"] = sig_cap_z - sig_art_z

    # Combined score: higher → more capillary-like
    combined = -dist_z + obs["sig_diff_z"]
    obs["pericyte_combined_score"] = combined

    # Assign states only within pericytes
    state = pd.Series(np.nan, index=obs.index, dtype=object)
    state.loc[pericyte_mask & combined.notna()] = np.where(
        combined.loc[pericyte_mask & combined.notna()] >= 0,
        "capillary_like", "arteriolar_like",
    )
    obs["pericyte_state"] = state

    # Donor balance
    state_donor_counts: Dict[str, int] = {}
    if donor_key in obs.columns:
        for st in ["capillary_like", "arteriolar_like"]:
            mask = (obs["pericyte_state"] == st)
            donors = obs.loc[mask & pericyte_mask, donor_key].dropna().unique()
            state_donor_counts[st] = len(donors)
            logging.info("State %s has %d donors", st, len(donors))

        for st, nd in state_donor_counts.items():
            if nd < min_donors_per_state:
                logging.warning(
                    "State %s has only %d donors (< %d); donor balance criterion FAILED.",
                    st, nd, min_donors_per_state
                )
    else:
        logging.warning("No donor_key '%s' in obs; cannot assess donor balance.", donor_key)

    return state_donor_counts


def run_leiden_on_pericytes(
    adata: AnnData, use_rep: str, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", resolution: float = 0.5,
    neighbors_key: str | None = None, leiden_key: str = "leiden_pericytes",
) -> AnnData:
    """
    Run Leiden clustering on pericytes only, using given representation.
    """
    per = adata[adata.obs[cluster_key] == pericyte_label].copy()
    sc.pp.neighbors(
        per, use_rep=use_rep, n_neighbors=75, random_state=13,
    )
    sc.tl.leiden(per, key_added=leiden_key, resolution=resolution)

    # Map cluster labels back to full AnnData
    adata.obs[leiden_key] = np.nan
    adata.obs.loc[per.obs_names, leiden_key] = per.obs[leiden_key].astype(str).values

    return adata


def analyze_leiden_vs_pericyte_state(
    adata: AnnData, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", state_key: str = "pericyte_state",
    leiden_key: str = "leiden_pericytes",
) -> None:
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()

    # Drop cells without state or cluster
    mask = per[state_key].notna() & per[leiden_key].notna()
    per = per.loc[mask].copy()

    # Contingency table
    tab = pd.crosstab(per[leiden_key], per[state_key], normalize="index")
    tab.to_csv(outdir / "leiden_vs_pericyte_state.csv")

    # Metrics of alignment
    ari = adjusted_rand_score(per[state_key], per[leiden_key])
    nmi = normalized_mutual_info_score(per[state_key], per[leiden_key])

    with open(outdir / "leiden_vs_pericyte_state_summary.txt", "w") as fh:
        fh.write(f"Adjusted Rand Index (ARI): {ari:.3f}\n")
        fh.write(f"Normalized Mutual Information (NMI): {nmi:.3f}\n\n")
        fh.write("Row-normalized contingency table (leiden x pericyte_state):\n")
        fh.write(tab.to_string())
        fh.write("\n")

    # Quick plot
    plt.figure(figsize=(max(6, 1.2 * tab.shape[0]), 4))
    tab.plot(kind="bar", stacked=True, ax=plt.gca())
    plt.ylabel("Proportion within Leiden cluster")
    plt.xlabel("Leiden cluster")
    plt.legend(title="Pericyte state", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(outdir / "leiden_vs_pericyte_state.png", dpi=300)
    plt.savefig(outdir / "leiden_vs_pericyte_state.pdf")
    plt.close()


def analyze_leiden_vs_airspace(
    adata: AnnData, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", leiden_key: str = "leiden_pericytes",
    score_key: str = "airspace_dist_z",
) -> None:
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    per = per[per[leiden_key].notna() & per[score_key].notna()].copy()

    # Cluster-level stats
    stats = (
        per.groupby(leiden_key)[score_key]
        .agg(["mean", "std", "count"])
        .rename(columns={"count": "n_cells"})
        .reset_index()
    )
    stats.to_csv(outdir / f"{score_key}_by_leiden.csv", index=False)

    # Violin plot
    plt.figure(figsize=(max(6, 1.2 * stats.shape[0]), 4))
    sns.violinplot(data=per, x=leiden_key, y=score_key, cut=0)
    plt.ylabel(score_key)
    plt.xlabel("Leiden cluster")
    plt.tight_layout()
    plt.savefig(outdir / f"{score_key}_by_leiden_violin.png", dpi=300)
    plt.savefig(outdir / f"{score_key}_by_leiden_violin.pdf")
    plt.close()


def statistical_tests_continuous_vs_leiden(
    adata: AnnData, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", leiden_key: str = "leiden_pericytes",
    score_key: str = "airspace_dist_z",
):
    """
    Test whether the continuous score differs among Leiden clusters.
    """
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    per = per[per[leiden_key].notna() & per[score_key].notna()].copy()

    groups = {
        cl: vals[score_key].dropna().values
        for cl, vals in per.groupby(leiden_key)
        if len(vals[score_key].dropna()) > 1
    }

    # Only compute tests if >= 2 groups with >1 cell
    if len(groups) < 2:
        print("Not enough valid Leiden clusters for statistical tests.")
        return

    # ANOVA
    try:
        anova_F, anova_p = f_oneway(*groups.values())
    except Exception:
        anova_F, anova_p = np.nan, np.nan

    # Kruskal-Wallis
    try:
        kw_H, kw_p = kruskal(*groups.values())
    except Exception:
        kw_H, kw_p = np.nan, np.nan

    # Save results
    out = outdir / f"{score_key}_leiden_stats.txt"
    with open(out, "w") as fh:
        fh.write(f"Statistical tests for {score_key} vs Leiden clusters\n")
        fh.write("------------------------------------------------------\n\n")
        fh.write("Groups included:\n")
        for cl, arr in groups.items():
            fh.write(f"  - Leiden {cl}: N={len(arr)}\n")
        fh.write("\n")

        fh.write(f"ANOVA F = {anova_F:.4f}, p = {anova_p:.4e}\n")
        fh.write(f"Kruskal-Wallis H = {kw_H:.4f}, p = {kw_p:.4e}\n")

    print(f"Saved statistical test results to: {out}")


def statistical_tests_state_vs_leiden(
    adata: AnnData, outdir: Path, cluster_key: str = "subclusters",
    pericyte_label: str = "Pericytes", state_key: str = "pericyte_state",
    leiden_key: str = "leiden_pericytes",
):
    """Chi-square test and Cramer's V"""
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    per = per[per[state_key].notna() & per[leiden_key].notna()].copy()

    tab = pd.crosstab(per[state_key], per[leiden_key])
    chi2, p, dof, expected = chi2_contingency(tab)

    # Cramer's V
    n = tab.sum().sum()
    r, k = tab.shape
    cramers_v = np.sqrt(chi2 / (n * (min(r, k) - 1)))

    out = outdir / "pericyte_state_vs_leiden_stats.txt"
    with open(out, "w") as fh:
        fh.write("Chi-square association test: pericyte_state ~ leiden\n")
        fh.write("------------------------------------------------------\n")
        fh.write(tab.to_string())
        fh.write("\n\n")
        fh.write(f"Chi^2 = {chi2:.4f}, df = {dof}, p = {p:.4e}\n")
        fh.write(f"Cramer's V = {cramers_v:.4f}\n")

    print(f"Saved categorical association results to: {out}")


def plot_violin_by_state(
    adata: AnnData, cluster_key: str, pericyte_label: str, outdir: Path,
) -> None:
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()

    def _violin(col: str, fname: str, ylabel: str):
        plt.figure(figsize=(5, 5))
        sns.violinplot(
            data=per, x="pericyte_state", y=col,
            order=["capillary_like", "arteriolar_like"], cut=0,
        )
        sns.stripplot(
            data=per, x="pericyte_state", y=col,
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
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    results: Dict[str, pd.DataFrame] = {}

    for key in batch_keys:
        if key not in per.columns:
            continue
        tab = pd.crosstab(per[key], per[state_key], normalize="index")
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


def write_qc_report(
    outdir: Path, adata: AnnData, cluster_key: str, pericyte_label: str,
    donor_counts: Dict[str, int], batch_tables: Dict[str, pd.DataFrame],
) -> None:
    per = adata.obs[adata.obs[cluster_key] == pericyte_label].copy()
    n_per = per.shape[0]
    state_counts = per["pericyte_state"].value_counts(dropna=False)

    lines: List[str] = []
    lines.append("# Pericyte state QC")
    lines.append("")
    lines.append(f"- Total pericytes: **{n_per}**")
    lines.append("- Pericyte states (cell counts):")
    for st, cnt in state_counts.items():
        lines.append(f"  - `{st}`: {cnt}")
    lines.append("")
    lines.append("## Donor balance")
    if donor_counts:
        for st, nd in donor_counts.items():
            status = "OK" if nd >= 5 else "FAIL"
            lines.append(f"- State `{st}`: {nd} donors (**{status}**)")
    else:
        lines.append("- Donor information not available; donor balance not assessed.")

    lines.append("")
    lines.append("## Batch / dataset composition")
    if not batch_tables:
        lines.append("No batch or dataset covariates were available for QC.")
    else:
        for key, tab in batch_tables.items():
            lines.append(f"### {key}")
            lines.append("")
            lines.append("State composition by " + key + " (rows sum to 1):")
            lines.append("")
            lines.append(tab.to_markdown())
            lines.append("")

    lines.append("")
    lines.append("## Notes")
    lines.append("- Airspace distance (`airspace_dist_z`) is z-scored among pericytes.")
    lines.append("- Signature difference (`sig_diff_z`) = capillary_sig_z - arteriolar_sig_z.")
    lines.append("- Combined score = -airspace_dist_z + sig_diff_z, thresholded at 0 for state assignment.")

    (outdir / "pericyte_state_qc.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    configure_logging()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)

    # Neighbor graph
    ensure_neighbors(adata, use_rep=args.use_rep,
                     n_neighbors=args.neighbors, seed=args.seed)

    # Graph distance to EC
    dist_to_ec = compute_graph_distance_to_ec(
        adata, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label,
        ec_labels=args.ec_labels,
    )
    adata.obs["airspace_dist"] = dist_to_ec

    # Signature scores
    layer = "logcounts" if "logcounts" in adata.layers else None
    score_signatures(adata, layer=layer)

    # Assign states
    donor_counts = assign_pericyte_states(
        adata, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label,
        min_donors_per_state=5, donor_key="donor_id",
    )

    # QC plots and report
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

    write_qc_report(
        outdir=qc_dir, adata=adata, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, donor_counts=donor_counts,
        batch_tables=batch_tables,
    )

    # De novo clustering on pericytes
    adata = run_leiden_on_pericytes(
        adata, use_rep=args.use_rep,
        cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label,
        resolution=1.0, leiden_key="leiden_pericytes",
    )

    # Compare pericyte_state vs Leiden
    analyze_leiden_vs_pericyte_state(
        adata, outdir=qc_dir, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, state_key="pericyte_state",
        leiden_key="leiden_pericytes",
    )

    # Compare airspace distance / score vs Leiden
    analyze_leiden_vs_airspace(
        adata, outdir=qc_dir, cluster_key=args.cluster_key,
        pericyte_label=args.pericyte_label, leiden_key="leiden_pericytes",
        score_key="airspace_dist_z",
    )

    # Save AnnData with pericyte_state annotations
    adata.write(outdir / "airspace_pericyte_states.h5ad")
    logging.info("Pericyte state analysis complete.")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
