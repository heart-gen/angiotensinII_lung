"""
De novo pericyte state discovery (NVU pattern).

States are STABLE, data-driven Leiden clusters on the study-integrated embedding
(X_pca_harmony) selected by a bootstrap stability sweep -- they are NOT defined by
marker-score argmax. Curated marker panels and Wilcoxon DE then ANNOTATE the
clusters (characterize, not define), and per-state donor/source support is
reported so no state rests on a single donor or study. Because the embedding is
already study-integrated, study is accounted for exactly once (at clustering);
module scores are taken at face value for annotation.

Mirrors nvu-neuroimmune-reference/03_state_discovery/{01_cluster_stability,
02_annotate_states}.py:
  1. stability sweep: neighbors x resolutions x bootstraps -> per-cluster Jaccard
     + ARI/NMI; keep solutions passing median-Jaccard / min-cells / min-donors;
     choose the most stable.
  2. annotate: score 6 curated programs (mean per cluster), Wilcoxon markers,
     per-cluster dominant program -> state_program label.

The basement_membrane panel was added after basement_membrane/_h/01.state_gate.py
showed that clusters 1/3/5 (4,200 cells, 36% of pericytes) are dominantly enriched
for basement-membrane deposition rather than fibrillar ECM. Under the previous
5-panel set those clusters were labelled `fibroblast_like` by an argmax over
programs that were ALL negatively enriched -- a least-negative default, not a
positive identity. Clustering is unaffected (it runs on X_pca_harmony, which does
not depend on the panels); only the annotation changes.

Outputs:
  - pericyte_states.h5ad             (obs: pericyte_state [stable cluster], state_program)
  - pericytes_states_metadata.tsv.gz (per-cell obs for donor-aware R stats)
  - stability/cluster_stability_grid.tsv, chosen_cluster_solution.tsv,
    cluster_size_donor_grid.tsv
  - annotations/state_markers.tsv.gz, state_summary.tsv, state_program_map.tsv
  - figures/umap_{pericyte_state,state_program,<score>}.{png,pdf}
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from scipy import sparse
from anndata import AnnData
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sns.set_context("talk")
sns.set_style("whitegrid")

# Curated human marker panels per functional program. Used to ANNOTATE the
# data-driven clusters, never to define them. Genes absent from the object are
# dropped at runtime with a warning.
STATE_PANELS = {
    "vascular_stabilizing": [
        "RGS5", "PDGFRB", "NOTCH3", "KCNJ8", "ABCC9", "HIGD1B",
        "COX4I2", "NDUFA4L2", "GJA4", "CSPG4", "PLXDC1",
    ],
    "inflammatory": [
        "IL6", "CCL2", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CXCL6",
        "CXCL8", "CXCL10", "IL1A", "IL1B", "MIF", "ICAM1", "VCAM1",
        "SELE", "NFKBIA",
    ],
    "synthetic_contractile": [
        "ACTA2", "MYH11", "TAGLN", "CNN1", "MYL9", "DES", "VIM", "PDGFRB",
    ],
    "activated_migratory": [
        "ADAMTS1", "THBS1", "TIMP1", "MMP2", "MMP3", "MMP9", "SERPINE1",
        "POSTN",
    ],
    "fibroblast_like": [
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "FN1", "LUM", "DCN",
        "PDGFA", "FBLN1",
    ],
    # Basement-membrane deposition. Kept SEPARATE from fibroblast_like because
    # the two are biologically distinct matrices and, empirically, near-orthogonal
    # in these cells (panel-score r = 0.05 once the shared COL4A1 is removed).
    # Panel definition and provenance: basement_membrane/_h/bm_panels.py.
    "basement_membrane": [
        "COL4A1", "COL4A2", "COL18A1", "LAMA3", "LAMA4", "LAMA5",
        "LAMB1", "LAMB2", "LAMC1", "NID1", "NID2", "HSPG2", "AGRN",
    ],
}

EXTRA_GENES = ["AGTR1", "AGTR2", "ACTA2"]


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="Pericyte AnnData with integrated embedding")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--use-rep", default="X_pca_harmony",
                   help="Study-integrated embedding to cluster on")
    p.add_argument("--neighbors", default="15,30",
                   help="Comma-separated n_neighbors settings")
    p.add_argument("--resolutions", default="0.3,0.5,0.7,0.9",
                   help="Comma-separated Leiden resolutions")
    p.add_argument("--n-bootstrap", type=int, default=30)
    p.add_argument("--resample-fraction", type=float, default=0.8)
    p.add_argument("--stability-threshold", type=float, default=0.6,
                   help="Min median best-Jaccard for a solution to 'pass'")
    p.add_argument("--min-cluster-cells", type=int, default=50)
    p.add_argument("--min-cluster-donors", type=int, default=3)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def save_figure(fig, base_path: Path):
    fig.savefig(base_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def load_anndata(path: Path) -> AnnData:
    adata = sc.read_h5ad(path)
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        if not np.array_equal(symbols.to_numpy(), adata.var_names.to_numpy()):
            adata.var["ensembl_id"] = adata.var_names
            adata.var_names = adata.var.index = symbols
    adata.raw = None
    adata.var_names_make_unique()
    return adata


def add_gene_expression(adata: AnnData, gene: str):
    if gene not in adata.var_names:
        logging.warning(f"{gene} not in var_names")
        return
    expr = adata[:, gene].layers["logcounts"]
    expr = expr.toarray().ravel() if sparse.issparse(expr) else np.asarray(expr).ravel()
    adata.obs[f"{gene}_expr"] = expr
    adata.obs[f"{gene}_detect"] = (expr > 0).astype(int)


def add_depth(adata: AnnData):
    """Per-cell depth for the downstream total-counts partial correlation."""
    for key in ("total_counts", "n_counts"):
        if key in adata.obs.columns:
            adata.obs["total_counts"] = np.asarray(adata.obs[key], dtype=float)
            break
    else:
        adata.obs["total_counts"] = np.asarray(
            adata.layers["counts"].sum(axis=1)).ravel().astype(float)
    for key in ("n_genes", "n_genes_by_counts"):
        if key in adata.obs.columns:
            adata.obs["n_genes"] = np.asarray(adata.obs[key], dtype=float)
            break


# ---- stability sweep (NVU pattern) --------------------------------------
def cluster_embedding(emb: np.ndarray, n_neighbors: int, resolution: float,
                      seed: int) -> np.ndarray:
    """Leiden on a minimal AnnData carrying only the integrated embedding."""
    tmp = AnnData(np.zeros((emb.shape[0], 1), dtype="float32"))
    tmp.obsm["rep"] = emb
    sc.pp.neighbors(tmp, use_rep="rep", n_neighbors=min(n_neighbors, emb.shape[0] - 1),
                    metric="cosine", random_state=seed)
    sc.tl.leiden(tmp, resolution=resolution, random_state=seed,
                 key_added="leiden", flavor="leidenalg")
    return tmp.obs["leiden"].astype(str).to_numpy()


def best_jaccard(original: np.ndarray, bootstrap: np.ndarray) -> list:
    """Best bootstrap overlap (Jaccard) for each original cluster."""
    out = []
    boot_sets = {b: set(np.flatnonzero(bootstrap == b).tolist())
                 for b in np.unique(bootstrap)}
    for cluster in np.unique(original):
        a = set(np.flatnonzero(original == cluster).tolist())
        best = 0.0
        for b in boot_sets.values():
            denom = len(a | b)
            if denom:
                best = max(best, len(a & b) / denom)
        out.append(best)
    return out


def stability_sweep(adata, rep, neighbors, resolutions, n_bootstrap,
                    resample_fraction, thr, min_cells, min_donors, seed, outdir):
    emb = np.asarray(adata.obsm[rep])
    donors = adata.obs["donor_id"].to_numpy()
    rng = np.random.default_rng(seed)
    records, cluster_records, labels_by_key = [], [], {}

    for nn in neighbors:
        for res in resolutions:
            key = f"leiden_n{nn}_r{res}"
            logging.info(f"Clustering {key}")
            original = cluster_embedding(emb, nn, res, seed)
            labels_by_key[key] = original
            clusters = np.unique(original)

            jacc, aris, nmis = [], [], []
            for i in range(n_bootstrap):
                keep = np.sort(rng.choice(emb.shape[0],
                               size=int(emb.shape[0] * resample_fraction),
                               replace=False))
                boot = cluster_embedding(emb[keep], nn, res, seed + i + 1)
                orig_sub = original[keep]
                jacc.extend(best_jaccard(orig_sub, boot))
                aris.append(adjusted_rand_score(orig_sub, boot))
                nmis.append(normalized_mutual_info_score(orig_sub, boot))

            sizes = pd.Series(original).value_counts()
            dcounts = (pd.DataFrame({"c": original, "d": donors})
                       .groupby("c")["d"].nunique())
            median_j = float(np.median(jacc)) if jacc else 0.0
            records.append({
                "n_neighbors": nn, "resolution": res, "cluster_key": key,
                "n_clusters": len(clusters), "median_jaccard": median_j,
                "mean_ari": float(np.mean(aris)), "mean_nmi": float(np.mean(nmis)),
                "min_cluster_cells": int(sizes.min()),
                "min_cluster_donors": int(dcounts.min()),
                "passes": bool(median_j >= thr and sizes.min() >= min_cells
                               and dcounts.min() >= min_donors),
            })
            for c, n in sizes.items():
                cluster_records.append({"cluster_key": key, "cluster": c,
                                        "n_cells": int(n),
                                        "n_donors": int(dcounts.loc[c])})

    summary = (pd.DataFrame(records)
               .sort_values(["passes", "median_jaccard", "n_clusters"],
                            ascending=[False, False, False]))
    sdir = outdir / "stability"; sdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(sdir / "cluster_stability_grid.tsv", sep="\t", index=False)
    pd.DataFrame(cluster_records).to_csv(
        sdir / "cluster_size_donor_grid.tsv", sep="\t", index=False)

    chosen = summary.iloc[0]
    if not chosen["passes"]:
        logging.warning("No solution passed all gates; taking the most stable.")
    chosen_key = str(chosen["cluster_key"])
    pd.DataFrame([chosen.to_dict()]).to_csv(
        sdir / "chosen_cluster_solution.tsv", sep="\t", index=False)
    logging.info(f"Chosen {chosen_key}: {chosen['n_clusters']} clusters, "
                 f"median Jaccard {chosen['median_jaccard']:.2f}, "
                 f"passes={chosen['passes']}")
    # Relabel clusters by descending size for stable, readable ids.
    labels = labels_by_key[chosen_key]
    order = pd.Series(labels).value_counts().index.tolist()
    remap = {old: str(i) for i, old in enumerate(order)}
    return pd.Categorical([remap[x] for x in labels],
                          categories=[str(i) for i in range(len(order))])


# ---- annotation (NVU pattern: markers characterize the clusters) --------
def score_panels(adata, seed):
    adata.X = adata.layers["logcounts"]
    score_cols = []
    for state, genes in STATE_PANELS.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            logging.warning(f"[{state}] dropping absent genes: {missing}")
        if not present:
            logging.warning(f"[{state}] no genes present; skipping")
            continue
        col = f"{state}_score"
        sc.tl.score_genes(adata, present, score_name=col, random_state=seed,
                          use_raw=False)
        score_cols.append(col)
    return score_cols


def annotate_states(adata, score_cols, outdir):
    adir = outdir / "annotations"; adir.mkdir(parents=True, exist_ok=True)

    # Wilcoxon DE markers per stable cluster (markers characterize the cluster).
    sc.tl.rank_genes_groups(adata, groupby="pericyte_state", method="wilcoxon",
                            layer="logcounts", use_raw=False,
                            key_added="rank_genes_state")
    res = adata.uns["rank_genes_state"]
    rows = []
    for grp in res["names"].dtype.names:
        for gene, sc_, p, lfc in zip(res["names"][grp][:50], res["scores"][grp][:50],
                                     res["pvals_adj"][grp][:50],
                                     res["logfoldchanges"][grp][:50]):
            rows.append({"pericyte_state": grp, "gene": gene, "score": sc_,
                         "pval_adj": p, "logfoldchange": lfc})
    pd.DataFrame(rows).to_csv(adir / "state_markers.tsv.gz", sep="\t", index=False)

    # Per-cluster program annotation by RELATIVE enrichment, not raw magnitude.
    # Raw mean scores are not comparable across programs: contractile markers
    # (ACTA2/PDGFRB/VIM) are uniformly high in pericytes, so a raw-magnitude argmax
    # collapses every cluster to "synthetic_contractile". Following the NVU pattern
    # (02_annotate_states.py keeps scores descriptive and lets markers/relative
    # enrichment drive labels), we z-score each program across cells so programs
    # share a scale, then label each cluster by the program it is MOST enriched for
    # relative to the average pericyte. Raw means are retained in the summary.
    programs = [c.replace("_score", "") for c in score_cols]
    mean_scores = adata.obs.groupby("pericyte_state", observed=True)[score_cols].mean()
    zscored = adata.obs[score_cols].apply(
        lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9))
    zscored["pericyte_state"] = adata.obs["pericyte_state"].to_numpy()
    rel_enrich = zscored.groupby("pericyte_state", observed=True)[score_cols].mean()
    rel_enrich.columns = [f"{p}_relenrich" for p in programs]
    dominant = rel_enrich.to_numpy().argmax(axis=1)
    prog_map = pd.Series([programs[i] for i in dominant], index=rel_enrich.index,
                         name="state_program")
    pd.concat([prog_map, rel_enrich], axis=1).reset_index().to_csv(
        adir / "state_program_map.tsv", sep="\t", index=False)
    adata.obs["state_program"] = pd.Categorical(
        adata.obs["pericyte_state"].map(prog_map).astype(str))

    # Per-state summary with donor / source support (the rigor check).
    summary = (adata.obs.groupby("pericyte_state", observed=True)
               .agg(n_cells=("pericyte_state", "size"),
                    n_donors=("donor_id", "nunique"),
                    n_studies=("study", "nunique"),
                    n_datasets=("dataset", "nunique"),
                    AGTR1_pos_frac=("AGTR1_detect", "mean"),
                    AGTR1_mean=("AGTR1_expr", "mean"))
               .reset_index())
    summary = summary.merge(mean_scores.reset_index(), on="pericyte_state", how="left")
    summary = summary.merge(rel_enrich.reset_index(), on="pericyte_state", how="left")
    summary = summary.merge(prog_map.reset_index(), on="pericyte_state", how="left")
    summary.to_csv(adir / "state_summary.tsv", sep="\t", index=False)
    logging.info("State summary (donor / source support):\n" +
                 summary[["pericyte_state", "n_cells", "n_donors", "n_studies",
                          "state_program"]].to_string(index=False))
    return prog_map


def plot_umap(adata, col, outdir: Path, categorical=False):
    coords = adata.obsm["X_umap"]
    vals = adata.obs[col]
    fig, ax = plt.subplots(figsize=(6, 5))
    if categorical:
        cats = pd.Categorical(vals)
        palette = sns.color_palette("tab10", len(cats.categories))
        for color, cat_val in zip(palette, cats.categories):
            m = cats == cat_val
            ax.scatter(coords[m, 0], coords[m, 1], s=6, linewidths=0,
                       color=color, label=str(cat_val))
        ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize="small")
    else:
        sca = ax.scatter(coords[:, 0], coords[:, 1], c=vals, s=6, linewidths=0,
                         cmap="viridis")
        fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(col)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    save_figure(fig, outdir / f"umap_{col}")


def main():
    args = parse_args()
    configure_logging()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    adata = load_anndata(args.adata)
    logging.info(f"Loaded {adata.shape[0]} pericytes x {adata.shape[1]} genes")
    if args.use_rep not in adata.obsm:
        raise ValueError(f"{args.use_rep} not in obsm: {list(adata.obsm)}")

    for gene in EXTRA_GENES:
        add_gene_expression(adata, gene)
    add_depth(adata)

    neighbors = [int(x) for x in args.neighbors.split(",") if x.strip()]
    resolutions = [float(x) for x in args.resolutions.split(",") if x.strip()]
    adata.obs["pericyte_state"] = stability_sweep(
        adata, args.use_rep, neighbors, resolutions, args.n_bootstrap,
        args.resample_fraction, args.stability_threshold,
        args.min_cluster_cells, args.min_cluster_donors, args.seed, outdir)

    score_cols = score_panels(adata, seed=args.seed)
    annotate_states(adata, score_cols, outdir)

    # Plots
    fig_dir = outdir / "figures"; fig_dir.mkdir(exist_ok=True)
    plot_umap(adata, "pericyte_state", fig_dir, categorical=True)
    plot_umap(adata, "state_program", fig_dir, categorical=True)
    for col in score_cols:
        plot_umap(adata, col, fig_dir, categorical=False)
    if "leiden_pericytes" in adata.obs:
        plot_umap(adata, "leiden_pericytes", fig_dir, categorical=True)

    # Crosswalk: stable state x prior Leiden (transparency)
    if "leiden_pericytes" in adata.obs:
        cross = pd.crosstab(adata.obs["leiden_pericytes"], adata.obs["pericyte_state"])
        cross.to_csv(outdir / "leiden_state_crosswalk.tsv", sep="\t")

    adata.obs["pericyte_state"].value_counts().to_csv(
        outdir / "state_counts.tsv", sep="\t")

    # Per-cell metadata for donor-aware R stats
    keep_cols = [c for c in [
        "pericyte_state", "state_program", "leiden_pericytes", "donor_id",
        "disease", "lung_condition", "subject_type", "smoking_status", "sex",
        "self_reported_ethnicity", "age_or_mean_of_age_range", "BMI",
        "dataset", "study", "anatomical_region_ccf_score",
        "total_counts", "n_genes",
    ] if c in adata.obs.columns]
    keep_cols += score_cols
    keep_cols += [f"{g}_expr" for g in EXTRA_GENES if f"{g}_expr" in adata.obs]
    keep_cols += [f"{g}_detect" for g in EXTRA_GENES if f"{g}_detect" in adata.obs]
    adata.obs[keep_cols].to_csv(outdir / "pericytes_states_metadata.tsv.gz", sep="\t")

    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(adata.var.index.to_numpy(), adata.var[col].to_numpy()):
            adata.var.index.name = None
    adata.write(outdir / "pericyte_states.h5ad")

    session_info.show()


if __name__ == "__main__":
    main()
