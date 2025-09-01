"""
Compare pericyte subclusters between two models (core subset vs full).
"""
import session_info
from json import dumps
from math import log as mlog
from os import path, makedirs
from typing import Tuple, Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.stats import chisquare
from statsmodels.stats.multitest import multipletests
from matplotlib.colors import ListedColormap, to_hex
from sklearn.metrics import normalized_mutual_info_score

def read_adata(path: str) -> sc.AnnData:
    return sc.read_h5ad(path)


def align_cells(adata_core: sc.AnnData, adata_full: sc.AnnData,
                cluster_key: str) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """
    Build a dataframe with aligned cell labels for intersecting cells.
    """
    intersect = adata_core.obs_names.intersection(adata_full.obs_names).tolist()
    if len(intersect) == 0:
        raise ValueError("No overlapping cell barcodes between core and full datasets.")

    core_ser = adata_core.obs.loc[intersect, cluster_key].astype("category")
    full_ser = adata_full.obs.loc[intersect, cluster_key].astype("category")

    # preserve order of appearance for nicer plotting
    core_levels = list(pd.unique(core_ser))
    full_levels = list(pd.unique(full_ser))

    df_labels = pd.DataFrame(
        {"core": pd.Categorical(core_ser, categories=core_levels, ordered=False),
         "full": pd.Categorical(full_ser, categories=full_levels, ordered=False)},
        index=intersect,
    )
    return df_labels, core_levels, full_levels


def _infer_is_disease(series: pd.Series, controls: List[str]) -> pd.Series:
    """
    Converts series to boolean for masking disease from control
    """
    s = series.astype(str).str.lower()
    return ~s.isin(controls)


def contingency_and_props(df_labels: pd.DataFrame, core_levels: List[str],
                          full_levels: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.crosstab(
        pd.Categorical(df_labels["core"], categories=core_levels),
        pd.Categorical(df_labels["full"], categories=full_levels),
        dropna=False,
    )
    props = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return counts, props


def compute_nmi(df_labels: pd.DataFrame) -> float:
    return normalized_mutual_info_score(df_labels["core"].astype(str), df_labels["full"].astype(str))


def entropy_and_purity(rowprops: pd.Series) -> Tuple[float, float]:
    """Shannon entropy (base e) and normalized purity = 1 - H/log(C)."""
    p = rowprops.values.astype(float)
    p = p[p > 0]
    if p.size == 0:
        return 0.0, 0.0
    H = -(p * np.log(p)).sum()
    C = len(rowprops)
    Hmax = mlog(C) if C > 1 else 1.0
    purity = 1.0 - (H / Hmax)
    return H, purity


def per_cluster_stats(counts: pd.DataFrame) -> pd.DataFrame:
    """
    For each core cluster (row), compute entropy, purity, chi-square vs expected
    (expected based on column marginals among overlapping cells).
    """
    col_totals = counts.sum(axis=0).values.astype(float)
    grand_total = col_totals.sum()
    
    rows = []
    for core, obs in counts.iterrows():
        obs = obs.astype(float)
        row_sum = float(obs.sum())
        if row_sum == 0:
            H, purity = 0.0, 0.0
            chi2, p = np.nan, np.nan
        else:
            # Entropy & purity from proportions
            rowprops = obs / row_sum
            H, purity = entropy_and_purity(rowprops)

            # Chi-square vs expected from column marginals
            expected = row_sum * (col_totals / grand_total)
            # Avoid zero
            eps = 1e-12
            expected = np.where((expected == 0) & (obs.values > 0), eps, expected)
            chi2, p = chisquare(f_obs=obs.values, f_exp=expected)
        rows.append({"core_cluster": core, "entropy": H, "purity": purity,
                     "chi2": chi2, "pval": p})

    df = pd.DataFrame(rows).set_index("core_cluster")
    # BH-FDR across core clusters
    mask = df["pval"].notna()
    if mask.any():
        _, q, _, _ = multipletests(df.loc[mask, "pval"], method="fdr_bh")
        df.loc[mask, "fdr_bh"] = q
    else:
        df["fdr_bh"] = np.nan
    return df


def top_matches(rowprops: pd.Series, k=2) -> List[Tuple[str, float]]:
    order = rowprops.sort_values(ascending=False)
    return list(zip(order.index[:k], order.values[:k]))


def save_heatmap(props: pd.DataFrame, full_disease_frac: pd.Series, out_png: str,
                 out_pdf: str):
    """
    Heatmap of row-wise proportions with column labels augmented by disease %.
    """
    # augment column labels
    col_labels = [
        f"{col}{'*' if (full_disease_frac.get(col, 0) >= 0.5) else ''}\n({full_disease_frac.get(col, np.nan)*100:.0f}% dis.)"
        for col in props.columns
    ]

    fig, ax = plt.subplots(figsize=(1.4*props.shape[1] + 3, 0.5*props.shape[0] + 3))
    im = ax.imshow(props.values, aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(props.shape[1]))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(range(props.shape[0]))
    ax.set_yticklabels(props.index)
    ax.set_xlabel("Full subclusters  (percent disease in full, all cells)")
    ax.set_ylabel("Core subclusters")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Proportion of core cluster cells")

    # annotate with percentages
    for i in range(props.shape[0]):
        for j in range(props.shape[1]):
            val = props.values[i, j]
            if val >= 0.02:  # annotate only non-trivial
                ax.text(j, i, f"{val*100:.0f}%", ha="center", va="center",
                        fontsize=8, color="white" if val > 0.5 else "black")

    fig.tight_layout()
    for path in [out_png, out_pdf]:
        fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _interp_color(a: float, c0=(0, 0, 1), c1=(1, 0, 0)) -> str:
    """Linear interpolate between blue (control) and red (disease) in RGB, return hex."""
    a = float(np.clip(a, 0, 1))
    r = c0[0] + a * (c1[0] - c0[0])
    g = c0[1] + a * (c1[1] - c0[1])
    b = c0[2] + a * (c1[2] - c0[2])
    return to_hex((r, g, b))


def make_sankey(counts: pd.DataFrame, full_disease_frac: pd.Series, out_html: str,
                out_png: str = None):
    core_nodes = list(counts.index)
    full_nodes = list(counts.columns)

    node_labels = (
        [f"core:{c}" for c in core_nodes] +
        [f"full:{f}{'*' if (full_disease_frac.get(f, 0) >= 0.5) else ''}" for f in full_nodes]
    )
    # node colors: core=gray, full=blue->red by disease frac
    node_colors = (
        ["#bdbdbd"] * len(core_nodes) +
        [
            "#800000" if (full_disease_frac.get(f, 0) >= 0.5) else _interp_color(full_disease_frac.get(f, 0.0))
            for f in full_nodes
        ]
    )

    # links
    src, tgt, val, link_colors = [], [], [], []
    for i, c in enumerate(core_nodes):
        for j, f in enumerate(full_nodes):
            v = counts.loc[c, f]
            if v <= 0:
                continue
            src.append(i)                           # core node index
            tgt.append(len(core_nodes) + j)         # full node index
            val.append(float(v))
            # color like destination node, semi-transparent
            link_colors.append(_interp_color(full_disease_frac.get(f, 0.0)) + "80")  # add alpha

    fig = go.Figure(
        data=[go.Sankey(
            arrangement="snap",
            node=dict(label=node_labels, color=node_colors, pad=14, thickness=16),
            link=dict(source=src, target=tgt, value=val, color=link_colors),
        )]
    )
    fig.update_layout(
        title="Core -> Full pericyte subcluster mapping",
        font=dict(size=12)
    )
    fig.write_html(out_html)
    if out_png:
        try:
            fig.write_image(out_png, scale=2)
        except Exception:
            print("Kaleido not available; saved HTML only.")


def plot_umaps(adata_full: sc.AnnData, df_labels: pd.DataFrame, cluster_key: str,
               disease_bool: pd.Series, outdir: str):
    umap_key = "X_umap"
    X = adata_full.obsm[umap_key]
    df_full = adata_full.obs[[cluster_key]].copy()
    df_full["is_disease"] = disease_bool.values if disease_bool.index.equals(adata_full.obs.index) \
        else _infer_is_disease(adata_full.obs.loc[:, disease_bool.name])
    df_full["disease"] = adata_full.obs.disease.values
    
    # color by full cluster
    cats = pd.Categorical(df_full[cluster_key])
    n_cat = len(cats.categories)
    cmap = plt.get_cmap("tab20", max(n_cat, 3))
    colors = [cmap(i % cmap.N) for i in range(n_cat)]
    color_map = dict(zip(cats.categories, colors))
    point_colors = [color_map[c] for c in cats]

    # UMAP by full cluster
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(X[:, 0], X[:, 1], s=6, c=point_colors, linewidths=0, alpha=0.9)
    handles = [plt.Line2D([], [], marker='o', linestyle='', color=color_map[c],
                          label=c, markersize=6) for c in cats.categories]
    ax.legend(handles=handles, loc='upper right', title=full_cluster_key,
              fontsize=8, frameon=False, ncol=1)
    ax.set_title("UMAP -- Full pericyte subclusters")

    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.tight_layout()
    fig.savefig(path.join(outdir, "umap_full_by_fullcluster.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(path.join(outdir, "umap_full_by_fullcluster.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Overlay core cells as outlines on top of full (by full cluster colors)
    is_core = adata_full.obs_names.isin(df_labels.index)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(X[~is_core, 0], X[~is_core, 1], s=5,
               c=np.array(point_colors, dtype=object)[~is_core], alpha=0.25, linewidths=0)
    ax.scatter(X[is_core, 0], X[is_core, 1], s=12, facecolors='none',
               edgecolors='black', linewidths=0.5, alpha=0.9)
    ax.set_title("UMAP -- Full subclusters with Core cells overlaid")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.tight_layout()
    fig.savefig(path.join(outdir, "umap_full_with_core_overlay.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(path.join(outdir, "umap_full_with_core_overlay.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # UMAP by disease (general)
    fig, ax = plt.subplots(figsize=(7, 6))
    c = np.array(["#377eb8", "#e41a1c"])  # control=blue, disease=red
    disease_arr = df_full["is_disease"].astype(bool).values
    ax.scatter(X[:, 0], X[:, 1], s=6, c=c[disease_arr.astype(int)], linewidths=0, alpha=0.9)
    handles = [
        plt.Line2D([], [], marker='o', linestyle='', color="#377eb8",
                   label="control", markersize=6),
        plt.Line2D([], [], marker='o', linestyle='', color="#e41a1c",
                   label="disease", markersize=6),
    ]
    ax.legend(handles=handles, loc='upper right', title="condition",
              fontsize=8, frameon=False, ncol=1)
    ax.set_title("UMAP -- Disease highlighting")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.tight_layout()
    fig.savefig(path.join(outdir, "umap_full_by_disease.png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.join(outdir, "umap_full_by_disease.pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    # UMAP by disease (specific)
    fig, ax = plt.subplots(figsize=(7, 6))
    disease_cat = pd.Categorical(df_full["disease"])
    categories = disease_cat.categories
    # Pick a colormap with enough distinct colors
    cmap = plt.get_cmap("tab20", len(categories))
    colors = [cmap(i) for i in range(len(categories))]
    # Map colors
    color_map = dict(zip(categories, colors))
    point_colors = [color_map[val] for val in disease_cat]
    ax.scatter(X[:, 0], X[:, 1], s=6, c=point_colors, linewidths=0, alpha=0.9)
    handles = [
        plt.Line2D([], [], marker='o', linestyle='', color=color_map[cat],
                   label=cat, markersize=6)
        for cat in categories
    ]
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
              title="Disease", fontsize=7, frameon=False, ncol=2)
    ax.set_title("UMAP -- Disease")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    fig.tight_layout()
    fig.savefig(path.join(outdir, "umap_full_by_disease_categories.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(path.join(outdir, "umap_full_by_disease_categories.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def run(core_path: str, full_path: str, cluster_key: str, outdir: str):
    # Create output directory
    makedirs(outdir, exist_ok=True)

    adata_core = read_adata(core_path)
    adata_full = read_adata(full_path)

    # Align labels on overlapping cells
    df_labels, core_levels, full_levels = align_cells(
        adata_core, adata_full, cluster_key
    )

    # Disease (computed across ALL full cells, not just overlap)
    disease_bool_full = _infer_is_disease(adata_full.obs["disease"], ["normal"])
    # Disease fraction per full cluster (all full cells)
    disease_frac_by_full = (
        pd.DataFrame({
            cluster_key: adata_full.obs[cluster_key].astype(str).values,
            "is_disease": disease_bool_full.values,
        })
        .groupby(cluster_key)["is_disease"]
        .mean()
    )
    disease_specific_clusters = disease_frac_by_full[disease_frac_by_full >= 0.5].index.tolist()
    print("Disease-specific full clusters (>50% disease):",
          disease_specific_clusters)
    disease_frac_by_full.to_csv(path.join(outdir,
                                          "full_clusters_disease_fraction.csv"))

    # Contingency + proportions (on overlapping cells only)
    counts, props = contingency_and_props(df_labels, core_levels, full_levels)

    # Global metric
    nmi = compute_nmi(df_labels)
    with open(path.join(outdir, "global_metrics.txt"), "w") as fh:
        fh.write(dumps({"NMI_core_vs_full": float(nmi)}, indent=2))
    print(f"NMI (core vs full on overlapping cells): {nmi:.3f}")

    # Per-core cluster stats
    stats_df = per_cluster_stats(counts)
    stats_df.to_csv(path.join(outdir, "per_cluster_metrics.csv"))

    # Top matches table
    top_rows = []
    for core in props.index:
        tm = top_matches(props.loc[core], k=2)
        row = {"core_cluster": core}
        for i, (f, p) in enumerate(tm, 1):
            row[f"top{i}_full_cluster"] = f
            row[f"top{i}_prop"] = float(p)
        top_rows.append(row)
    pd.DataFrame(top_rows).set_index("core_cluster").to_csv(
        path.join(outdir, "mapping_top_matches.csv")
    )

    # Save matrices
    counts.to_csv(path.join(outdir, "contingency_counts.csv"))
    props.to_csv(path.join(outdir, "contingency_rowprop.csv"))

    # Heatmap with disease fractions in column labels
    save_heatmap(
        props,
        full_disease_frac=disease_frac_by_full,
        out_png=path.join(outdir, "heatmap_proportions.png"),
        out_pdf=path.join(outdir, "heatmap_proportions.pdf"),
    )

    # Sankey (full node colors by disease fraction)
    make_sankey(
        counts,
        full_disease_frac=disease_frac_by_full,
        out_html=path.join(outdir, "sankey_core_to_full.html"),
        out_png=path.join(outdir, "sankey_core_to_full.png"),
    )

    # UMAPs (full dataset), disease highlighting, and core overlay
    # Build disease bool aligned to full.obs
    disease_bool_aligned = pd.Series(disease_bool_full.values,
                                     index=adata_full.obs_names,
                                     name="disease")
    plot_umaps(adata_full=adata_full, df_labels=df_labels, cluster_key=cluster_key,
               disease_bool=disease_bool_aligned, outdir=outdir)

    print(f"Done. Results written to: {outdir}")


def main():
    # Main analysis
    core_path = "../../_m/pericyte.hlca_core.subclustered.h5ad"
    full_path = "../../_m/pericyte.hlca_full.subclustered.h5ad"
    run(core_path=core_path, full_path=full_path, cluster_key='leiden',
        outdir="./")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
