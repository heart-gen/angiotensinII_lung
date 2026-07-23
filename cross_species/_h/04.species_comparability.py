"""
Cross-species comparability diagnostic: is the mouse mural set a valid
counterpart to the human pericyte analysis?

Motivation. The human analysis (`pericyte_states/`) clusters 11,680 *pericytes*
and labels clusters by RELATIVE enrichment (z-score each program across cells,
average within cluster, argmax). The mouse pipeline (`02.conserved_states.py`)
clusters the *whole lung* and labels clusters by RAW score magnitude. Before any
conservation claim is made, we quantify how far apart those two designs are.

This script is read-only and additive: it does not modify the mouse state
assignment or any published output. It writes a set of diagnostic tables that
the manuscript's cross-species Methods/Limitations cite directly.

Outputs (all under --outdir):
  species_comparability_celltype_by_dataset.tsv  cell type x dataset (the confound)
  species_comparability_state_by_celltype.tsv    published state x cell type
  species_comparability_agtr1a_by_layer.tsv      Agtr1a per cell type, per layer
  species_comparability_agtr1a_within_dataset.tsv  within-dataset pericyte vs SMC
  species_comparability_agtr1a_tests.tsv         Fisher / CMH on detection
  species_comparability_relenrich_human_algo.tsv  human algorithm applied to mouse
  species_comparability_mural_cells.tsv           per-cell table, all mural cells
  species_comparability_pericyte_cells.tsv        per-cell table, pericytes only
  species_comparability_pericyte_by_donor.tsv     per-donor pericyte detection
  species_comparability_summary.tsv              headline numbers for the text
"""
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from scipy import sparse, stats

MURAL_TYPES = ["pericyte", "vascular associated smooth muscle cell",
               "smooth muscle cell of the pulmonary artery"]
PERICYTE = "pericyte"
# Layer preference for an *expression* claim. `scvi_corrected` is a dense
# denoised reconstruction: every cell is non-zero, so detection fractions
# computed on it are meaningless (they are trivially 1.0). Raw claims must use
# `lognorm`.
EXPR_LAYERS = ["lognorm", "scvi_corrected"]
RECEPTORS = ["Agtr1a", "Agtr1b", "Agtr2"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="mouse_integrated.h5ad")
    p.add_argument("--metadata", required=True, type=Path,
                   help="mouse_states_metadata.tsv.gz from 02.conserved_states.py")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def write_tsv(df, path, index=False):
    df.to_csv(path, sep="\t", index=index)
    logging.info(f"wrote {path.name}")


def relative_enrichment(obs, score_cols, cluster_key):
    """Reproduce the HUMAN labelling rule (pericyte_states/_h/00.state_discovery.py
    `annotate_states`): z-score each program across cells, average within cluster,
    then argmax. Raw magnitudes are not comparable across panels of different size
    and baseline expression, which is why the human side does not use them."""
    progs = [c.replace("_score", "") for c in score_cols]
    z = obs[score_cols].apply(lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9))
    z = z.assign(_clu=obs[cluster_key].to_numpy())
    rel = z.groupby("_clu", observed=True)[score_cols].mean()
    rel.columns = [f"{p}_relenrich" for p in progs]
    rel["assigned"] = [progs[i] for i in rel.to_numpy().argmax(axis=1)]
    return rel


def fetch_gene(adata, gene, layer):
    x = adata[:, gene].layers[layer]
    return x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)
    mural = meta[meta["is_mural"] == True].copy()  # noqa: E712
    mural["cell_type"] = mural["cell_type"].astype(str)
    mural["dataset_id"] = mural["dataset_id"].astype(str)
    logging.info(f"mural cells: {len(mural)}  pericytes: "
                 f"{int((mural.cell_type == PERICYTE).sum())}")

    # ---- 1. The confound: cell type x dataset -----------------------------
    ct_ds = pd.crosstab(mural["cell_type"], mural["dataset_id"])
    write_tsv(ct_ds.reset_index(),
              args.outdir / "species_comparability_celltype_by_dataset.tsv")
    # A dataset is "informative" for a pericyte-vs-SMC contrast only if it
    # contains pericytes AND at least one SMC type. Otherwise cell type is
    # aliased with batch and the contrast is not estimable within dataset.
    informative = [d for d in ct_ds.columns
                   if ct_ds.loc[PERICYTE, d] > 0 and
                   ct_ds[d].drop(PERICYTE).sum() > 0]
    logging.info(f"datasets supporting a within-dataset pericyte-vs-SMC "
                 f"contrast: {len(informative)} of {ct_ds.shape[1]}")

    # ---- 2. What the published states are actually made of ----------------
    state_ct = pd.crosstab(mural["pericyte_state"], mural["cell_type"])
    write_tsv(state_ct.reset_index(),
              args.outdir / "species_comparability_state_by_celltype.tsv")

    # ---- 3. Human labelling rule applied to the mouse mural cells ---------
    score_cols = [c for c in mural.columns if c.endswith("_score")]
    rel = relative_enrichment(mural, score_cols, "leiden_scvi")
    rel["n_mural"] = mural.groupby("leiden_scvi", observed=True).size()
    write_tsv(rel.reset_index(),
              args.outdir / "species_comparability_relenrich_human_algo.tsv")
    human_algo_label = mural["leiden_scvi"].map(rel["assigned"])
    comp = pd.concat([
        mural["pericyte_state"].value_counts().rename("raw_magnitude_published"),
        human_algo_label.value_counts().rename("relative_enrichment_human_rule"),
    ], axis=1).fillna(0).astype(int)

    # ---- 4. Agtr1a by cell type, raw vs denoised --------------------------
    adata = sc.read_h5ad(args.adata, backed="r")
    keep = np.where(adata.obs["cell_type"].astype(str).isin(MURAL_TYPES).to_numpy())[0]
    sub = adata[keep].to_memory()
    layers = [l for l in EXPR_LAYERS if l in sub.layers]
    logging.info(f"layers available: {layers}")

    rows, within = [], []
    expr = {}
    for gene in [g for g in RECEPTORS if g in sub.var_names]:
        for layer in layers:
            v = fetch_gene(sub, gene, layer)
            expr[(gene, layer)] = v
            d = pd.DataFrame({"v": v,
                              "cell_type": sub.obs["cell_type"].astype(str).to_numpy(),
                              "dataset_id": sub.obs["dataset_id"].astype(str).to_numpy()})
            g_all = d.groupby("cell_type")["v"].agg(
                n="size", mean="mean", detect=lambda s: float((s > 0).mean()))
            rows.append(g_all.assign(gene=gene, layer=layer).reset_index())
            for ds in informative:
                g_ds = d[d.dataset_id == ds].groupby("cell_type")["v"].agg(
                    n="size", mean="mean", detect=lambda s: float((s > 0).mean()))
                within.append(g_ds.assign(gene=gene, layer=layer,
                                          dataset_id=ds).reset_index())
    write_tsv(pd.concat(rows, ignore_index=True),
              args.outdir / "species_comparability_agtr1a_by_layer.tsv")
    write_tsv(pd.concat(within, ignore_index=True),
              args.outdir / "species_comparability_agtr1a_within_dataset.tsv")

    # ---- 5. Formal test on DETECTION, stratified by dataset ---------------
    # Detection (not mean) because one arm is exactly zero, so a linear model is
    # completely separated. Fisher per dataset + Mantel-Haenszel across them.
    tests = []
    if ("Agtr1a", "lognorm") in expr:
        d = pd.DataFrame({
            "det": (expr[("Agtr1a", "lognorm")] > 0).astype(int),
            "cell_type": sub.obs["cell_type"].astype(str).to_numpy(),
            "dataset_id": sub.obs["dataset_id"].astype(str).to_numpy()})
        tables = []
        for ds in informative:
            g = d[d.dataset_id == ds]
            per = g[g.cell_type == PERICYTE]["det"]
            smc = g[g.cell_type != PERICYTE]["det"]
            tab = np.array([[int(per.sum()), int((1 - per).sum())],
                            [int(smc.sum()), int((1 - smc).sum())]])
            tables.append(tab)
            _, p = stats.fisher_exact(tab)
            tests.append({"test": "fisher_exact_detection", "dataset_id": ds,
                          "n_pericyte": int(per.size), "n_smc": int(smc.size),
                          "detect_pericyte": float(per.mean()),
                          "detect_smc": float(smc.mean()), "p_value": float(p)})
        if len(tables) > 1:
            # Mantel-Haenszel chi-square with continuity correction; pools the
            # per-dataset 2x2 tables so batch cannot drive the result. The
            # descriptive columns are computed from the SAME strata that enter
            # the test (informative datasets only) -- reporting whole-cohort
            # counts next to a stratified p-value would misdescribe the test.
            num = sum(t[0, 0] - (t[0].sum() * t[:, 0].sum()) / t.sum() for t in tables)
            den = sum((t[0].sum() * t[1].sum() * t[:, 0].sum() * t[:, 1].sum()) /
                      (t.sum() ** 2 * (t.sum() - 1)) for t in tables)
            chi = (abs(num) - 0.5) ** 2 / den if den > 0 else np.nan
            n_p = sum(int(t[0].sum()) for t in tables)
            n_s = sum(int(t[1].sum()) for t in tables)
            tests.append({"test": "mantel_haenszel_detection_stratified_by_dataset",
                          "dataset_id": "|".join(informative),
                          "n_pericyte": n_p, "n_smc": n_s,
                          "detect_pericyte": sum(int(t[0, 0]) for t in tables) / n_p,
                          "detect_smc": sum(int(t[1, 0]) for t in tables) / n_s,
                          "p_value": float(stats.chi2.sf(chi, 1))})
    write_tsv(pd.DataFrame(tests),
              args.outdir / "species_comparability_agtr1a_tests.tsv")

    # ---- 5b. Cell-level view of every mouse pericyte ----------------------
    # n is small enough (tens of cells) to inspect exhaustively rather than
    # trusting a mean. Reports every layer side by side so a denoised non-zero
    # cannot be mistaken for an observed count.
    allc = pd.DataFrame(index=sub.obs_names)
    allc["cell_type"] = sub.obs["cell_type"].astype(str).to_numpy()
    for c in ["dataset_id", "donor_id", "sex", "disease"]:
        if c in sub.obs.columns:
            allc[c] = sub.obs[c].astype(str).to_numpy()
    if "counts" in sub.layers:
        allc["total_counts_cell"] = np.asarray(
            sub.layers["counts"].sum(axis=1)).ravel()
    for gene in [g for g in RECEPTORS if g in sub.var_names]:
        for layer in [l for l in ["counts"] + EXPR_LAYERS if l in sub.layers]:
            allc[f"{gene}_{layer}"] = fetch_gene(sub, gene, layer)
    allc = allc.sort_values([c for c in ["dataset_id", "cell_type", "donor_id"]
                             if c in allc])
    # Whole mural compartment (drives the supplement figure) ...
    write_tsv(allc.reset_index().rename(columns={"index": "cell_id"}),
              args.outdir / "species_comparability_mural_cells.tsv")
    # ... and the pericyte-only view cited in the manuscript text.
    cells = allc[allc.cell_type == PERICYTE].drop(columns=["cell_type"])
    write_tsv(cells.reset_index().rename(columns={"index": "cell_id"}),
              args.outdir / "species_comparability_pericyte_cells.tsv")

    # Per-donor rollup: donor is the unit of inference elsewhere in the repo, so
    # report how many donors show ANY Agtr1a+ pericyte, not just the cell count.
    if "Agtr1a_lognorm" in cells.columns and "donor_id" in cells.columns:
        donor = (cells.assign(det=(cells["Agtr1a_lognorm"] > 0).astype(int))
                 .groupby(["dataset_id", "donor_id"])
                 .agg(n_pericytes=("det", "size"), n_agtr1a_pos=("det", "sum"),
                      mean_lognorm=("Agtr1a_lognorm", "mean"))
                 .reset_index())
        donor["any_pos"] = donor["n_agtr1a_pos"] > 0
        write_tsv(donor, args.outdir / "species_comparability_pericyte_by_donor.tsv")
        logging.info(f"donors with >=1 Agtr1a+ pericyte: "
                     f"{int(donor.any_pos.sum())} of {len(donor)}")

    # ---- 6. Headline summary ---------------------------------------------
    n_per = int((mural.cell_type == PERICYTE).sum())
    top_clu = mural.groupby("leiden_scvi", observed=True).size().max()
    summary = pd.DataFrame([
        {"metric": "n_mural_cells", "value": len(mural)},
        {"metric": "n_pericytes", "value": n_per},
        {"metric": "pct_mural_that_are_pericytes", "value": round(100 * n_per / len(mural), 2)},
        {"metric": "n_pericyte_donors", "value": int(mural[mural.cell_type == PERICYTE].donor_id.nunique())},
        {"metric": "n_pericyte_datasets", "value": int(mural[mural.cell_type == PERICYTE].dataset_id.nunique())},
        {"metric": "n_datasets_supporting_pericyte_vs_smc", "value": len(informative)},
        {"metric": "pct_mural_in_single_largest_whole_lung_cluster",
         "value": round(100 * top_clu / len(mural), 2)},
        {"metric": "n_human_pericytes_for_comparison", "value": 11680},
    ])
    write_tsv(summary, args.outdir / "species_comparability_summary.tsv")

    logging.info("\n--- published (raw magnitude) vs human rule (relative enrichment) ---\n"
                 + comp.to_string())
    logging.info("\n--- headline ---\n" + summary.to_string(index=False))
    session_info.show()


if __name__ == "__main__":
    main()
