"""
Figure 5B -- decompose the pericyte "airspace score" into its four compartments.

localization/airspace_analysis/_h/01.airspace_analysis.py scores each pericyte as
the MEAN of four cosine similarities, in X_pca_harmony, to the centroids of AT1,
AT2, EC aerocyte capillary and EC general capillary, and discards the four values.
Its AGTR1 association is a null (beta -0.0020, P 0.79). A mean of four can hide
compartment-specific structure, so this script recomputes the four similarities
EXACTLY as that script does and keeps them.

Vocabulary: this is transcriptional similarity in an integrated latent space --
"niche affinity" -- not anatomical proximity. Never call it proximity.

Two reference axes are added, labelled as such and kept out of the four-axis
family: similarity to Smooth muscle (the mural neighbour; AGTR1 is a mural
compartment label, so a positive slope here is the expected direction and serves
as a positive control) and to EC arterial.

Also emitted, from pericyte_states.h5ad (the object every other AGTR1 analysis
uses): donor-level pseudobulk (expm1 -> mean -> log1p) and detection for AGTR1
and for EVERY gene whose pericyte detection lies within +/- --tol of AGTR1's.
02.niche_affinity_stats.R refits the same models with each matched gene in
place of AGTR1, which gives the detection-matched null that a single-gene
exposure needs (memory: gene-set-score-null-not-zero -- take every gene inside a
tolerance, not the K nearest, or the null's mean detection is biased).

Outputs (to --outdir):
  pericyte_niche_affinity.tsv.gz      per pericyte: ids, covariates, 4+2 affinities
  niche_affinity_gene_pseudobulk.tsv.gz  donor x gene: expr, detect, n_cells, role
  niche_affinity_centroid_sizes.tsv   cells per reference class
"""
import argparse
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

FOUR = {"AT1": "AT1", "AT2": "AT2",
        "EC_aerocyte": "EC aerocyte capillary",
        "EC_gcap": "EC general capillary"}
REFERENCE = {"SMC": "Smooth muscle", "EC_arterial": "EC arterial"}


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--airspace", type=Path, required=True,
                   help="localization/airspace_analysis/_m/pericytes_with_airspace_score.h5ad")
    p.add_argument("--pericytes", type=Path, required=True,
                   help="pericyte_states/_m/pericyte_states.h5ad")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--rep", default="X_pca_harmony")
    p.add_argument("--key", default="subclusters")
    p.add_argument("--tol", type=float, default=0.03,
                   help="detection tolerance around AGTR1 for the matched null")
    return p.parse_args()


def unit_rows(M):
    return M / np.linalg.norm(M, axis=1, keepdims=True)


def main():
    a = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    a.outdir.mkdir(parents=True, exist_ok=True)

    # ---- affinities: obs + obsm only; X stays on disk ---------------------
    air = ad.read_h5ad(a.airspace, backed="r")
    obs = air.obs.copy()
    X = np.asarray(air.obsm[a.rep])
    lab = obs[a.key].astype(str).to_numpy()
    logging.info(f"airspace object: {air.n_obs} cells, {a.rep} {X.shape}")

    classes = {**FOUR, **REFERENCE}
    sizes, cent = [], {}
    for short, cls in classes.items():
        m = lab == cls
        sizes.append({"axis": short, "class": cls, "n_cells": int(m.sum()),
                      "family": "four_axis" if short in FOUR else "reference"})
        if m.sum() == 0:
            raise KeyError(f"no cells labelled '{cls}' in obs['{a.key}']")
        # Mean of the RAW latent vectors, then normalise -- identical to
        # compute_centroids() + compute_airspace_scores().
        cent[short] = X[m].mean(axis=0)
    pd.DataFrame(sizes).to_csv(a.outdir / "niche_affinity_centroid_sizes.tsv",
                               sep="\t", index=False)

    peri = lab == "Pericytes"
    V = unit_rows(X[peri])
    C = unit_rows(np.stack([cent[k] for k in classes]))
    S = V @ C.T                                         # cells x axes
    aff = pd.DataFrame(S, columns=[f"affinity_{k}" for k in classes],
                       index=obs.index[peri])

    # ---- regression gate: the mean of the four must BE the published score ---
    mean4 = aff[[f"affinity_{k}" for k in FOUR]].mean(axis=1).to_numpy()
    pub = obs.loc[peri, "airspace_score"].to_numpy(dtype=float)
    dev = np.nanmax(np.abs(mean4 - pub))
    logging.info(f"gate: max |mean(4 affinities) - airspace_score| = {dev:.2e}")
    if not np.allclose(mean4, pub, atol=1e-5, equal_nan=False):
        raise AssertionError(
            f"the four affinities do not reproduce airspace_score (max dev {dev:.3g}). "
            "This object is not the one the published null was computed on, or "
            "the centroid definition has drifted. Stop.")

    keep = ["donor_id", "study", "dataset", "sex", "age_or_mean_of_age_range",
            "lung_condition", "log10_total_counts", "AGTR1_expr", "AGTR1_detect",
            "airspace_score"]
    keep = [c for c in keep if c in obs.columns]
    out = pd.concat([obs.loc[peri, keep], aff], axis=1)
    out.index.name = "index"

    # ---- same cells in pericyte_states? ------------------------------------
    pst = sc.read_h5ad(a.pericytes)
    if "logcounts" not in pst.layers:
        raise KeyError("pericyte_states.h5ad has no 'logcounts' layer")
    if "feature_name" in pst.var.columns:
        sym = pst.var["feature_name"].astype(str)
        if not (sym.values == pst.var_names.values).all():
            pst.var_names = sym
            pst.var_names_make_unique()
    shared = out.index.intersection(pst.obs_names)
    frac = len(shared) / len(out)
    logging.info(f"pericyte barcodes matched to pericyte_states: "
                 f"{len(shared)}/{len(out)} ({100 * frac:.2f}%)")
    if frac < 0.99:
        raise AssertionError("fewer than 99% of airspace pericytes are in "
                             "pericyte_states.h5ad; the two objects disagree")
    out = out.loc[shared]
    if "pericyte_state" in pst.obs.columns:
        out["pericyte_state"] = pst.obs.loc[shared, "pericyte_state"].astype(str)
    out.to_csv(a.outdir / "pericyte_niche_affinity.tsv.gz", sep="\t")
    logging.info(f"wrote {len(out)} pericytes, {out['donor_id'].nunique()} donors")

    # ---- AGTR1 + detection-matched genes, donor pseudobulk -----------------
    pst = pst[shared].copy()
    L = pst.layers["logcounts"]
    L = L.tocsc() if sparse.issparse(L) else sparse.csc_matrix(L)
    detect = np.asarray((L > 0).mean(axis=0)).ravel()
    det = pd.Series(detect, index=pst.var_names)
    if "AGTR1" not in det.index:
        raise KeyError("AGTR1 absent from pericyte_states.h5ad")
    d0 = float(det["AGTR1"])
    pool = det[(det - d0).abs() <= a.tol].index.difference(["AGTR1"])
    logging.info(f"AGTR1 pericyte detection {d0:.4f}; {len(pool)} genes within "
                 f"+/-{a.tol} (achieved mean {det[pool].mean():.4f})")

    genes = ["AGTR1"] + list(pool)
    gidx = [pst.var_names.get_loc(g) for g in genes]
    sub = L[:, gidx]
    donors = pst.obs["donor_id"].astype(str).to_numpy()
    rows = []
    for dnr in pd.unique(donors):
        m = donors == dnr
        blk = sub[m]
        n = int(m.sum())
        lin = np.asarray(blk.expm1().mean(axis=0)).ravel()
        dt = np.asarray((blk > 0).mean(axis=0)).ravel()
        for g, e, f in zip(genes, np.log1p(lin), dt):
            rows.append((dnr, g, e, f, n))
    pb = pd.DataFrame(rows, columns=["donor_id", "gene", "expr", "detect", "n_cells"])
    pb["role"] = np.where(pb["gene"] == "AGTR1", "observed", "null")
    pb["pericyte_detect_all"] = pb["gene"].map(det)
    pb.to_csv(a.outdir / "niche_affinity_gene_pseudobulk.tsv.gz", sep="\t", index=False)
    logging.info(f"gene pseudobulk: {pb['donor_id'].nunique()} donors x {len(genes)} genes")


if __name__ == "__main__":
    main()
