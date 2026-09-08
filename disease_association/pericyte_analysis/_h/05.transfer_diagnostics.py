"""Diagnostics for the scANVI label transfer (defects P2-22 and P2-23).

Two questions this module never answered about its own output, both computable
from artifacts already on disk:

  1. LAYER AUDIT (P2-22, P2-23). What is actually stored in `X` and in
     `layers["counts"]` at each stage? The pipeline's intent was `X` = log1p
     CP10K and `counts` = raw integers. Neither holds.

  2. PARTITION CONCORDANCE (P2-22, "related, unmeasured"). Does the transferred
     label agree with the query's own de novo Leiden partition? If a label
     transfer does not reproduce structure the data finds by itself, the labels
     are the model's opinion rather than a measurement.

Writes `transfer_diagnostics/` and prints everything it writes. Read-only with
respect to every other artifact.
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--clustered", default="_m/results/clustered_data.h5ad")
ap.add_argument("--query-hvg", default="_m/query_hvg.h5ad")
ap.add_argument("--ref-hvg", default="_m/ref_hvg.h5ad")
ap.add_argument("--upstream", default="../ipf_analysis/_m/ipf_dataset.h5ad")
ap.add_argument("--outdir", default="_m/transfer_diagnostics")
args = ap.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=True)


SAMPLE_N = 300_000


def audit(path, label):
    """Range/integrality of X and every layer, without loading the matrix.

    Reads the first SAMPLE_N stored non-zeros rather than the whole matrix, so
    `min`/`max` are BOUNDS on the true range, not the range itself (the full
    `clustered_data` X reaches 1290.76 where this sample sees 844.11). That is
    ample for the question being asked -- whether a matrix can be log1p CP10K,
    which a single value above ~10 already settles -- but the columns are named
    `sampled_*` so the numbers are not mistaken for exact.
    """
    rows = []
    if not Path(path).exists():
        return rows
    with h5py.File(path, "r") as f:
        names = ["X"] + [f"layers/{k}" for k in f.get("layers", {})]
        for p in names:
            node = f[p]
            data = node["data"] if isinstance(node, h5py.Group) and "data" in node else node
            v = data[:SAMPLE_N]
            v = v[v != 0]
            if not v.size:
                continue
            rows.append(dict(
                object=label, matrix=p, n_sampled=int(v.size),
                sampled_min=float(v.min()), sampled_median=float(np.median(v)),
                sampled_max=float(v.max()),
                integers_in_sample=bool(np.allclose(v, np.round(v))),
                # A log1p CP10K matrix cannot exceed ~10: log1p(1e4) = 9.21, so
                # one sampled value above that is conclusive on its own.
                consistent_with_log1p_cp10k=bool(v.max() <= 10.0),
            ))
    return rows


layers = pd.DataFrame(
    audit(args.upstream, "ipf_dataset (upstream input)")
    + audit(args.ref_hvg, "ref_hvg (scANVI training reference)")
    + audit(args.query_hvg, "query_hvg (scANVI query)")
    + audit(args.clustered, "clustered_data (shipped output)"))
layers.to_csv(out / "layer_audit.tsv", sep="\t", index=False)
print("=== layer audit (P2-22, P2-23) ===")
print(layers.to_string(index=False))

a = ad.read_h5ad(args.clustered)
o = a.obs

ct = pd.crosstab(o["leiden"], o["predicted_labels"])
ct.to_csv(out / "leiden_vs_transferred_label.tsv", sep="\t")
print("\n=== de novo Leiden x transferred label ===")
print(ct.to_string())
print("\nrow %:")
print((100 * ct.div(ct.sum(1), axis=0)).round(1).to_string())

conc = dict(
    n_cells=int(a.n_obs),
    n_leiden=int(o["leiden"].nunique()),
    n_transferred=int(o["predicted_labels"].nunique()),
    adjusted_rand_index=float(adjusted_rand_score(
        o["leiden"].astype(str), o["predicted_labels"].astype(str))),
    normalized_mutual_info=float(normalized_mutual_info_score(
        o["leiden"].astype(str), o["predicted_labels"].astype(str))),
    median_confidence=float(o["prediction_confidence"].median()),
    frac_below_half=float((o["prediction_confidence"] < 0.5).mean()),
)
(out / "concordance.json").write_text(json.dumps(conc, indent=2))
print("\n=== concordance ===")
for k, v in conc.items():
    print(f"  {k:26s} {v}")

lab = (o.groupby("predicted_labels", observed=True)["prediction_confidence"]
        .agg(n="size", median_confidence="median").reset_index())
lab["pct_of_cells"] = (100 * lab["n"] / len(o)).round(2)
lab.to_csv(out / "label_confidence.tsv", sep="\t", index=False)
print("\n=== per transferred label ===")
print(lab.round(3).to_string(index=False))

dis = o["disease"].value_counts().rename_axis("disease").reset_index(name="n_cells")
dis.to_csv(out / "disease_composition.tsv", sep="\t", index=False)
print("\n=== disease composition ===")
print(dis.to_string(index=False))
print(f"\nwrote {out}/")
