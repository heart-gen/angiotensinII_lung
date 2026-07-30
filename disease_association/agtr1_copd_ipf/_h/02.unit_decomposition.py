"""Why did the GSE136831 AGTR1 fibroblast/myofibroblast direction change?

Reconstructs the OLD quantity from disease_association/ipf_analysis/_h/01
(per-patient mean of per-CELL log-normalised AGTR1, no cell floor, no
covariates, in both the all_cells and filter_cells flavours) and puts it next to
the NEW quantity from agtr1_copd_ipf (donor x compartment pseudobulk, sum raw
counts -> CP10K -> log1p, >=5-cell floor).

Per-cell normalisation here is log1p(count / nUMI * 1e4). The old script used
scran::computeSumFactors + logNormCounts; deconvolution factors are proportional
to library size up to a composition correction, so this reproduces the old
DIRECTION but not its third decimal place. That is enough to attribute a sign
flip, which is all this is for.

X in the h5ad is stored gene-major (indptr spans genes), so a single gene is one
contiguous slice -- no need to touch the other 45k rows.
"""
import h5py
import numpy as np
import pandas as pd

H5 = ("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung"
      "/disease_association/ipf_analysis/_m/ipf_dataset.h5ad")
PB = ("/ocean/projects/bio260021p/kbenjamin/projects/angiotensinII_lung"
      "/disease_association/agtr1_copd_ipf/_m/gse136831_ras_pseudobulk.tsv.gz")
COMPS = ["Fibroblast", "Myofibroblast", "Pericyte"]
GENES = ["AGTR1", "AGTR2"]


def h5str(f, path):
    """Read an anndata string column, categorical or plain."""
    node = f[path]
    if isinstance(node, h5py.Group):          # categorical
        cats = node["categories"][:].astype(str)
        return cats[node["codes"][:]]
    v = node[:]
    return v.astype(str) if v.dtype.kind != "O" else np.array(
        [x.decode() if isinstance(x, bytes) else x for x in v])


with h5py.File(H5, "r") as f:
    var_names = h5str(f, "var/_index")
    obs = pd.DataFrame({
        "donor_id": h5str(f, "obs/patient"),
        "compartment": h5str(f, "obs/Manuscript_Identity"),
        "disease": h5str(f, "obs/Disease_Identity"),
        "nUMI": f["obs/nUMI"][:].astype(np.float64),
    })
    indptr = f["X/indptr"][:]
    for g in GENES:
        (gi,) = np.where(var_names == g)
        gi = int(gi[0])
        lo, hi = indptr[gi], indptr[gi + 1]
        cell_idx = f["X/indices"][lo:hi]
        vals = f["X/data"][lo:hi]
        col = np.zeros(len(obs), dtype=np.float64)
        col[cell_idx] = vals
        obs[g] = col

print(f"cells: {len(obs)}   nonzero AGTR1 cells: {(obs.AGTR1 > 0).sum()}")

# Old quantity: per-cell CP10K log1p, then mean within (donor, compartment).
for g in GENES:
    obs[f"{g}_pc"] = np.log1p(obs[g] / obs["nUMI"] * 1e4)

d = obs[obs.compartment.isin(COMPS)].copy()
d["disease"] = d.disease.replace({"Control": "Control"})


def donor_means(frame, tag):
    out = (frame.groupby(["compartment", "disease", "donor_id"])
                .agg(n_cells=("AGTR1", "size"), val=("AGTR1_pc", "mean"))
                .reset_index())
    out["variant"] = tag
    return out


old_all = donor_means(d, "old: all cells")
# filter_cells/: the old script kept only cells with AGTR1 > 0 | AGTR2 > 0, so
# the mean is CONDITIONAL ON DETECTION and drops every dropout zero.
old_pos = donor_means(d[(d.AGTR1 > 0) | (d.AGTR2 > 0)], "old: filter_cells")

new = pd.read_csv(PB, sep="\t")
new = new[new.compartment.isin(COMPS)][
    ["compartment", "disease", "donor_id", "n_cells", "AGTR1__expr"]]
new = new.rename(columns={"AGTR1__expr": "val"})

rows = []
for tag, frame, floor in [("old: all cells", old_all, 1),
                          ("old: filter_cells", old_pos, 1),
                          ("new: pseudobulk", new.assign(variant="x"), 1),
                          ("new: pseudobulk", new.assign(variant="x"), 5)]:
    fr = frame[frame.n_cells >= floor]
    label = f"{tag} (floor {floor})"
    for comp in COMPS:
        sub = fr[fr.compartment == comp]
        m = sub.groupby("disease").val.agg(["size", "mean"])
        if "Control" not in m.index:
            continue
        ctl = m.loc["Control", "mean"]
        for arm in ["COPD", "IPF"]:
            if arm not in m.index:
                continue
            rows.append(dict(quantity=label, compartment=comp, arm=arm,
                             n_ctl=int(m.loc["Control", "size"]),
                             n_dis=int(m.loc[arm, "size"]),
                             ctl_mean=ctl, dis_mean=m.loc[arm, "mean"],
                             diff=m.loc[arm, "mean"] - ctl))

res = pd.DataFrame(rows)
pd.set_option("display.width", 200, "display.max_rows", 200)
print("\n=== unadjusted disease-minus-Control AGTR1, by quantity ===")
print(res.to_string(index=False, float_format=lambda x: f"{x:8.3f}"))

# Detection rate per arm -- the other thing the old dotplots showed.
print("\n=== fraction of cells with AGTR1 > 0 (what DotPlot's dot size shows) ===")
det = (d.groupby(["compartment", "disease"])
        .agg(cells=("AGTR1", "size"), pct_pos=("AGTR1", lambda s: 100 * (s > 0).mean()))
        .reset_index())
print(det.to_string(index=False, float_format=lambda x: f"{x:8.2f}"))

res.to_csv("unit_decomp.tsv", sep="\t", index=False)
