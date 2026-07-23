"""
Donor x compartment pseudobulk of the BM panel in GSE136831 (Adams/Kaminski).

This is the only dataset in the project with a real COPD arm alongside Control
and IPF. `X` is raw CSC counts with no layers, so pseudobulk is built the proper
way for counts: sum counts within a (donor x compartment) unit, divide by that
unit's total counts, scale to CP10K, log1p. That is a genuine pseudobulk library
rather than an average of per-cell normalized values.

Compartments group `Manuscript_Identity`. Pericyte is kept separate from SMC (so
the pericyte-specific report stays pericyte-specific) AND pooled as `Mural`,
because at a 5-cell floor the pooled mural compartment is the only mural-adjacent
COPD contrast with donors on both arms.

POWER, verified up front and carried into the outputs so no downstream step can
forget it: at a 5-cell floor GSE136831 has 6 COPD but only ONE Control donor with
>=5 pericytes. A pericyte-specific COPD-vs-Control test is therefore not
estimable. Powered compartments are Endothelial, Fibroblast and Myofibroblast.
"""
import numpy as np
import pandas as pd
import anndata as ad
import session_info
import logging, argparse
from pathlib import Path

COMPARTMENTS = {
    "Pericyte": "Pericyte",
    "SMC": "SMC",
    "Fibroblast": "Fibroblast",
    "Myofibroblast": "Myofibroblast",
    "VE_Capillary_A": "Endothelial", "VE_Capillary_B": "Endothelial",
    "VE_Arterial": "Endothelial", "VE_Venous": "Endothelial",
    "VE_Peribronchial": "Endothelial", "Lymphatic": "Endothelial",
    "ATI": "ATI", "ATII": "ATII",
}
# Pericyte + SMC are additionally emitted pooled under this label.
MURAL_MEMBERS = ("Pericyte", "SMC")


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="ipf_dataset.h5ad")
    p.add_argument("--genes", required=True, type=Path)
    p.add_argument("--demo", required=True, type=Path, help="sample_demo.csv")
    p.add_argument("--outfile", required=True, type=Path)
    p.add_argument("--chunk", type=int, default=40000)
    return p.parse_args()


def build(adata, gcols, genes, labels, donors, chunk=40000):
    """Sum raw counts per (donor x compartment) unit -> CP10K -> log1p.

    `labels` spans ALL cells, with None for cells outside the compartment map.
    A backed AnnData cannot be sliced twice (no view of a view), so the
    compartment filter is applied inside each chunk rather than by pre-slicing.
    """
    labels = np.asarray(labels, dtype=object)
    valid = np.array([isinstance(l, str) for l in labels])
    unit = np.where(valid,
                    labels.astype(str) + "||" + np.asarray(donors, dtype=object).astype(str),
                    None)
    codes_full = np.full(len(labels), -1, dtype=np.int64)
    lev_valid, levels = pd.factorize(pd.Series(unit[valid]))
    codes_full[valid] = lev_valid
    n_units = len(levels)

    sum_counts = np.zeros((n_units, len(genes)), dtype=np.float64)
    sum_det = np.zeros((n_units, len(genes)), dtype=np.float64)
    sum_total = np.zeros(n_units, dtype=np.float64)
    n_cells = np.zeros(n_units, dtype=np.int64)

    for start in range(0, adata.n_obs, chunk):
        stop = min(start + chunk, adata.n_obs)
        cb = codes_full[start:stop]
        keep = cb >= 0
        if not keep.any():
            continue
        block = adata[start:stop].to_memory()
        Xb = block.X
        tot = np.asarray(Xb.sum(axis=1)).ravel()[keep]
        sub = Xb[:, gcols]
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        sub = sub[keep]
        cbk = cb[keep]
        np.add.at(sum_counts, cbk, sub)
        np.add.at(sum_det, cbk, (sub > 0).astype(np.float64))
        np.add.at(sum_total, cbk, tot)
        np.add.at(n_cells, cbk, 1)
        logging.info(f"  {stop}/{adata.n_obs} cells ({keep.sum()} in map)")

    with np.errstate(divide="ignore", invalid="ignore"):
        cp10k = sum_counts / sum_total[:, None] * 1e4
    cp10k[~np.isfinite(cp10k)] = 0.0
    out = pd.DataFrame(np.log1p(cp10k), columns=[f"{g}__expr" for g in genes])
    for j, g in enumerate(genes):
        out[f"{g}__detect"] = sum_det[:, j] / n_cells
    meta = pd.DataFrame({
        "compartment": [s.split("||")[0] for s in levels],
        "donor_id": [s.split("||")[1] for s in levels],
        "n_cells": n_cells,
        "unit_total_counts": sum_total,
    })
    meta["mean_log10_counts"] = np.log10(
        np.where(n_cells > 0, sum_total / np.maximum(n_cells, 1), np.nan))
    return pd.concat([meta, out], axis=1)


def main():
    configure_logging()
    args = parse_args()
    args.outfile.parent.mkdir(parents=True, exist_ok=True)

    want = sorted(set(pd.read_csv(args.genes, sep="\t")["gene"].astype(str)))
    adata = ad.read_h5ad(args.adata, backed="r")
    logging.info(f"Opened {args.adata} ({adata.n_obs} x {adata.n_vars})")

    genes = [g for g in want if g in adata.var_names]
    missing = sorted(set(want) - set(genes))
    if missing:
        logging.warning(f"absent: {missing}")
    gcols = [adata.var_names.get_loc(g) for g in genes]

    obs = adata.obs
    comp = obs["Manuscript_Identity"].astype(str).map(COMPARTMENTS)
    donors = obs["Subject_Identity"].astype(str).to_numpy()
    disease = obs["Disease_Identity"].astype(str).to_numpy()

    lab = comp.where(comp.notna(), None).to_numpy().astype(object)
    logging.info(f"Cells in mapped compartments: "
                 f"{sum(isinstance(l, str) for l in lab)}")

    # Pass 1: named compartments. Pass 2: Pericyte+SMC pooled as Mural, so the
    # pericyte-specific report stays pericyte-specific while the pooled mural
    # compartment (the only mural contrast with donors on both COPD arms) is
    # still available.
    lab_mural = np.array(
        [("Mural" if l in MURAL_MEMBERS else None) for l in lab], dtype=object)

    frames = []
    for tag, labels in (("named", lab), ("mural", lab_mural)):
        if not any(isinstance(l, str) for l in labels):
            continue
        frames.append(build(adata, gcols, genes, labels, donors))
        logging.info(f"built {tag} pass")

    res = pd.concat(frames, ignore_index=True)

    # Attach donor-level disease + demographics. `Ever Smoker` is the single most
    # important covariate for a COPD contrast and is available here but NOT in
    # HLCA -- a genuine advantage of this dataset.
    dmap = pd.DataFrame({"donor_id": donors, "disease": disease}).drop_duplicates()
    res = res.merge(dmap, on="donor_id", how="left")

    demo = pd.read_csv(args.demo)
    demo.columns = [c.strip().strip('"') for c in demo.columns]
    demo = demo.rename(columns={"Subject": "donor_id", "Ever Smoker": "ever_smoker",
                                "Sex": "sex", "Age": "age", "Race": "race"})
    cols = [c for c in ["donor_id", "sex", "age", "race", "ever_smoker"]
            if c in demo.columns]
    res = res.merge(demo[cols].drop_duplicates("donor_id"), on="donor_id", how="left")
    logging.info(f"demographics matched for {res['sex'].notna().sum()}/{len(res)} units")

    res.to_csv(args.outfile, sep="\t", index=False)
    logging.info(f"Wrote {args.outfile}: {res.shape}")

    for floor in (5, 10, 20):
        tab = (res[res.n_cells >= floor]
               .groupby(["compartment", "disease"])["donor_id"]
               .nunique().unstack(fill_value=0))
        logging.info(f"\nDonors per compartment at >={floor} cells:\n{tab}")

    session_info.show()


if __name__ == "__main__":
    main()
