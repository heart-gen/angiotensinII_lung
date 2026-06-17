"""
Pericyte-endothelial niche-stability index (donor-level).

Links transcriptomic pericyte states to the in-vivo vascular phenotype by
combining, per donor:

  Niche-stability components (vascular support):
    + mean vascular_stabilizing module score
    + mean airspace-proximity score (closeness to AT1/AT2/EC niche)

  Injury-stromal components (destabilization), PRIMARY composite:
    + mean inflammatory + fibroblast_like + activated_migratory module scores
    + injury-program state fraction (fraction of pericytes in stable clusters
      annotated as injury programs)

  AGTR1+ pericyte fraction is deliberately EXCLUDED from the primary injury
  composite for two reasons: (i) AGTR1 detection is dropout-prone in pericytes --
  as a low-abundance GPCR, a binary AGTR1+ donor fraction largely tracks per-cell
  sequencing depth rather than true receptor abundance (the same dropout that made
  binary AGTR1+/- receivers unusable in the CCC analysis), so it is a technically
  unreliable composite component; and (ii) AGTR1 is the study's focal receptor, so
  folding its detection rate into the disease-association outcome risks a
  circular/confounded readout. It is retained as a separate column and combined
  into a clearly-labelled SENSITIVITY composite (injury_stromal_score_sens_agtr1 /
  niche_index_sens_agtr1) for the supplement only -- not the main narrative/figures.

Components use the curated module scores at face value (the NVU-pattern states
from 00.state_discovery.py are clustered on the study-integrated embedding, so
study is handled once at clustering -- the scores are not additionally batch-
standardized). Each component is z-scored across donors so stability and injury
terms combine on a common scale; the composite scores are the mean of their
z-scored components, and niche_index = stability - injury. Tests of the disease
effect live in 01.niche_disease_stats.R.

Inputs (all light): pericytes_states_metadata.tsv.gz (Module 2),
airspace_donor_summary.csv (existing airspace analysis).
"""
import numpy as np
import pandas as pd
import logging, argparse
from pathlib import Path

INJURY_PROGRAMS = ["inflammatory", "fibroblast_like", "activated_migratory"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--states-meta", required=True, type=Path,
                   help="pericytes_states_metadata.tsv.gz")
    p.add_argument("--airspace-summary", required=True, type=Path,
                   help="airspace_donor_summary.csv")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-cells", type=int, default=20,
                   help="Min pericytes per donor to include")
    return p.parse_args()


def zscore(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / (s.std(ddof=0) + 1e-9)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.states_meta, sep="\t")
    score_cols = [c for c in df.columns if c.endswith("_score")]

    # Per-donor composition + scores
    g = df.groupby("donor_id")
    # Injury-program fraction: pericytes in stable clusters annotated as an
    # injury program (state_program), donor-normalized.
    df["_is_injury"] = df["state_program"].isin(INJURY_PROGRAMS).astype(float)
    donor = pd.DataFrame({"n_cells": g.size()})
    donor["injury_frac"] = g["_is_injury"].mean()
    donor["AGTR1_pos_frac"] = g["AGTR1_detect"].mean()
    # Donor means of the curated module scores (used to build the index).
    for c in score_cols:
        donor[f"mean_{c}"] = g[c].mean()
    # carry covariates
    meta = g.agg(disease=("disease", "first"),
                 lung_condition=("lung_condition", "first"),
                 smoking_status=("smoking_status", "first"),
                 sex=("sex", "first"),
                 ethnicity=("self_reported_ethnicity", "first"),
                 age=("age_or_mean_of_age_range", "mean"),
                 BMI=("BMI", "first"))
    donor = donor.join(meta)
    donor = donor[donor["n_cells"] >= args.min_cells].copy()

    # Airspace proximity per donor
    air = pd.read_csv(args.airspace_summary)
    air = air[["donor_id", "mean_airspace_score"]].dropna(subset=["mean_airspace_score"])
    donor = donor.merge(air, on="donor_id", how="left")

    # ---- Build composite scores from z-scored components --------------------
    # Curated module-score components (face value; states already integration-defined).
    stab_components = ["mean_vascular_stabilizing_score", "mean_airspace_score"]
    # PRIMARY injury composite: no AGTR1+ fraction (avoid circularity w/ focal receptor).
    inj_components = ["mean_inflammatory_score", "mean_fibroblast_like_score",
                      "mean_activated_migratory_score", "injury_frac"]
    # SENSITIVITY-only injury composite adds the AGTR1+ fraction (supplement).
    inj_components_sens = inj_components + ["AGTR1_pos_frac"]
    stab_components = [c for c in stab_components if c in donor.columns]
    inj_components = [c for c in inj_components if c in donor.columns]
    inj_components_sens = [c for c in inj_components_sens if c in donor.columns]

    for c in set(stab_components + inj_components + inj_components_sens):
        donor[f"z_{c}"] = zscore(donor[c])

    donor["niche_stability_score"] = donor[[f"z_{c}" for c in stab_components]].mean(axis=1)
    # primary (no AGTR1)
    donor["injury_stromal_score"] = donor[[f"z_{c}" for c in inj_components]].mean(axis=1)
    donor["niche_index"] = donor["niche_stability_score"] - donor["injury_stromal_score"]
    # sensitivity (with AGTR1+ fraction) -- supplement only
    donor["injury_stromal_score_sens_agtr1"] = donor[[f"z_{c}" for c in inj_components_sens]].mean(axis=1)
    donor["niche_index_sens_agtr1"] = donor["niche_stability_score"] - donor["injury_stromal_score_sens_agtr1"]

    donor.to_csv(args.outdir / "niche_index_per_donor.tsv.gz", sep="\t", index=False)
    # Record component definitions for transparency
    with open(args.outdir / "niche_index_components.txt", "w") as fh:
        fh.write("stability_components: " + ", ".join(stab_components) + "\n")
        fh.write("injury_components (PRIMARY): " + ", ".join(inj_components) + "\n")
        fh.write("injury_components_sens_agtr1 (SUPPLEMENT only): "
                 + ", ".join(inj_components_sens) + "\n")
    logging.info(f"Wrote niche index for {donor.shape[0]} donors")


if __name__ == "__main__":
    main()
