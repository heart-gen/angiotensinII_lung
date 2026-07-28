"""
Donor -> study/dataset + lung_condition map for the stroma HLCA subset.

`01.disease_association.R` builds `mean_expr/donor_metadata.tsv` (donor x cell
type AGTR1 summaries) but keeps only the columns it needed for a Kruskal-Wallis
test, so it carries no `dataset`. The three-group AGTR1 model in
`05.agtr1_celltype_disease.R` treats STUDY as a modelled factor exactly as
`03.disease_forest.R` does, which needs that column back.

Reads obs only (no X), so this is cheap despite the 5 GB file.
"""
import h5py
import argparse
import numpy as np
import pandas as pd
import session_info
from pathlib import Path

# Columns pulled out of obs. `patient` is the donor key `01.disease_association.R`
# aggregates on; `donor_id` is kept alongside because the pericyte modules key on
# it and the two are not always identical.
WANT = ["patient", "donor_id", "dataset", "study", "lung_condition", "disease"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outfile", required=True, type=Path)
    return p.parse_args()


def read_obs_col(obs, key):
    """Decode an h5ad obs column, categorical or plain array."""
    node = obs[key]
    if isinstance(node, h5py.Group):          # categorical
        cats = node["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        codes = node["codes"][:]
        out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        return out
    arr = node[:]
    if arr.dtype.kind in "SO":
        arr = np.array([a.decode() if isinstance(a, bytes) else a for a in arr])
    return arr


def main():
    args = parse_args()
    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.adata, "r") as f:
        obs = f["obs"]
        present = [k for k in WANT if k in obs.keys()]
        missing = sorted(set(WANT) - set(present))
        if missing:
            print(f"absent from obs (skipped): {missing}")
        df = pd.DataFrame({k: read_obs_col(obs, k) for k in present})

    # One row per donor. Most HLCA donors sit in exactly one `dataset`, but a
    # handful were sequenced by two sub-datasets of the same study (e.g. the
    # 3'/5' Meyer_2021 arms). Resolve those to the dataset contributing the most
    # cells rather than dropping them, and carry `n_datasets` so a downstream
    # model can see which donors were collapsed.
    n_ds = df.groupby("patient")["dataset"].nunique().rename("n_datasets")
    if (n_ds > 1).any():
        print(f"{int((n_ds > 1).sum())} donors span >1 dataset; "
              "assigning each its majority (most cells) dataset")

    modal = (df.groupby(["patient", "dataset"]).size()
             .rename("n_cells").reset_index()
             .sort_values(["patient", "n_cells"], ascending=[True, False])
             .drop_duplicates("patient")[["patient", "dataset", "n_cells"]])
    other = (df.drop(columns=["dataset"])
             .drop_duplicates(subset="patient"))
    out = (other.merge(modal, on="patient", how="left")
           .merge(n_ds, left_on="patient", right_index=True, how="left")
           .reset_index(drop=True))
    out.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote {args.outfile}: {out.shape[0]} donors, "
          f"{out['dataset'].nunique()} datasets")
    print(out["lung_condition"].value_counts())
    session_info.show()


if __name__ == "__main__":
    main()
