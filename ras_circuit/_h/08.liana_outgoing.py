"""
Figure 5F -- what do pericytes SIGNAL TO? (outgoing LIANA edges)

cell_communication/_h/01.run_liana.py already scored every source x target pair
(rank_aggregate with no source/target restriction), but exported only edges INTO
the pericyte/AT2 receivers. This script re-filters the existing per-disease
result tables -- no LIANA re-run -- to edges whose SOURCE is Pericytes, into the
alveolar-capillary neighbours, and puts the matching incoming edges beside them
so panel F can show both directions on one scale.

Flags, not deletions:
  - `not_a_ligand`: COPA, MMP14, SIRPB2 (the same three prior-network entries
    Figures 4, S8 and S9A exclude; see figures/_h/_fig_common.R).
  - `substrate_not_ligand`: AGT. Angiotensinogen is renin's substrate, not the
    AT1R ligand; an AGT -> AGTR1 "edge" is the direct edge the Figure 5 DAG
    exists to forbid (agt_axis/_h/01.ras_landscape_stats.R, header).
AT2 is split upstream on AGTR2 detectability (a receiver definition, not
biology here), so an `AT2 (either stratum)` row is added per edge carrying the
BETTER-supported of the two strata; both source rows are kept.

Outputs (to --outdir):
  liana_pericyte_outgoing.tsv.gz   pericyte -> neighbour edges, all diseases
  liana_pericyte_incoming.tsv.gz   neighbour -> pericyte edges, same targets
  liana_direction_summary.tsv      supported-edge counts by direction/partner/disease
"""
import argparse
import logging
from pathlib import Path

import pandas as pd

NON_LIGANDS = {"COPA", "MMP14", "SIRPB2"}
SUBSTRATES = {"AGT"}
PARTNERS = ["EC aerocyte capillary", "EC general capillary", "AT1",
            "AT2_AGTR2det", "AT2_AGTR2undet", "Alveolar fibroblasts",
            "Adventitial fibroblasts", "Alveolar macrophages",
            "Interstitial macrophages", "Vascular smooth muscle"]
DISEASES = ["Healthy", "COPD", "Fibrotic_ILD", "Other"]


def flag(df):
    lig = df["ligand_complex"].astype(str)
    df["not_a_ligand"] = lig.isin(NON_LIGANDS)
    df["substrate_not_ligand"] = lig.isin(SUBSTRATES)
    df["displayable"] = ~(df["not_a_ligand"] | df["substrate_not_ligand"])
    return df


def add_at2_either(df, partner_col):
    at2 = df[df[partner_col].isin(["AT2_AGTR2det", "AT2_AGTR2undet"])]
    if at2.empty:
        return df
    keys = ["source", "target", "ligand_complex", "receptor_complex", "disease_group"]
    keys = [k for k in keys if k != partner_col]
    best = (at2.sort_values("magnitude_rank")
               .groupby(keys, as_index=False, observed=True).first())
    best[partner_col] = "AT2 (either stratum)"
    best["at2_recombined"] = True
    df = df.assign(at2_recombined=False)
    return pd.concat([df, best], ignore_index=True)


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--liana-dir", type=Path, required=True,
                   help="cell_communication/_m")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--magnitude-thr", type=float, default=0.05,
                   help="magnitude_rank at or below which an edge is 'supported'")
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    a.outdir.mkdir(parents=True, exist_ok=True)

    frames = []
    for dz in DISEASES:
        f = a.liana_dir / f"liana_res_main_{dz}.tsv.gz"
        if not f.exists():
            logging.warning(f"{f.name} absent -- {dz} not produced upstream; skipped")
            continue
        d = pd.read_csv(f, sep="\t")
        if "disease_group" not in d.columns:
            d["disease_group"] = dz
        frames.append(d)
        logging.info(f"{dz}: {len(d)} edges, sources include Pericytes: "
                     f"{(d['source'] == 'Pericytes').any()}")
    if not frames:
        raise FileNotFoundError("no liana_res_main_*.tsv.gz under " + str(a.liana_dir))
    res = pd.concat(frames, ignore_index=True)
    if not (res["source"] == "Pericytes").any():
        raise AssertionError("no edge has source == 'Pericytes'; the upstream run "
                             "did not score pericytes as senders")

    out = res[(res["source"] == "Pericytes") & res["target"].isin(PARTNERS)].copy()
    inc = res[(res["target"] == "Pericytes") & res["source"].isin(PARTNERS)].copy()
    out = add_at2_either(flag(out), "target")
    inc = add_at2_either(flag(inc), "source")
    out["direction"] = "outgoing"
    inc["direction"] = "incoming"
    out["supported"] = out["magnitude_rank"] <= a.magnitude_thr
    inc["supported"] = inc["magnitude_rank"] <= a.magnitude_thr
    out.sort_values(["disease_group", "target", "magnitude_rank"]).to_csv(
        a.outdir / "liana_pericyte_outgoing.tsv.gz", sep="\t", index=False)
    inc.sort_values(["disease_group", "source", "magnitude_rank"]).to_csv(
        a.outdir / "liana_pericyte_incoming.tsv.gz", sep="\t", index=False)

    rows = []
    for dirn, df, pcol in (("outgoing", out, "target"), ("incoming", inc, "source")):
        for (dz, partner), g in df.groupby(["disease_group", pcol], observed=True):
            gg = g[g["displayable"]]
            rows.append({"direction": dirn, "disease_group": dz, "partner": partner,
                         "n_edges": len(gg), "n_supported": int(gg["supported"].sum()),
                         "n_flagged": int((~g["displayable"]).sum()),
                         "magnitude_thr": a.magnitude_thr})
    summ = pd.DataFrame(rows)
    summ.to_csv(a.outdir / "liana_direction_summary.tsv", sep="\t", index=False)
    logging.info("\n" + summ[summ["disease_group"] == "Healthy"].to_string(index=False))


if __name__ == "__main__":
    main()
