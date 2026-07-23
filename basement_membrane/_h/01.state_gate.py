"""
Pre-specified gate: would a basement-membrane panel change the pericyte state model?

`state_program` is assigned per CLUSTER by a relative-enrichment argmax over the
curated panels (pericyte_states/_h/00.state_discovery.py:279-300). Leiden
clustering itself runs on X_pca_harmony and is independent of the marker panels,
so adding a sixth panel cannot move the clusters -- it can only change which
program each cluster is labelled with. The gate is therefore deterministic and
cheap: recompute the argmax over 6 panels instead of 5 against the frozen
`pericyte_state` labels, and compare.

Decision rule, fixed in advance:
  metric A = % of cells whose state_program changes
  metric B = whether any cluster flips its dominant program
  Because programs are assigned per cluster, B == none implies A == 0; the two
  are not independent. Operative rule: escalate to the 6-panel model as canonical
  (and re-run pericyte_cogaps / cell_communication / niche_index /
  pathway_balance / figures) if any cluster flips AND the affected cells are >=5%
  of the total. A flip affecting <5% keeps the frozen model canonical, reports BM
  additively, and surfaces the alternative model for a human decision rather than
  silently choosing.

VALIDATION: the 5-panel arm recomputed here must reproduce the committed
pericyte_states/_m/annotations/state_program_map.tsv exactly. If it does not, the
gate logic has diverged from annotate_states and its 6-panel result cannot be
trusted -- the script fails rather than reporting a number built on drifted logic.

Outputs:
  - state_gate_summary.tsv                  (metrics A and B, escalation verdict)
  - state_gate_crosstab.tsv                 (cluster x program, both models)
  - state_gate_relenrich.tsv                (per-cluster relative enrichment, both)
  - annotations/state_program_map_6panel.tsv
"""
import numpy as np
import pandas as pd
import session_info
import logging, argparse, sys
from pathlib import Path

# Panel order must match STATE_PANELS in pericyte_states/_h/00.state_discovery.py:
# the argmax breaks ties by first column, so ordering is part of the logic.
BASE_PROGRAMS = ["vascular_stabilizing", "inflammatory", "synthetic_contractile",
                 "activated_migratory", "fibroblast_like"]
BM_PROGRAM = "basement_membrane"
ESCALATE_PCT = 5.0


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--state-meta", required=True, type=Path,
                   help="pericyte_states/_m/pericytes_states_metadata.tsv.gz")
    p.add_argument("--bm-meta", required=True, type=Path,
                   help="bm_metadata.tsv.gz from 00.bm_score.py")
    p.add_argument("--reference-map", required=True, type=Path,
                   help="pericyte_states/_m/annotations/state_program_map.tsv")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--escalate-pct", type=float, default=ESCALATE_PCT)
    return p.parse_args()


def relative_enrichment(df, programs):
    """Reproduces annotate_states: z-score each program across cells, average
    within cluster, label each cluster by its most-enriched program."""
    score_cols = [f"{p}_score" for p in programs]
    z = df[score_cols].apply(lambda s: (s - s.mean()) / (s.std(ddof=0) + 1e-9))
    z["pericyte_state"] = df["pericyte_state"].to_numpy()
    rel = z.groupby("pericyte_state", observed=True)[score_cols].mean()
    rel.columns = [f"{p}_relenrich" for p in programs]
    dominant = rel.to_numpy().argmax(axis=1)
    prog = pd.Series([programs[i] for i in dominant], index=rel.index,
                     name="state_program")
    return prog, rel


def main():
    configure_logging()
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "annotations").mkdir(parents=True, exist_ok=True)

    state = pd.read_csv(args.state_meta, sep="\t")
    bm = pd.read_csv(args.bm_meta, sep="\t")
    key = "index"
    if key not in state.columns:
        state = state.rename(columns={state.columns[0]: key})
    # Bring across every BM panel score: the primary gate needs only
    # basement_membrane, but the robustness variants below need the sub-panels
    # and the leave-COL4A1-out fibroblast panel.
    bm_scores = [c for c in bm.columns if c.endswith("_score")]
    merged = state.merge(bm[[key] + bm_scores], on=key, how="inner",
                         validate="one_to_one")
    if len(merged) != len(state):
        raise ValueError(f"join lost cells: {len(state)} -> {len(merged)}")
    logging.info(f"Merged {len(merged)} cells")

    prog5, rel5 = relative_enrichment(merged, BASE_PROGRAMS)
    prog6, rel6 = relative_enrichment(merged, BASE_PROGRAMS + [BM_PROGRAM])

    # --- validation: the 5-panel arm must match the committed canonical map ---
    ref = pd.read_csv(args.reference_map, sep="\t")
    ref_prog = (ref.set_index("pericyte_state")["state_program"]
                .astype(str))
    got = prog5.copy()
    got.index = got.index.astype(ref_prog.index.dtype)
    ref_prog = ref_prog.reindex(got.index)
    if not (got.astype(str).values == ref_prog.values).all():
        cmp = pd.DataFrame({"recomputed": got.astype(str), "committed": ref_prog})
        logging.error("5-panel arm does NOT reproduce the canonical map:\n"
                      + cmp.to_string())
        sys.exit("Gate logic has diverged from annotate_states; refusing to "
                 "report a 6-panel result built on it.")
    logging.info("Validation passed: 5-panel arm reproduces the canonical map.")

    # --- metrics ---
    sizes = merged["pericyte_state"].value_counts()
    sizes.index = sizes.index.astype(prog5.index.dtype)
    flipped = prog5.index[prog5.values != prog6.reindex(prog5.index).values]
    n_flip_cells = int(sizes.reindex(flipped).sum()) if len(flipped) else 0
    pct = 100.0 * n_flip_cells / len(merged)
    escalate = bool(len(flipped) > 0 and pct >= args.escalate_pct)

    crosstab = pd.DataFrame({
        "pericyte_state": prog5.index,
        "n_cells": sizes.reindex(prog5.index).to_numpy(),
        "program_5panel": prog5.to_numpy(),
        "program_6panel": prog6.reindex(prog5.index).to_numpy(),
    })
    crosstab["flipped"] = crosstab["program_5panel"] != crosstab["program_6panel"]
    crosstab.to_csv(args.outdir / "state_gate_crosstab.tsv", sep="\t", index=False)

    rel = rel5.join(rel6[[f"{BM_PROGRAM}_relenrich"]], how="left")
    rel.reset_index().to_csv(args.outdir / "state_gate_relenrich.tsv",
                             sep="\t", index=False)

    prog6.reset_index().join(
        rel6.reset_index(drop=True)).to_csv(
        args.outdir / "annotations" / "state_program_map_6panel.tsv",
        sep="\t", index=False)

    # --- robustness: is a flip an artifact of COL4A1 sitting in BOTH panels? ---
    # COL4A1 is a member of the frozen `fibroblast_like` panel and of the BM panel,
    # so a naive comparison double-counts it. Each variant below re-runs the argmax
    # with that shared gene removed, or with only part of the BM panel, to show
    # whether the verdict rests on one gene or on the whole programme.
    variants = {
        "fibroblast_like_without_COL4A1": (
            [p for p in BASE_PROGRAMS if p != "fibroblast_like"]
            + ["fibroblast_like_noCOL4A1", BM_PROGRAM]),
        "bm_laminin_only": BASE_PROGRAMS + ["bm_laminin"],
        "bm_linker_only": BASE_PROGRAMS + ["bm_linker"],
        "bm_collagen_iv_only": BASE_PROGRAMS + ["bm_collagen_iv"],
    }
    rob_rows = []
    for name, progs in variants.items():
        missing = [f"{p}_score" for p in progs
                   if f"{p}_score" not in merged.columns]
        if missing:
            logging.warning(f"[{name}] skipped; missing {missing}")
            continue
        vprog, _ = relative_enrichment(merged, progs)
        for cluster in prog5.index:
            rob_rows.append({
                "variant": name, "pericyte_state": cluster,
                "n_cells": int(sizes.get(cluster, 0)),
                "program_frozen": prog5[cluster],
                "program_variant": vprog[cluster],
                "flipped": prog5[cluster] != vprog[cluster],
            })
    robust = pd.DataFrame(rob_rows)
    robust.to_csv(args.outdir / "state_gate_robustness.tsv", sep="\t", index=False)

    # Correlation of the BM axis with the fibrillar axes: if BM merely restates
    # the fibroblast-like axis, these are high and the flip is uninformative.
    corr_rows = []
    for other in ["fibroblast_like", "fibroblast_like_noCOL4A1", "fibrillar_ecm",
                  "vascular_stabilizing"]:
        col = f"{other}_score"
        if col in merged.columns:
            corr_rows.append({
                "panel_a": BM_PROGRAM, "panel_b": other,
                "pearson_r": float(np.corrcoef(merged[f"{BM_PROGRAM}_score"],
                                               merged[col])[0, 1]),
            })
    pd.DataFrame(corr_rows).to_csv(args.outdir / "state_gate_axis_correlation.tsv",
                                   sep="\t", index=False)
    if len(robust):
        logging.info("\nRobustness (does the verdict survive each variant?):\n"
                     + robust.to_string(index=False))
    logging.info("\nBM axis correlations:\n"
                 + pd.DataFrame(corr_rows).to_string(index=False))

    summary = pd.DataFrame([{
        "n_cells": len(merged),
        "n_clusters": len(prog5),
        "metric_A_pct_cells_changed": round(pct, 4),
        "metric_B_n_clusters_flipped": len(flipped),
        "clusters_flipped": ",".join(map(str, flipped)) if len(flipped) else "",
        "escalate_threshold_pct": args.escalate_pct,
        "escalate": escalate,
        "verdict": ("ESCALATE: adopt 6-panel model as canonical and re-run "
                    "downstream modules" if escalate else
                    "HOLD: frozen 5-panel model stays canonical; BM reported "
                    "additively"),
    }])
    summary.to_csv(args.outdir / "state_gate_summary.tsv", sep="\t", index=False)

    logging.info("\n" + crosstab.to_string(index=False))
    logging.info(f"metric A (% cells changed): {pct:.3f}%")
    logging.info(f"metric B (clusters flipped): {len(flipped)}"
                 + (f" -> {list(flipped)}" if len(flipped) else ""))
    logging.info(summary["verdict"].iloc[0])
    if len(flipped) and not escalate:
        logging.warning(
            f"{len(flipped)} cluster(s) flipped but only {pct:.2f}% of cells "
            f"(<{args.escalate_pct}%). Frozen model held; the 6-panel map is "
            "written for review -- this is a decision point, not a silent pass.")

    session_info.show()


if __name__ == "__main__":
    main()
