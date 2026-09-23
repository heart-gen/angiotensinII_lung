"""
Figure 5D -- score every lung pericyte for the AT1R/AngII response.

PRIMARY: the frozen perturbation-derived signature (_h/signatures/, derived by
00b/00c from an in-vivo AngII-vs-vehicle single-cell experiment). Scored with
sc.tl.score_genes; `at1r_response_score` = up-score minus down-score (up only if
fewer than 5 down genes survive).
SECONDARY: decoupler ULM activities on cached networks (00.fetch_networks.py):
PROGENy MAPK, NFkB, TGFb, EGFR, JAK-STAT (all 14 pathways written, 5
pre-specified) and CollecTRI JUN, FOS, NFKB1, RELA, SMAD2, SMAD3. These are
response programmes, not proof of receptor activation.
LABELLED COMPARISON ONLY: pathway_balance's curated AT1R_PROGRAM, which its own
docstring flags as non-specific (it re-measures injury intensity).

DISJOINTNESS (README step 8). A response score that shares genes with the state,
matrix, tracer, TGF-beta or RAS panels would make every downstream association
partly arithmetic. Such genes are PRUNED from the signature and every pruning is
logged; >= --min-genes must remain or the script stops. AGTR1 itself is in the
RAS set, so the score can never contain the exposure it is tested against.

State panels are read by parsing pericyte_states/_h/00.state_discovery.py with
`ast` rather than importing it: importing would execute its module-level setup.
"""
import argparse
import ast
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "basement_membrane" / "_h"))
import bm_panels  # noqa: E402

RAS_GENES = {"AGT", "REN", "ACE", "ACE2", "CMA1", "CTSG", "CTSD", "ENPEP", "MME",
             "AGTR1", "AGTR2", "LRP2", "MAS1"}
PROGENY_PRESPEC = ["MAPK", "NFkB", "TGFb", "EGFR", "JAK-STAT"]
TF_PRESPEC = ["JUN", "FOS", "NFKB1", "RELA", "SMAD2", "SMAD3"]
AT1R_PROGRAM = [  # pathway_balance/_h/00.pathway_balance.py, verbatim
    "TGFB1", "CCN2", "CTGF", "CCN1", "CYR61", "SERPINE1", "COL1A1", "COL3A1",
    "FN1", "ACTA2", "EDN1", "NOX4", "IL6", "CCL2", "NFKB1", "MYC", "EGR1",
    "FOS", "AGT", "ACE", "PLAU"]


def state_panels():
    src = (ROOT / "pericyte_states" / "_h" / "00.state_discovery.py").read_text()
    tree = ast.parse(src)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "STATE_PANELS" for t in node.targets):
            out = {}
            for k, v in zip(node.value.keys, node.value.values):
                key = ast.literal_eval(k)
                try:
                    out[key] = list(ast.literal_eval(v))
                except ValueError:
                    if key != "basement_membrane":
                        raise
                    out[key] = list(bm_panels.BM_PANEL)
            return out
    raise KeyError("STATE_PANELS not found in 00.state_discovery.py")


def excluded():
    ex = {}
    for name, genes in state_panels().items():
        for g in genes:
            ex.setdefault(g, f"pericyte_states.{name}")
    for name, genes in bm_panels.PANELS.items():
        for g in genes:
            ex.setdefault(g, f"bm_panels.{name}")
    for g in RAS_GENES:
        ex.setdefault(g, "ras_panel")
    return ex


def load_null_module():
    spec = importlib.util.spec_from_file_location(
        "tgfb_null", ROOT / "basement_membrane" / "_h" / "11.tgfb_null_panels.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", type=Path, required=True)
    p.add_argument("--signature", type=Path, required=True)
    p.add_argument("--networks", type=Path, required=True)
    p.add_argument("--continuum", type=Path, required=True)
    p.add_argument("--bm", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--min-genes", type=int, default=15)
    p.add_argument("--seed", type=int, default=13)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    (a.outdir / "stats_data").mkdir(parents=True, exist_ok=True)

    nm = load_null_module()
    adata = nm.load_anndata(a.adata)
    det = nm.detection_frac(adata)
    logging.info(f"pericytes: {adata.n_obs} x {adata.n_vars}")

    sig = pd.read_csv(a.signature, sep="\t")
    ex = excluded()
    audit = []
    for _, r in sig.iterrows():
        g = r["gene"]
        present = g in adata.var_names
        reason = ex.get(g, "")
        audit.append({"gene": g, "direction": r["direction"],
                      "log2FoldChange": r.get("log2FoldChange", np.nan),
                      "present_in_object": present,
                      "pericyte_detect": float(det[g]) if present else np.nan,
                      "pruned_reason": reason if reason else ("absent" if not present else ""),
                      "kept": present and not reason,
                      "in_legacy_AT1R_PROGRAM": g in AT1R_PROGRAM})
    audit = pd.DataFrame(audit)
    audit.to_csv(a.outdir / "stats_data" / "at1r_signature_audit.tsv", sep="\t", index=False)
    up = audit.query("kept and direction == 'up'")["gene"].tolist()
    dn = audit.query("kept and direction == 'down'")["gene"].tolist()
    logging.info(f"signature: {len(sig)} genes -> kept {len(up)} up, {len(dn)} down; "
                 f"pruned {int((audit['pruned_reason'] != '').sum())}")
    if len(up) + len(dn) < a.min_genes:
        raise SystemExit(f"only {len(up) + len(dn)} signature genes survive pruning "
                         f"(< {a.min_genes}); README step 6 says fall back to the next "
                         "cell type -- delete the frozen file and re-derive.")
    if len(up) < 5:
        raise SystemExit("fewer than 5 up-regulated genes survive; the score is undefined")

    sc.tl.score_genes(adata, up, score_name="at1r_up", random_state=a.seed, use_raw=False)
    if len(dn) >= 5:
        sc.tl.score_genes(adata, dn, score_name="at1r_down", random_state=a.seed,
                          use_raw=False)
        adata.obs["at1r_response_score"] = adata.obs["at1r_up"] - adata.obs["at1r_down"]
        mode = "up_minus_down"
    else:
        adata.obs["at1r_down"] = np.nan
        adata.obs["at1r_response_score"] = adata.obs["at1r_up"]
        mode = "up_only"
    leg = [g for g in AT1R_PROGRAM if g in adata.var_names]
    sc.tl.score_genes(adata, leg, score_name="at1r_legacy_curated_score",
                      random_state=a.seed, use_raw=False)

    # ---- decoupler (cached networks only) ------------------------------------
    import decoupler as dc
    acts = {}
    for fname, prefix, keep in (("progeny_human_top500.tsv", "progeny_", None),
                                ("collectri_human.tsv", "tf_", TF_PRESPEC)):
        f = a.networks / fname
        if not f.exists():
            raise FileNotFoundError(f"{f} -- run 00.fetch_networks.py first")
        net = pd.read_csv(f, sep="\t")
        if keep is not None:
            net = net[net["source"].isin(keep)]
        dc.mt.ulm(data=adata, net=net, tmin=5)
        sc_ = adata.obsm["score_ulm"]
        sc_ = sc_ if isinstance(sc_, pd.DataFrame) else pd.DataFrame(sc_, index=adata.obs_names)
        for c in sc_.columns:
            acts[prefix + str(c)] = sc_[c].to_numpy()
        missing = sorted(set(keep or PROGENY_PRESPEC) - set(map(str, sc_.columns)))
        if missing:
            logging.warning(f"{fname}: pre-specified regulators pruned by tmin: {missing}")
        del adata.obsm["score_ulm"]
        if "padj_ulm" in adata.obsm:
            del adata.obsm["padj_ulm"]
    acts = pd.DataFrame(acts, index=adata.obs_names)
    if {"tf_JUN", "tf_FOS"} <= set(acts.columns):
        acts["tf_AP1"] = acts[["tf_JUN", "tf_FOS"]].mean(axis=1)

    obs = adata.obs[["donor_id", "pericyte_state", "at1r_up", "at1r_down",
                     "at1r_response_score", "at1r_legacy_curated_score"]].copy()
    obs = pd.concat([obs, acts], axis=1)
    obs.index.name = "index"
    cont = pd.read_csv(a.continuum, sep="\t", index_col=0)
    bm = pd.read_csv(a.bm, sep="\t", index_col=0)
    cc = [c for c in ["dpt_pseudotime", "AGTR1_expr", "AGTR1_scvi",
                      "vascular_stabilizing_score", "inflammatory_score",
                      "synthetic_contractile_score", "activated_migratory_score",
                      "fibroblast_like_score"] if c in cont.columns]
    bc = [c for c in ["study", "dataset", "state_program", "log10_total_counts",
                      "total_counts", "basement_membrane_score", "fibrillar_ecm_score",
                      "fibrillar_collagen_score", "ambient_tracer_score"] if c in bm.columns]
    out = obs.join(cont[cc], how="left").join(bm[bc], how="left")
    miss = int(out["dpt_pseudotime"].isna().sum() + out["basement_membrane_score"].isna().sum())
    if len(out) != adata.n_obs or miss:
        raise AssertionError(f"barcode join lost cells: rows {len(out)} vs {adata.n_obs}, "
                             f"{miss} unmatched fields")
    out["score_mode"] = mode
    out.to_csv(a.outdir / "at1r_response_metadata.tsv.gz", sep="\t")
    logging.info(f"wrote {len(out)} pericytes; score mode {mode}; activities "
                 f"{[c for c in acts.columns]}")


if __name__ == "__main__":
    main()
