"""
Supplementary Table S3 -- curated pericyte program gene sets.

Only the basement-membrane panel was ever machine-readable
(basement_membrane/_m/bm_panel_genes.tsv). The other five programs existed solely
as the STATE_PANELS dict inside pericyte_states/_h/00.state_discovery.py, and no
table anywhere recorded which panel genes are actually present in the analysis
object -- which is how the activated_migratory panel came to be described as 7
genes in writings/STATES_SUMMARY.md while the script lists 8.

The panels are IMPORTED from their defining modules rather than copied, so this
table cannot silently fork from the code that built the scores.

Outputs: tsv/tableS03.tsv
"""
import argparse
import importlib.util
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import session_info


def load_module(path: Path, name: str):
    """Import a module by file path (module names here start with digits)."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def register_part(registry: Path, part: str, df: pd.DataFrame, file: str,
                  title: str, supports: str, sources: str, notes: str,
                  status: str = "complete"):
    """Append this part to the shared manifest that 08.assemble_tables.R reads.

    Mirrors write_part() in _tab_common.R -- same columns, same replace-on-rerun
    behaviour. Without this the Python-built table would exist on disk but never
    reach the workbook, because the manifest is the sole authority for numbering.
    """
    row = pd.DataFrame([{
        "part": part, "title": title, "supports": supports,
        "n_rows": len(df), "n_cols": df.shape[1], "status": status,
        "notes": notes, "file": file, "sources": sources,
        "built": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}])
    if registry.exists():
        old = pd.read_csv(registry, sep="\t", dtype=str)
        old = old[old["part"] != part]
        row = pd.concat([old, row], ignore_index=True)
    row = row.assign(_n=row["part"].str.len()).sort_values(["_n", "part"]) \
             .drop(columns="_n")
    row.to_csv(registry, sep="\t", index=False, na_rep="")
    logging.info(f"registered part S{part} in {registry}")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--root", type=Path, default=Path("../.."),
                   help="Repository root")
    p.add_argument("--adata", type=Path, default=None,
                   help="pericyte_states.h5ad (default: <root>/pericyte_states/_m/...)")
    p.add_argument("--outdir", type=Path, default=Path("./tsv"))
    return p.parse_args()


# Biological rationale per program. The BM subclasses come from the comments in
# basement_membrane/_h/bm_panels.py; the rest are the categories the panels were
# curated under. Free text nowhere else in the repo is machine-readable.
PROGRAM_RATIONALE = {
    "vascular_stabilizing": "Canonical pericyte identity and vessel-stabilizing "
                            "signalling (PDGFRB/NOTCH3 axis, KATP subunits, "
                            "mitochondrial COX4I2/NDUFA4L2)",
    "inflammatory": "NF-kB-driven cytokine, chemokine and leukocyte-adhesion "
                    "program",
    "synthetic_contractile": "Smooth-muscle-like contractile apparatus "
                             "(actin/myosin/calponin)",
    "activated_migratory": "Matrix-remodelling proteases, their inhibitors and "
                           "matricellular activation markers",
    "fibroblast_like": "Fibrillar interstitial matrix and fibroblast identity "
                       "markers",
    "basement_membrane": "Vascular basement-membrane deposition: collagen IV "
                         "network, laminins, and nidogen/perlecan/agrin linkers",
}

GENE_CATEGORY = {  # structural subclass where the panel defines one
    **{g: "collagen IV" for g in ["COL4A1", "COL4A2"]},
    **{g: "laminin" for g in ["LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMC1"]},
    **{g: "linker/proteoglycan" for g in ["NID1", "NID2", "HSPG2", "AGRN"]},
    "COL18A1": "collagen XVIII (BM-associated)",
}


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    root = args.root.resolve()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    sd = load_module(root / "pericyte_states" / "_h" / "00.state_discovery.py",
                     "state_discovery")
    bm = load_module(root / "basement_membrane" / "_h" / "bm_panels.py", "bm_panels")

    rows = []
    for program, genes in sd.STATE_PANELS.items():
        for g in genes:
            rows.append({"panel_set": "state_program", "program": program, "gene": g})
    for g in sd.EXTRA_GENES:
        rows.append({"panel_set": "extra", "program": "receptor/control genes",
                     "gene": g})

    # Basement-membrane subpanels and contrast panels, from their own module.
    for name, genes in bm.BM_SUBPANELS.items():
        for g in genes:
            rows.append({"panel_set": "bm_subpanel", "program": name, "gene": g})
    # `bm_v1` is the frozen 13-gene BM panel every published BM number was
    # computed on; it is listed so a reader can see exactly which seven genes the
    # 2026-09-01 expansion added. The TGF-beta panels are listed because the
    # matrix-vs-TGF-beta association is only interpretable if the reader can
    # verify the response signature shares no gene with the matrix panels.
    for name, genes in [("fibrillar_ecm", bm.FIBRILLAR_CONTRAST),
                        ("fibrillar_collagen", bm.FIBRILLAR_COLLAGEN),
                        ("bm_v1_frozen_13gene", bm.BM_PANEL_V1),
                        ("tgfb_response", bm.TGFB_RESPONSE),
                        ("tgfb_receptor", bm.TGFB_RECEPTOR),
                        ("fibroblast_like_noCOL4A1", bm.FIBROBLAST_LIKE_NO_COL4A1)]:
        for g in genes:
            rows.append({"panel_set": "contrast_panel", "program": name, "gene": g})

    df = pd.DataFrame(rows)

    # Multi-panel membership, counted over the PRIMARY state programs only --
    # subpanels are by construction subsets of basement_membrane and would inflate
    # the count into meaninglessness.
    primary = df[df["panel_set"] == "state_program"]
    n_programs = primary.groupby("gene")["program"].nunique()
    programs_of = (primary.groupby("gene")["program"]
                   .apply(lambda s: "; ".join(sorted(set(s)))))
    df["n_state_programs"] = df["gene"].map(n_programs).fillna(0).astype(int)
    df["state_programs"] = df["gene"].map(programs_of).fillna("")
    df["multi_panel"] = df["n_state_programs"] > 1

    df["program_rationale"] = df["program"].map(PROGRAM_RATIONALE).fillna("")
    df["gene_category"] = df["gene"].map(GENE_CATEGORY).fillna("")

    # ---- detection in the analysis object -------------------------------
    adata_path = args.adata or (root / "pericyte_states" / "_m" / "pericyte_states.h5ad")
    if adata_path.exists():
        backed = ad.read_h5ad(adata_path, backed="r")
        var_names = set(map(str, backed.var_names))
        wanted = sorted(set(df["gene"]) & var_names)
        logging.info(f"genes present in object: {len(wanted)}/{df['gene'].nunique()}")

        idx = [backed.var_names.get_loc(g) for g in wanted]
        n = backed.n_obs
        det = np.zeros(len(wanted))
        tot = np.zeros(len(wanted))
        chunk = 50000
        for start in range(0, n, chunk):
            stop = min(start + chunk, n)
            X = backed[start:stop].to_memory().X
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            sub = X[:, idx]
            det += (sub > 0).sum(axis=0)
            tot += sub.sum(axis=0)
        stats = pd.DataFrame({"gene": wanted,
                              "detect_frac": det / n,
                              "mean_expr": tot / n})
        df = df.merge(stats, on="gene", how="left")
        df["detected_in_object"] = df["gene"].isin(var_names)
        df["n_cells_object"] = n
    else:
        logging.warning(f"{adata_path} not found -- detection columns left empty")
        df["detect_frac"] = np.nan
        df["mean_expr"] = np.nan
        df["detected_in_object"] = pd.NA
        df["n_cells_object"] = pd.NA

    df = df.sort_values(["panel_set", "program", "gene"]).reset_index(drop=True)
    out = outdir / "tableS03.tsv"
    df.to_csv(out, sep="\t", index=False)
    logging.info(f"wrote {out} ({len(df)} rows)")

    n_absent = int((df["detected_in_object"] == False).sum())  # noqa: E712
    register_part(
        outdir.parent / "manifest_parts.tsv", part="03", df=df, file=out.name,
        title="Curated pericyte program gene sets",
        supports="Figure 2; Figure 3; Methods",
        sources=("pericyte_states/_h/00.state_discovery.py (STATE_PANELS); "
                 "basement_membrane/_h/bm_panels.py; "
                 "pericyte_states/_m/pericyte_states.h5ad"),
        notes=("Panels are IMPORTED from the modules that define them, not copied, "
               "so this table cannot fork from the code that built the scores. "
               "`n_state_programs` counts membership across the six primary "
               "programs only; BM subpanels are subsets of basement_membrane and "
               "would inflate the count. "
               f"{n_absent} of {df['gene'].nunique()} curated genes are absent "
               "from the analysis object."))

    # Resolve the documented 8-vs-7 activated_migratory discrepancy explicitly.
    am = df[(df["panel_set"] == "state_program") &
            (df["program"] == "activated_migratory")]
    if am["detected_in_object"].notna().any():
        absent = am.loc[am["detected_in_object"] == False, "gene"].tolist()  # noqa: E712
        logging.info(f"activated_migratory: {len(am)} curated, "
                     f"{int((am['detected_in_object'] == True).sum())} present, "  # noqa: E712
                     f"absent = {absent or 'none'}")

    print("\nReproducibility information:")
    session_info.show()


if __name__ == "__main__":
    main()
