"""
Program x protein-category enrichment for the pericyte marker panel.

Replaces the (uninformative) state->category->role alluvial flow with a
quantitative program x protein-category readout. For each pericyte the six
continuous program scores are z-scored across cells and the dominant program is
taken by argmax (the same relative-enrichment logic used to annotate the stable
clusters in pericyte_states/00.state_discovery.py); this keeps all SIX programs
as groups (the discrete clustering only yields three). We then compute, per
program group, the mean fraction of cells expressing each marker, aggregated to
the eight protein categories. A second table gives per-gene detection by program
(for the curated gene -> category -> role taxonomy alluvial).

Outputs (small TSVs; h5ad read routed through SLURM, step_3b.sh):
  - program_category_enrichment.tsv   program x category mean detection (+ n_cells)
  - gene_program_detection.tsv        gene x program mean detection + category/role
                                      + an `overall` prevalence row per gene
"""
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from scipy import sparse

# gene -> (category, candidate_role); identical to 03a.marker_fractions.py
PANEL = {
    "TGFB1": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "TGFB2": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "TGFB3": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "WNT5A": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "NODAL": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "CXCL1": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL2": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL3": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL6": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL8": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL10": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CCL2": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CCL20": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "IL1A": ("Cytokines", "Inflammatory amplification"),
    "IL1B": ("Cytokines", "Inflammatory amplification"),
    "IL6": ("Cytokines", "Inflammatory amplification"),
    "MIF": ("Cytokines", "Inflammatory amplification"),
    "ICAM1": ("Adhesion molecules", "Leukocyte trafficking"),
    "VCAM1": ("Adhesion molecules", "Leukocyte trafficking"),
    "SELE": ("Adhesion molecules", "Leukocyte trafficking"),
    "COL1A1": ("ECM structural", "Basement membrane / vascular stability"),
    "COL4A1": ("ECM structural", "Basement membrane / vascular stability"),
    "FN1": ("ECM structural", "Basement membrane / vascular stability"),
    "LAMA2": ("ECM structural", "Basement membrane / vascular stability"),
    "LAMB1": ("ECM structural", "Basement membrane / vascular stability"),
    "MMP2": ("Matrix-remodeling", "ECM degradation / sprouting"),
    "MMP3": ("Matrix-remodeling", "ECM degradation / sprouting"),
    "MMP9": ("Matrix-remodeling", "ECM degradation / sprouting"),
    "ACTA2": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "PDGFA": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "PDGFRB": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "ADAMTS1": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "VIM": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "RGS5": ("Mural identity", "Vascular support"),
}

PROGRAM_SCORES = {
    "vascular_stabilizing": "vascular_stabilizing_score",
    "inflammatory": "inflammatory_score",
    "synthetic_contractile": "synthetic_contractile_score",
    "activated_migratory": "activated_migratory_score",
    "fibroblast_like": "fibroblast_like_score",
    "basement_membrane": "basement_membrane_score",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    X = adata.layers.get("logcounts", adata.X)

    # dominant program by z-scored argmax (relative enrichment, not raw magnitude)
    score_cols = list(PROGRAM_SCORES.values())
    missing = [c for c in score_cols if c not in adata.obs]
    if missing:
        raise KeyError(f"missing program score columns: {missing}")
    Z = adata.obs[score_cols].to_numpy(dtype=float)
    Z = (Z - Z.mean(0)) / Z.std(0)
    prog_names = np.array(list(PROGRAM_SCORES.keys()))
    program = prog_names[np.argmax(Z, axis=1)]
    adata.obs["dominant_program"] = program
    logging.info("dominant-program cell counts:\n%s",
                 pd.Series(program).value_counts().to_string())

    genes = [g for g in PANEL if g in adata.var_names]
    logging.info(f"markers present: {len(genes)}/{len(PANEL)}")
    Xg = X[:, [adata.var_names.get_loc(g) for g in genes]]
    detect = (Xg > 0)
    detect = np.asarray(detect.todense()) if sparse.issparse(detect) else np.asarray(detect)
    det = pd.DataFrame(detect.astype(float), columns=genes, index=adata.obs_names)
    det["program"] = program

    cat = {g: PANEL[g][0] for g in genes}
    role = {g: PANEL[g][1] for g in genes}

    # (1) gene x program mean detection (+ overall prevalence)
    gp = det.groupby("program")[genes].mean().T            # gene x program
    gp["overall"] = det[genes].mean()
    gp = gp.reset_index().rename(columns={"index": "gene"}).melt(
        id_vars="gene", var_name="program", value_name="mean_detect")
    gp["category"] = gp["gene"].map(cat)
    gp["role"] = gp["gene"].map(role)
    gp.to_csv(args.outdir / "gene_program_detection.tsv", sep="\t", index=False)

    # (2) program x category mean detection: per cell average detection over the
    # genes in each category, then mean across the program's cells.
    cats = sorted(set(cat.values()))
    catmean = pd.DataFrame(index=det.index)
    for c in cats:
        cg = [g for g in genes if cat[g] == c]
        catmean[c] = det[cg].mean(axis=1)
    catmean["program"] = program
    pc = catmean.groupby("program")[cats].mean()
    n = catmean.groupby("program").size().rename("n_cells")
    pc = pc.reset_index().melt(id_vars="program", var_name="category",
                               value_name="mean_detect").merge(n, on="program")
    pc.to_csv(args.outdir / "program_category_enrichment.tsv", sep="\t", index=False)
    logging.info("wrote program_category_enrichment.tsv (%d rows) and "
                 "gene_program_detection.tsv (%d rows)", len(pc), len(gp))
    session_info.show()


if __name__ == "__main__":
    main()
