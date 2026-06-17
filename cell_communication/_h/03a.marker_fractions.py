"""
Export per-cluster / per-state marker detection fractions for the alluvial plot.

For each pericyte Leiden subcluster and each functional state, compute the
fraction of cells expressing (>0) each marker in the signaling/role panels.
Feeds 03.figures.R (ggalluvial: cluster/state -> protein category -> role).
"""
import numpy as np
import pandas as pd
import scanpy as sc
import logging, argparse
from pathlib import Path
from scipy import sparse

# Marker -> (category, candidate_role) from the study's signaling tables.
PANEL = {
    # Signaling ligands -> pericyte-endothelial crosstalk
    "TGFB1": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "TGFB2": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "TGFB3": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "WNT5A": ("Signaling ligands", "Pericyte-EC crosstalk"),
    "NODAL": ("Signaling ligands", "Pericyte-EC crosstalk"),
    # Chemokines -> neutrophil/monocyte recruitment
    "CXCL1": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL2": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL3": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL6": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL8": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CXCL10": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CCL2": ("Chemokines", "Neutrophil/monocyte recruitment"),
    "CCL20": ("Chemokines", "Neutrophil/monocyte recruitment"),
    # Cytokines -> inflammatory amplification
    "IL1A": ("Cytokines", "Inflammatory amplification"),
    "IL1B": ("Cytokines", "Inflammatory amplification"),
    "IL6": ("Cytokines", "Inflammatory amplification"),
    "MIF": ("Cytokines", "Inflammatory amplification"),
    # Adhesion -> leukocyte trafficking
    "ICAM1": ("Adhesion molecules", "Leukocyte trafficking"),
    "VCAM1": ("Adhesion molecules", "Leukocyte trafficking"),
    "SELE": ("Adhesion molecules", "Leukocyte trafficking"),
    # ECM structural -> basement membrane integrity
    "COL1A1": ("ECM structural", "Basement membrane / vascular stability"),
    "COL4A1": ("ECM structural", "Basement membrane / vascular stability"),
    "FN1": ("ECM structural", "Basement membrane / vascular stability"),
    "LAMA2": ("ECM structural", "Basement membrane / vascular stability"),
    "LAMB1": ("ECM structural", "Basement membrane / vascular stability"),
    # Matrix-remodeling enzymes -> ECM degradation / sprouting
    "MMP2": ("Matrix-remodeling", "ECM degradation / sprouting"),
    "MMP3": ("Matrix-remodeling", "ECM degradation / sprouting"),
    "MMP9": ("Matrix-remodeling", "ECM degradation / sprouting"),
    # Fibrotic mediators -> myofibroblast transition / fibrosis
    "ACTA2": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "PDGFA": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "PDGFRB": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "ADAMTS1": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
    "RGS5": ("Mural identity", "Vascular support"),
    "VIM": ("Fibrotic mediators", "Myofibroblast transition / fibrosis"),
}


def main():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    X = adata.layers.get("logcounts", adata.X)
    genes = [g for g in PANEL if g in adata.var_names]
    logging.info(f"markers present: {len(genes)}/{len(PANEL)}")

    rows = []
    for key in ["leiden_pericytes", "pericyte_state"]:
        if key not in adata.obs:
            continue
        for grp in adata.obs[key].astype(str).unique():
            m = (adata.obs[key].astype(str) == grp).to_numpy()
            sub = X[m][:, [adata.var_names.get_loc(g) for g in genes]]
            frac = (np.asarray((sub > 0).mean(axis=0)).ravel()
                    if sparse.issparse(sub) else (sub > 0).mean(axis=0))
            for g, f in zip(genes, frac):
                cat, role = PANEL[g]
                rows.append({"grouping": key, "cluster": grp, "gene": g,
                             "category": cat, "role": role, "frac_expr": float(f)})
    out = pd.DataFrame(rows)
    out.to_csv(args.outdir / "marker_fractions_by_cluster.tsv.gz", sep="\t", index=False)
    logging.info(f"wrote {out.shape[0]} rows")


if __name__ == "__main__":
    main()
