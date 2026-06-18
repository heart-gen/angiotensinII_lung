"""
Export the pericyte UMAP coordinates (and dpt pseudotime if present) from
pericyte_states.h5ad as a small per-cell TSV, so the integrated pericyte-layer
figure (pericyte_layer_figure.R) can overlay localization "where" markers and
the state/continuum "what/why" readouts on the SAME embedding in ggplot/theme_ms.

The embedding is the one localization/pericyte_analysis built (obsm["X_umap"]),
which pericyte_states reads without recomputing -- so localization gene-UMAPs and
the state/pseudotime UMAPs are identical coordinates. All other per-cell fields
(program scores, state_program, AGTR1/ACTA2, etc.) already live in
pericyte_states/_m/pericytes_states_metadata.tsv.gz and are joined in R by barcode.

Output (figures/_m, h5ad read routed through SLURM, step_figures.sh):
  - pericyte_umap_coords.tsv.gz   barcode, UMAP1, UMAP2 (+ dpt_pseudotime if in obs)
"""
import argparse
import logging
from pathlib import Path

import pandas as pd
import scanpy as sc
import session_info


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
    if "X_umap" not in adata.obsm:
        raise KeyError(f"X_umap not in obsm: {list(adata.obsm)}")
    um = adata.obsm["X_umap"]
    out = pd.DataFrame({"barcode": adata.obs_names,
                        "UMAP1": um[:, 0], "UMAP2": um[:, 1]})
    for pt in ("dpt_pseudotime", "pseudotime"):
        if pt in adata.obs:
            out["dpt_pseudotime"] = adata.obs[pt].to_numpy()
            logging.info("included pseudotime column %r from obs", pt)
            break
    fout = args.outdir / "pericyte_umap_coords.tsv.gz"
    out.to_csv(fout, sep="\t", index=False)
    logging.info("wrote %s (%d cells, cols=%s)", fout, len(out), list(out.columns))
    session_info.show()


if __name__ == "__main__":
    main()
