"""
Score basement-membrane (BM) deposition programs on the pericyte subset.

Scores the BM panel, its three structural sub-panels, and a fibrillar-ECM
contrast panel on pericyte_states/_m/pericyte_states.h5ad, using the same
sc.tl.score_genes idiom as pericyte_states/_h/00.state_discovery.py so the new
scores are directly comparable to the existing five program scores.

The canonical state model is NOT touched here: this step only adds scores.
Whether the BM panel would change `state_program` if it were allowed into the
relative-enrichment argmax is decided separately by 01.state_gate.py.

Per-gene expression and detection are emitted alongside the panel scores because
the selectivity and disease steps work gene-by-gene, and because a panel score
can hide a single dominant gene.

Outputs:
  - bm_metadata.tsv.gz   (per-cell scores + per-gene expr/detect; key `index`)
  - bm_panel_genes.tsv   (panel -> gene, the single source of truth for the R steps)
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from anndata import AnnData

import bm_panels


def configure_logging():
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def load_anndata(path: Path) -> AnnData:
    adata = sc.read_h5ad(path)
    if "logcounts" not in adata.layers:
        raise ValueError(f"{path} has no 'logcounts' layer")
    # var_names are already gene symbols in this object; guard anyway so a
    # future Ensembl-indexed input fails loudly instead of silently dropping
    # every gene as "absent".
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        if not (symbols.values == adata.var_names.values).all():
            adata.var_names = symbols
            adata.var_names_make_unique()
    return adata


def score_panels(adata, seed):
    """Mirrors pericyte_states/_h/00.state_discovery.py:score_panels."""
    adata.X = adata.layers["logcounts"]
    score_cols = []
    for panel, genes in bm_panels.PANELS.items():
        present = [g for g in genes if g in adata.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            logging.warning(f"[{panel}] dropping absent genes: {missing}")
        if not present:
            logging.warning(f"[{panel}] no genes present; skipping")
            continue
        col = f"{panel}_score"
        sc.tl.score_genes(adata, present, score_name=col, random_state=seed,
                          use_raw=False)
        score_cols.append(col)
        logging.info(f"[{panel}] scored {len(present)}/{len(genes)} genes")
    return score_cols


def per_gene_values(adata, genes):
    """Per-cell normalized expression and detection for each BM gene."""
    present = [g for g in genes if g in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if missing:
        logging.warning(f"per-gene: absent from object: {missing}")
    mat = adata[:, present].layers["logcounts"]
    mat = mat.toarray() if hasattr(mat, "toarray") else np.asarray(mat)
    out = {}
    for i, gene in enumerate(present):
        out[f"{gene}_expr"] = mat[:, i]
        out[f"{gene}_detect"] = (mat[:, i] > 0).astype(int)
    return pd.DataFrame(out, index=adata.obs_names), present


def main():
    configure_logging()
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    bm_panels.panel_table().to_csv(args.outdir / "bm_panel_genes.tsv",
                                   sep="\t", index=False)

    adata = load_anndata(args.adata)
    logging.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    score_cols = score_panels(adata, args.seed)
    gene_df, present = per_gene_values(adata, bm_panels.BM_PANEL)
    logging.info(f"Per-gene values for {len(present)}/{len(bm_panels.BM_PANEL)} "
                 "BM genes")

    # Carry the keys the R steps join on. Everything else already lives in
    # pericytes_states_metadata.tsv.gz and is merged there, not duplicated here.
    keep = [c for c in ["donor_id", "study", "dataset", "lung_condition",
                        "pericyte_state", "state_program", "total_counts",
                        "log10_total_counts"] if c in adata.obs.columns]
    meta = adata.obs[keep].copy()
    for col in score_cols:
        meta[col] = adata.obs[col].to_numpy()
    meta = pd.concat([meta, gene_df], axis=1)
    meta.index.name = "index"
    meta = meta.reset_index()

    out = args.outdir / "bm_metadata.tsv.gz"
    meta.to_csv(out, sep="\t", index=False)
    logging.info(f"Wrote {out} ({meta.shape[0]} rows x {meta.shape[1]} cols)")

    # Panel-level detection summary -- a fast read on which BM components
    # pericytes express at all, before any cross-cell-type comparison.
    det = pd.DataFrame({
        "gene": present,
        "detect_frac": [gene_df[f"{g}_detect"].mean() for g in present],
        "mean_expr": [gene_df[f"{g}_expr"].mean() for g in present],
        "mean_expr_in_pos": [
            gene_df.loc[gene_df[f"{g}_detect"] == 1, f"{g}_expr"].mean()
            for g in present],
    }).sort_values("detect_frac", ascending=False)
    det.to_csv(args.outdir / "bm_gene_detection_pericytes.tsv",
               sep="\t", index=False)
    logging.info("Pericyte BM detection:\n" + det.to_string(index=False))

    session_info.show()


if __name__ == "__main__":
    main()
