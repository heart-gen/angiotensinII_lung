"""
Marker dot-plot table for the state stability / annotation supplement (Figure S3).

00.state_discovery.py writes Wilcoxon markers with score / p / logFC, but a dot
plot also needs the per-cluster mean expression and detection rate, which it never
persists. This script recovers those WITHOUT re-running discovery: a full re-run
would rewrite pericyte_states.h5ad and the `pericyte_state` labels that
cell_communication, basement_membrane, niche_index and disease_association all key
on. Nothing downstream should move because a supplementary panel needed a number.

Scope note: the OTHER quantity that figure needs -- per-cluster bootstrap Jaccard
-- is produced by `00.state_discovery.py --stability-only` (driver: step_0c.sh),
which writes stability/cluster_bootstrap_jaccard.tsv. Do not duplicate it here.

The dot-plot table carries TWO kinds of gene block, distinguished by `block_type`:

  `marker`   top Wilcoxon genes per cluster -- what the data says separates P0-P5
  `curated`  genes from the STATE_PANELS used to annotate the clusters -- what the
             programs are defined by

Showing both is the honest version. The data-driven markers for the two largest
clusters are largely translational/housekeeping transcripts, so on their own they
would leave the reader unable to connect P0-P5 to the program labels; the curated
block supplies that link, and a gene appearing in both blocks is a validation
rather than a redundancy.

Output (pericyte_states/_m):
  - annotations/state_marker_dotplot.tsv      gene, block, block_type, block_rank,
                                              pericyte_state, mean_expr, frac_expr
"""
import argparse
import importlib.util
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from scipy import sparse

## Ribosomal-protein, mitochondrial and mitoribosomal genes dominate the Wilcoxon
## ranking of the two largest clusters (P0, P2) because those clusters differ
## partly in sequencing depth. They are excluded from the DISPLAYED markers only;
## nothing upstream is filtered, and state_markers.tsv.gz still holds the full
## unfiltered ranking.
HOUSEKEEPING_RE = re.compile(r"^(RP[LS]\d|RPLP|RPSA|MRP[LS]\d|MT-|MT[LR]?\d)")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="pericyte_states.h5ad written by 00.state_discovery.py")
    p.add_argument("--outdir", required=True, type=Path,
                   help="pericyte_states/_m (read annotations/, write back)")
    p.add_argument("--discovery", required=True, type=Path,
                   help="Path to 00.state_discovery.py, imported for STATE_PANELS")
    ## Equal block sizes: the figure lays the 12 blocks out in a facet grid with
    ## uniform column widths, so unequal counts would waste space and shrink the
    ## gene labels. Four per block keeps every label legible at print width.
    p.add_argument("--top-markers", type=int, default=4,
                   help="Displayed Wilcoxon markers per cluster, after filtering")
    p.add_argument("--curated-per-program", type=int, default=4,
                   help="Displayed curated panel genes per program")
    p.add_argument("--min-lfc", type=float, default=0.25)
    return p.parse_args()


def load_discovery(path: Path):
    """Import 00.state_discovery.py by path (the leading digits make it an
    invalid module name, so a plain import will not work)."""
    spec = importlib.util.spec_from_file_location("state_discovery", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pick_markers(markers: pd.DataFrame, top_n: int, min_lfc: float) -> pd.DataFrame:
    """Representative markers per cluster: significant, up, not housekeeping.

    A gene that tops more than one cluster is kept only for the cluster where it
    scores highest, so the dot plot has one block per cluster and no duplicated
    rows.
    """
    m = markers[(markers["pval_adj"] < 0.05) & (markers["logfoldchange"] >= min_lfc)]
    m = m[~m["gene"].astype(str).str.match(HOUSEKEEPING_RE)]
    m = m.sort_values("score", ascending=False)
    m = m.drop_duplicates(subset="gene", keep="first")
    ## m is already sorted by descending score, so head() takes the top n and
    ## cumcount() is the within-cluster rank.
    out = m.groupby("pericyte_state", sort=False).head(top_n).copy()
    out["block_rank"] = out.groupby("pericyte_state").cumcount() + 1
    out["block_type"] = "marker"
    return out.rename(columns={"pericyte_state": "block"})[
        ["gene", "block", "block_type", "block_rank"]]


def pick_curated(disc, adata, per_program: int) -> pd.DataFrame:
    """Representative genes from each curated program panel.

    Ranked by the SD of per-cluster mean expression, i.e. by how much the gene
    actually discriminates among P0-P5 -- most of the inflammatory panel is near
    zero in every pericyte cluster and would otherwise fill the block with blanks.
    Panels overlap (COL4A1 is in both fibroblast_like and basement_membrane;
    PDGFRB in both vascular_stabilizing and synthetic_contractile), so a gene is
    kept only for the first program that claims it, in STATE_PANELS order.
    """
    taken, rows = set(), []
    for program, genes in disc.STATE_PANELS.items():
        present = [g for g in genes if g in adata.var_names and g not in taken]
        if not present:
            logging.warning("[%s] no unclaimed panel genes present; skipping", program)
            continue
        st = expression_stats(adata, present)
        spread = st.groupby("gene")["mean_expr"].std().sort_values(ascending=False)
        for rank, gene in enumerate(spread.head(per_program).index, start=1):
            rows.append({"gene": gene, "block": program, "block_type": "curated",
                         "block_rank": rank})
            taken.add(gene)
    return pd.DataFrame(rows)


def expression_stats(adata, genes, layer="logcounts") -> pd.DataFrame:
    """Per-cluster mean expression and detection rate for the displayed genes."""
    present = [g for g in genes if g in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if missing:
        logging.warning("markers absent from var_names, dropped: %s", missing)
    sub = adata[:, present]
    mat = sub.layers[layer]
    mat = mat.toarray() if sparse.issparse(mat) else np.asarray(mat)
    expr = pd.DataFrame(mat, columns=present)
    expr["pericyte_state"] = adata.obs["pericyte_state"].astype(str).to_numpy()
    mean = expr.groupby("pericyte_state")[present].mean()
    frac = expr.groupby("pericyte_state")[present].agg(lambda s: (s > 0).mean())
    long = (mean.reset_index().melt("pericyte_state", var_name="gene",
                                    value_name="mean_expr")
            .merge(frac.reset_index().melt("pericyte_state", var_name="gene",
                                           value_name="frac_expr"),
                   on=["pericyte_state", "gene"]))
    return long


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    adir = args.outdir / "annotations"

    disc = load_discovery(args.discovery)
    adata = sc.read_h5ad(args.adata)
    logging.info("Loaded %d cells x %d genes", *adata.shape)

    # ---- dot-plot table: data-driven markers + curated panel genes --------
    markers = pd.read_csv(adir / "state_markers.tsv.gz", sep="\t")
    markers["pericyte_state"] = markers["pericyte_state"].astype(str)
    blocks = pd.concat([pick_markers(markers, args.top_markers, args.min_lfc),
                        pick_curated(disc, adata, args.curated_per_program)],
                       ignore_index=True)
    for kind, grp in blocks.groupby("block_type"):
        logging.info("Displayed %s genes:\n%s", kind,
                     grp.groupby("block")["gene"].apply(list).to_string())

    ## A gene can appear in both a marker block and a curated block (e.g. a
    ## cluster's top Wilcoxon hit that is also on its program's panel). Merging on
    ## `gene` alone would then fan out across blocks, so the stats are computed
    ## once on the unique gene set and joined back per block row.
    stats = expression_stats(adata, sorted(set(blocks["gene"])))
    dot = blocks.merge(stats, on="gene", how="inner")
    dot = dot.sort_values(["block_type", "block", "block_rank", "pericyte_state"])
    dot.to_csv(adir / "state_marker_dotplot.tsv", sep="\t", index=False)
    logging.info("wrote %s (%d rows, %d unique genes, %d blocks)",
                 adir / "state_marker_dotplot.tsv", len(dot),
                 dot["gene"].nunique(), dot["block"].nunique())

    session_info.show()


if __name__ == "__main__":
    main()
