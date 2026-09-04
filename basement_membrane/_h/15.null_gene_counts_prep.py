"""
Integer counts for the detection-matched null genes, from raw/X.

The AGTR1-versus-matrix REGRESSIONS all sit inside their detection-matched null
(14.agtr1_null_models.R). The count model is this module's designated arbiter
for those contrasts, so whether the biological claim survives depends on whether
the count model also sits inside its null -- and that has to be measured on the
same integer counts the arbiter uses, not on the soupX matrix.

This mirrors 09.agtr1_counts_prep.py exactly, but for the N matched genes chosen
by 13.agtr1_null_genes.py instead of for AGTR1 alone.

Output: null_gene_count_input.tsv.gz -- one column per null gene, plus the same
`raw_total_counts` library size the arbiter offsets on.
"""
import numpy as np
import pandas as pd
import h5py, logging, argparse
from pathlib import Path
from scipy import sparse
from anndata.io import read_elem


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--raw-adata", required=True, type=Path)
    p.add_argument("--null-genes", required=True, type=Path,
                   help="agtr1_null_genes.tsv from 13.agtr1_null_genes.py")
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--out", type=Path, default=Path("./null_gene_count_input.tsv.gz"))
    return p.parse_args()


def symbol_positions(var_raw, var_main, wanted):
    """Map gene symbols to raw/X column indices.

    raw/var is keyed on Ensembl IDs while the symbols live in the main var, so
    the mapping goes symbol -> main row -> positional index. 09's locate_gene
    asserts the two vars are the same length and aligned; the same assertion is
    the only thing that makes the positional lookup valid, so it is repeated
    here rather than assumed.
    """
    if len(var_main) != len(var_raw):
        raise ValueError("raw and main var have different lengths; the "
                         "positional mapping below would be silently wrong")
    sym = None
    for cand in ("feature_name", "gene_symbols", "symbol"):
        if cand in var_main.columns:
            sym = var_main[cand].astype(str).values
            break
    if sym is None:
        sym = var_main.index.astype(str).values
    lookup = {}
    for i, s in enumerate(sym):
        lookup.setdefault(s, i)
    out, missing = {}, []
    for g in wanted:
        if g in lookup:
            out[g] = lookup[g]
        else:
            missing.append(g)
    if missing:
        logging.warning("%d null genes absent from raw var, dropped: %s",
                        len(missing), missing[:5])
    return out


def main():
    configure_logging()
    args = parse_args()

    meta = pd.read_csv(args.metadata, sep="\t").set_index("index")
    genes = pd.read_csv(args.null_genes, sep="\t")["gene"].astype(str).tolist()
    logging.info("null genes requested: %d", len(genes))

    with h5py.File(args.raw_adata, "r") as fh:
        obs = read_elem(fh["obs"])
        var_main = read_elem(fh["var"])
        var_raw = read_elem(fh["raw"]["var"])
        pos = symbol_positions(var_raw, var_main, genes)
        X = sparse.csr_matrix(read_elem(fh["raw"]["X"]))

    if not np.allclose(X.data, np.round(X.data)):
        raise ValueError("raw/X is not integer-valued; the count model "
                         "requires true counts")

    cols = {"raw_total_counts": np.asarray(X.sum(axis=1)).ravel().astype(np.int64)}
    idx = sorted(pos.values())
    sub = np.asarray(X[:, idx].todense())
    for g, j in zip([g for g in genes if g in pos],
                    [idx.index(pos[g]) for g in genes if g in pos]):
        cols[f"null_{g}"] = sub[:, j].astype(np.int64)
    out = pd.DataFrame(cols, index=obs.index)

    missing = out.index.difference(meta.index)
    if len(missing):
        logging.warning("%d raw cells absent from metadata", len(missing))
    out = out.loc[out.index.intersection(meta.index)]
    logging.info("writing %d cells x %d null genes", out.shape[0], out.shape[1] - 1)
    out.reset_index(names="index").to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
