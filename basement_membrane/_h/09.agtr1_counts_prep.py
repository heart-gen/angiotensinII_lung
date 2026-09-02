"""Extract integer AGTR1 counts per pericyte for the count-model lens.

WHY THIS EXISTS. The project's AGTR1 readouts are all derived quantities:
`AGTR1_expr` and `AGTR1_detect` come from the soupX ambient-corrected matrix,
and `AGTR1_scvi` is a model output. The count-model lens asks the association
question without imputing anything -- AGTR1 counts are the RESPONSE of a
negative-binomial GLMM with a library-size offset, so dropout is handled inside
the likelihood rather than removed beforehand. That requires true integer
counts, which survive only in `raw/X` of pericyte.hlca_full.dataset.h5ad (the
`counts`/`soupX` layers are float, minimum nonzero 1.86e-4, and only ~2% of
their stored values are integers).

Emits one row per pericyte: the integer AGTR1 count, the integer library size
that becomes the offset, and the donor/study/cluster keys plus the matrix scores
that 10.agtr1_count_models.R models against.

Alignment is asserted, never assumed: raw var is keyed on Ensembl IDs while the
scored object is keyed on symbols, so the AGTR1 column is located through
var['ensembl_id'] and the resulting symbol is checked.
"""
import logging
import argparse
import h5py
import numpy as np
import pandas as pd
import session_info
from pathlib import Path
from scipy import sparse
from anndata.io import read_elem

GENE = "AGTR1"


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--raw-adata", type=Path,
                        default=Path("../../localization/pericyte_analysis/_m/"
                                     "pericyte.hlca_full.dataset.h5ad"),
                        help="Pericyte AnnData whose raw/X holds integer counts")
    parser.add_argument("--metadata", type=Path,
                        default=Path("./bm_metadata.tsv.gz"),
                        help="Per-cell BM metadata (scores, cluster, donor, study)")
    parser.add_argument("--out", type=Path, default=Path("./agtr1_count_input.tsv.gz"))
    return parser.parse_args()


def locate_gene(var_raw, var_main, gene=GENE):
    """Column index of `gene` in raw/X.

    The two HLCA-derived objects in this project are keyed differently -- one on
    symbols with an `ensembl_id` column, the other on Ensembl IDs with a
    `feature_name` column -- so resolve through whichever is present and then
    cross-check the answer against the other table.
    """
    def _find(frame, col):
        if col == "<index>":
            vals = frame.index.astype(str).values
        elif col in frame:
            vals = frame[col].astype(str).values
        else:
            return None
        hit = np.flatnonzero(vals == gene)
        return int(hit[0]) if hit.size else None

    pos = _find(var_raw, "feature_name")
    if pos is None:
        pos = _find(var_raw, "<index>")
    if pos is None and "ensembl_id" in var_main:
        ens_hit = _find(var_main, "<index>")
        if ens_hit is not None:
            ens = str(var_main["ensembl_id"].iloc[ens_hit])
            pos = int(np.flatnonzero(var_raw.index.astype(str).values == ens)[0])
    if pos is None:
        raise KeyError(f"{gene} could not be located in raw/var")

    ## Cross-check against the main var table, which must agree gene-for-gene
    if len(var_main) != len(var_raw):
        raise ValueError("raw and main var have different lengths")
    for col in ("feature_name", "original_gene_symbols"):
        if col in var_main:
            got = str(var_main[col].iloc[pos])
            if got != gene:
                raise ValueError(f"main var {col} at column {pos} is {got}, not {gene}")
            break
    logging.info("%s resolved to raw column %d (%s)", gene, pos, var_raw.index[pos])
    return pos


def main():
    args = parse_args()
    configure_logging()

    meta = pd.read_csv(args.metadata, sep="\t")
    meta = meta.set_index("index")
    logging.info("metadata: %d cells", len(meta))

    with h5py.File(args.raw_adata, "r") as fh:
        obs = read_elem(fh["obs"])
        var_main = read_elem(fh["var"])
        var_raw = read_elem(fh["raw"]["var"])
        pos = locate_gene(var_raw, var_main)
        X = sparse.csr_matrix(read_elem(fh["raw"]["X"]))

    if not np.allclose(X.data, np.round(X.data)):
        raise ValueError("raw/X is not integer-valued; the count model requires "
                         "true counts")

    counts = np.asarray(X[:, pos].todense()).ravel()
    libsize = np.asarray(X.sum(axis=1)).ravel()
    out = pd.DataFrame({f"{GENE}_count": counts.astype(np.int64),
                        "raw_total_counts": libsize.astype(np.int64)},
                       index=obs.index)

    ## Cells must line up with the metadata; report rather than silently inner-join
    missing = out.index.difference(meta.index)
    if len(missing):
        logging.warning("%d cells in raw object absent from metadata", len(missing))
    df = meta.join(out, how="inner")
    if len(df) != len(meta):
        raise ValueError(f"joined {len(df)} of {len(meta)} metadata cells")

    logging.info("%s detected in %d/%d pericytes (%.1f%%); median library %d",
                 GENE, int((df[f"{GENE}_count"] > 0).sum()), len(df),
                 100 * (df[f"{GENE}_count"] > 0).mean(), int(df["raw_total_counts"].median()))
    logging.info("count distribution: %s",
                 df[f"{GENE}_count"].value_counts().head(6).to_dict())

    keep = ["donor_id", "study", "dataset", "lung_condition", "pericyte_state",
            "state_program", "log10_total_counts", f"{GENE}_count", "raw_total_counts",
            "basement_membrane_score", "fibrillar_collagen_score", "fibrillar_ecm_score",
            "tgfb_response_score", "ambient_tracer_score"]
    keep = [c for c in keep if c in df]
    df.reset_index()[["index"] + keep].to_csv(args.out, sep="\t", index=False)
    logging.info("wrote %s (%d rows)", args.out, len(df))
    session_info.show()


if __name__ == "__main__":
    main()
