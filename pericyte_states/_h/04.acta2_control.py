"""
CONTROL / BENCHMARK (not a mechanistic pillar): does AGTR1 simply track canonical
ACTA2+ contractile mural identity?

This step only recomputes the one quantity that is NOT already in the per-cell
metadata: a LEAVE-ACTA2-OUT synthetic/contractile score. The full
synthetic_contractile panel is {ACTA2, MYH11, TAGLN, CNN1, MYL9, DES, VIM, PDGFRB}
(00.state_discovery.py); dropping ACTA2 leaves a broader contractile program so we
can ask, in 05.acta2_control.R, whether AGTR1 tracks contractile identity beyond the
single ACTA2 gene. Scored exactly like the curated panels (sc.tl.score_genes on the
logcounts layer, use_raw=False, seed=13) so it is comparable to the in-object scores.

Output is a minimal per-cell TSV (barcode + leave-out score) for a pure-TSV merge in
R; the h5ad read is routed through SLURM (step_4.sh) to keep the login node light.
"""
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import session_info

# synthetic/contractile panel from 00.state_discovery.py, minus the focal gene ACTA2
SYNTH_CONTRACTILE_NO_ACTA2 = ["MYH11", "TAGLN", "CNN1", "MYL9", "DES", "VIM", "PDGFRB"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--adata", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    adata.X = adata.layers["logcounts"]

    present = [g for g in SYNTH_CONTRACTILE_NO_ACTA2 if g in adata.var_names]
    missing = sorted(set(SYNTH_CONTRACTILE_NO_ACTA2) - set(present))
    if missing:
        logging.warning(f"leave-ACTA2-out contractile: dropping absent genes: {missing}")
    logging.info(f"leave-ACTA2-out contractile panel ({len(present)} genes): {present}")

    sc.tl.score_genes(adata, present, score_name="synth_contr_noACTA2_score",
                      random_state=args.seed, use_raw=False)

    out = pd.DataFrame(
        {"barcode": adata.obs_names,
         "synth_contr_noACTA2_score": np.asarray(adata.obs["synth_contr_noACTA2_score"])}
    )
    fn = args.outdir / "synth_contr_noACTA2.tsv.gz"
    out.to_csv(fn, sep="\t", index=False)
    logging.info(f"wrote {fn} ({len(out)} cells)")
    session_info.show()


if __name__ == "__main__":
    main()
