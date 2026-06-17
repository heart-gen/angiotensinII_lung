"""
Assign mouse pericyte/mural states on the integrated (scVI) representation and
quantify angiotensin receptors, for cross-species conservation.

Batch-robust design: state programs are scored on scVI batch-corrected
expression, and each cell's state is taken from its Leiden cluster (clustered on
the integrated latent) by the cluster's dominant program -- avoiding the
global per-cell z-score argmax that was confounded with dataset. A per-cell
argmax on corrected scores is also recorded for concordance.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse

STATE_PANELS_MOUSE = {
    "vascular_stabilizing": ["Rgs5", "Pdgfrb", "Notch3", "Kcnj8", "Abcc9",
                             "Higd1b", "Cox4i2", "Ndufa4l2", "Gja4", "Cspg4"],
    "inflammatory": ["Il6", "Ccl2", "Ccl20", "Cxcl1", "Cxcl2", "Cxcl3",
                     "Il1a", "Il1b", "Mif", "Icam1", "Vcam1", "Sele", "Nfkbia"],
    "synthetic_contractile": ["Acta2", "Myh11", "Tagln", "Cnn1", "Myl9", "Des", "Vim"],
    "activated_migratory": ["Adamts1", "Thbs1", "Timp1", "Mmp2", "Mmp3", "Serpine1", "Postn"],
    "fibroblast_like": ["Col1a1", "Col1a2", "Col3a1", "Col4a1", "Fn1", "Lum", "Dcn", "Pdgfa"],
}
MURAL_TYPES = ["pericyte", "vascular associated smooth muscle cell",
               "smooth muscle cell of the pulmonary artery"]


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="mouse_integrated.h5ad")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--cluster-key", default="leiden_scvi")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    # Score on batch-corrected expression
    if "scvi_corrected" in adata.layers:
        adata.X = adata.layers["scvi_corrected"]
    elif "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"]

    score_cols, mapping = [], []
    for state, genes in STATE_PANELS_MOUSE.items():
        present = [g for g in genes if g in adata.var_names]
        mapping.append({"state": state, "n_genes": len(genes), "n_present": len(present),
                        "present": ";".join(present)})
        if not present:
            continue
        col = f"{state}_score"
        sc.tl.score_genes(adata, present, score_name=col, random_state=args.seed, use_raw=False)
        score_cols.append(col)
    pd.DataFrame(mapping).to_csv(args.outdir / "ortholog_mapping_summary.tsv", sep="\t", index=False)
    states = [c.replace("_score", "") for c in score_cols]

    # Cluster-based assignment (batch-robust): cluster -> dominant program
    clu = adata.obs[args.cluster_key].astype(str)
    cluster_mean = adata.obs.groupby(clu, observed=True)[score_cols].mean()
    cluster_state = cluster_mean.apply(lambda r: states[int(np.argmax(r.values))], axis=1)
    adata.obs["pericyte_state"] = pd.Categorical(
        clu.map(cluster_state.to_dict()).values, categories=states)
    cluster_mean.assign(assigned=cluster_state).to_csv(
        args.outdir / "cluster_state_assignment.tsv", sep="\t")

    # Per-cell argmax on corrected scores (concordance reference)
    z = adata.obs[score_cols].to_numpy()
    z = (z - z.mean(0)) / (z.std(0) + 1e-9)
    adata.obs["pericyte_state_percell"] = pd.Categorical(
        [states[i] for i in z.argmax(1)], categories=states)
    concord = float((adata.obs["pericyte_state"] == adata.obs["pericyte_state_percell"]).mean())
    logging.info(f"cluster vs per-cell state concordance: {concord:.2f}")

    # Receptor expression (batch-corrected)
    for g in ["Agtr1a", "Agtr1b", "Agtr2"]:
        if g in adata.var_names:
            x = adata[:, g].X
            x = x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
            adata.obs[f"{g}_expr"] = x
            adata.obs[f"{g}_detect"] = (x > 0).astype(int)

    adata.obs["is_mural"] = adata.obs["cell_type"].astype(str).isin(MURAL_TYPES)
    logging.info(f"mural cells: {int(adata.obs['is_mural'].sum())}")

    keep = [c for c in ["cell_type", "dataset_id", "donor_id", "sex", "disease",
                        "is_mural", args.cluster_key, "pericyte_state",
                        "pericyte_state_percell"] + score_cols +
            [f"{g}_expr" for g in ["Agtr1a", "Agtr1b", "Agtr2"]] +
            [f"{g}_detect" for g in ["Agtr1a", "Agtr1b", "Agtr2"]]
            if c in adata.obs.columns]
    adata.obs[keep].to_csv(args.outdir / "mouse_states_metadata.tsv.gz", sep="\t")
    session_info.show()


if __name__ == "__main__":
    main()
