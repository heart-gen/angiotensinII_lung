"""
Build a donor-aware cell-cell communication (CCC) niche from the full HLCA
disease object.

The niche keeps the vascular / alveolar / stromal / recruited-immune cells that
plausibly signal to pericytes and alveolar epithelium. Three receiver grouping
schemes are emitted so the main analysis does NOT rest on binarizing a dropout-
prone GPCR transcript (which would contradict the manuscript's own AGTR1 dropout
result):

  ccc_group          (MAIN): pericytes whole ("Pericytes"); AT2 split by AGTR2
                     detectability ("AT2_AGTR2det" / "AT2_AGTR2undet").
  ccc_group_state    (STRATIFIED): pericytes split by functional state from
                     pericyte_states ("Pericyte_<program>"); AT2 as above.
  ccc_group_receptor (SENSITIVITY only): pericytes by AGTR1 detectability
                     ("Pericyte_AGTR1det" / "Pericyte_AGTR1undet"); AT2 as above.

Continuous AGTR1/AGTR2 expression is also carried (for continuous-score analyses
downstream) rather than used as a receiver identity in the main analysis.

To bound liana runtime while preserving donor structure, cells are capped per
(cell_type, donor). Disease is harmonized to {Healthy, COPD, Fibrotic_ILD, Other}
from lung_condition.

Outputs:
  - ccc_niche.h5ad           (logcounts in X; ccc_group{,_state,_receptor} + disease + donor)
  - ccc_niche_obs.tsv.gz     (cell metadata for auditing)
  - ccc_group_counts.tsv     (main group x disease), ccc_group_state_counts.tsv
"""
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import session_info
import logging, argparse
from pathlib import Path
from scipy import sparse

# cell_type labels (HLCA full) to retain as the signaling niche
NICHE_CELL_TYPES = [
    # mural / stroma
    "Pericytes", "Vascular smooth muscle",
    "Alveolar fibroblasts", "Adventitial fibroblasts", "Peribronchial fibroblasts",
    "Myofibroblasts", "Subpleural fibroblasts",
    # endothelium (vascular niche partners)
    "EC general capillary", "EC aerocyte capillary", "EC arterial",
    "EC venous systemic", "EC venous pulmonary", "Lymphatic EC",
    # alveolar epithelium
    "AT1", "AT2", "Transitional Club-AT2",
    # recruited / resident immune ligand sources
    "Alveolar macrophages", "Interstitial macrophages",
    "Classical monocytes", "Non-classical monocytes", "Mast cells", "DC2",
]


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path,
                   help="Full HLCA disease object (hlca_full.dataset.h5ad)")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--pericyte-states", type=Path, default=None,
                   help="pericytes_states_metadata.tsv.gz (barcode-indexed "
                        "state_program) for the state-stratified receiver scheme")
    p.add_argument("--celltype-key", default="cell_type")
    p.add_argument("--cap-per-group-donor", type=int, default=300,
                   help="Max cells per (cell_type, donor) to bound runtime")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def map_disease_group(lc: pd.Series) -> pd.Series:
    lc = lc.astype(str)
    out = np.where(lc.str.startswith("Healthy"), "Healthy",
          np.where(lc.eq("COPD"), "COPD",
          np.where(lc.str.contains("IPF|fibrosis|ILD|NSIP|Sarcoid|^HP$|Lymphangio|sclerosis",
                                    case=False, regex=True), "Fibrotic_ILD", "Other")))
    return pd.Series(out, index=lc.index)


def gene_value(adata, gene: str) -> np.ndarray:
    """Continuous log-normalized expression for a gene (0 if absent)."""
    if gene not in adata.var_names:
        logging.warning(f"{gene} absent; treating as 0")
        return np.zeros(adata.n_obs, dtype=float)
    sub = adata[:, gene]
    x = sub.layers["logcounts"] if "logcounts" in sub.layers else sub.X
    x = x.toarray().ravel() if sparse.issparse(x) else np.asarray(x).ravel()
    return np.asarray(x, dtype=float)


def load_pericyte_programs(path, obs_names) -> pd.Series:
    """Map cell barcode -> pericyte functional state (state_program) from the
    pericyte_states metadata. Returns a Series aligned to obs_names (NaN where
    unmatched). Logs the match fraction so a barcode mismatch is visible."""
    if path is None or not Path(path).exists():
        logging.warning("No pericyte-states metadata; state scheme will be unassigned")
        return pd.Series(index=obs_names, dtype=object)
    meta = pd.read_csv(path, sep="\t", index_col=0)
    col = "state_program" if "state_program" in meta.columns else None
    if col is None:
        logging.warning("state_program not in pericyte-states metadata; unassigned")
        return pd.Series(index=obs_names, dtype=object)
    prog = meta[col].astype(str)
    mapped = pd.Series(obs_names, index=obs_names).map(prog)
    logging.info(f"pericyte state_program matched {mapped.notna().sum()}/"
                 f"{len(meta)} state cells onto niche barcodes")
    return mapped


def main():
    args = parse_args()
    configure_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logging.info("Opening backed AnnData (obs-only subset first)...")
    backed = ad.read_h5ad(args.adata, backed="r")
    ct = backed.obs[args.celltype_key].astype(str)
    keep_mask = ct.isin(NICHE_CELL_TYPES).to_numpy()
    logging.info(f"Niche cells: {keep_mask.sum()} / {backed.n_obs}")

    # Cap per (cell_type, donor) on the backed obs before loading matrix
    obs = backed.obs.loc[keep_mask, [args.celltype_key, "donor_id"]].copy()
    rng = np.random.default_rng(args.seed)
    take_idx = (
        obs.assign(_i=np.arange(obs.shape[0]))
           .groupby([args.celltype_key, "donor_id"], observed=True)["_i"]
           .apply(lambda s: s.sample(min(len(s), args.cap_per_group_donor),
                                     random_state=args.seed))
           .explode().dropna().astype(int).to_numpy()
    )
    take_idx.sort()
    global_idx = np.where(keep_mask)[0][take_idx]
    logging.info(f"After capping {args.cap_per_group_donor}/group/donor: {len(global_idx)} cells")

    adata = backed[global_idx].to_memory()
    del backed

    # var_names -> symbols; logcounts into X for liana
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        adata.var["ensembl_id"] = adata.var_names
        adata.var_names = adata.var.index = symbols
    adata.raw = None
    adata.var_names_make_unique()
    # Full HLCA object stores log-normalized values in X (no 'logcounts' layer);
    # alias logcounts to X so downstream liana/nichenet find lognorm data.
    if "logcounts" not in adata.layers:
        adata.layers["logcounts"] = adata.X
    adata.X = adata.layers["logcounts"]

    # Receptor expression (continuous) + detectability (binary, sparse GPCR)
    base = adata.obs[args.celltype_key].astype(str)
    agtr1_v = gene_value(adata, "AGTR1"); agtr2_v = gene_value(adata, "AGTR2")
    agtr1 = (agtr1_v > 0).astype(int); agtr2 = (agtr2_v > 0).astype(int)
    is_peri = base.eq("Pericytes").to_numpy()
    is_at2 = base.eq("AT2").to_numpy()

    # AT2 split by AGTR2 detectability (cautious framing); used in all schemes.
    at2_label = np.where(agtr2 > 0, "AT2_AGTR2det", "AT2_AGTR2undet")

    # MAIN: pericytes whole; AT2 by AGTR2 detectability.
    g_main = base.to_numpy().astype(object)
    g_main[is_peri] = "Pericytes"
    g_main[is_at2] = at2_label[is_at2]
    adata.obs["ccc_group"] = pd.Categorical(g_main)

    # STATE: pericytes by functional state (state_program); AT2 as above.
    prog = load_pericyte_programs(args.pericyte_states, adata.obs_names)
    adata.obs["pericyte_program"] = prog.to_numpy()
    peri_state = np.where(prog.isna().to_numpy(), "Pericyte_unassigned",
                          "Pericyte_" + prog.astype(str).to_numpy())
    g_state = base.to_numpy().astype(object)
    g_state[is_peri] = peri_state[is_peri]
    g_state[is_at2] = at2_label[is_at2]
    adata.obs["ccc_group_state"] = pd.Categorical(g_state)

    # RECEPTOR (sensitivity): pericytes by AGTR1 detectability; AT2 as above.
    g_recep = base.to_numpy().astype(object)
    g_recep[is_peri] = np.where(agtr1[is_peri] > 0, "Pericyte_AGTR1det",
                                "Pericyte_AGTR1undet")
    g_recep[is_at2] = at2_label[is_at2]
    adata.obs["ccc_group_receptor"] = pd.Categorical(g_recep)

    adata.obs["disease_group"] = pd.Categorical(map_disease_group(adata.obs["lung_condition"]))
    adata.obs["AGTR1_detect"] = agtr1
    adata.obs["AGTR2_detect"] = agtr2
    adata.obs["AGTR1_expr"] = agtr1_v
    adata.obs["AGTR2_expr"] = agtr2_v

    # Summaries
    pd.crosstab(adata.obs["ccc_group"], adata.obs["disease_group"]).to_csv(
        args.outdir / "ccc_group_counts.tsv", sep="\t")
    pd.crosstab(adata.obs["ccc_group_state"], adata.obs["disease_group"]).to_csv(
        args.outdir / "ccc_group_state_counts.tsv", sep="\t")
    pd.crosstab(adata.obs["ccc_group_receptor"], adata.obs["disease_group"]).to_csv(
        args.outdir / "ccc_group_receptor_counts.tsv", sep="\t")
    audit_cols = [c for c in [args.celltype_key, "ccc_group", "ccc_group_state",
                              "ccc_group_receptor", "pericyte_program",
                              "disease_group", "donor_id", "disease",
                              "lung_condition", "smoking_status",
                              "AGTR1_detect", "AGTR2_detect",
                              "AGTR1_expr", "AGTR2_expr"]
                  if c in adata.obs.columns]
    adata.obs[audit_cols].to_csv(args.outdir / "ccc_niche_obs.tsv.gz", sep="\t")

    if adata.var.index.name in adata.var.columns:
        col = adata.var.index.name
        if not np.array_equal(adata.var.index.to_numpy(), adata.var[col].to_numpy()):
            adata.var.index.name = None
    adata.write(args.outdir / "ccc_niche.h5ad")
    logging.info(f"Wrote niche: {adata.shape}")

    session_info.show()


if __name__ == "__main__":
    main()
