"""
Disease-stratified cell-cell communication with liana (consensus rank_aggregate).

Runs liana separately within each disease_group so edges are comparable across
Healthy / COPD / Fibrotic_ILD, for each receiver grouping scheme built in
00.prepare_ccc_input.py:
    main     -- pericytes whole ("Pericytes"); AT2 by AGTR2 detectability
    state    -- pericytes by functional state ("Pericyte_<program>")
    receptor -- pericytes by AGTR1 detectability (SENSITIVITY only)
The main + state analyses are the biologically clean primaries; receptor is a
sensitivity check that does not anchor the story on binarizing a dropout-prone
GPCR transcript.

For each scheme it focuses on ligand->receptor edges INTO the pericyte / AT2
receivers and scores a curated injury-ligand panel (seeds 02.nichenet.R).

Outputs (per scheme):
  - liana_res_<scheme>_<group>.tsv.gz   full ranked LR table per disease group
  - liana_into_receivers_<scheme>.tsv.gz
  - ligand_panel_into_receivers_<scheme>.tsv.gz
  - expressed_fraction_<scheme>.tsv.gz  (main also -> expressed_fraction_per_group.tsv.gz)
  - figures/dotplot_into_<receiver>.{pdf,png}
"""
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import session_info
import logging, argparse
from pathlib import Path

# scheme -> obs grouping column
SCHEME_KEYS = {"main": "ccc_group", "state": "ccc_group_state",
               "receptor": "ccc_group_receptor",
               # orthogonal data-driven scheme: pericytes by dominant CoGAPS
               # pattern (00b.cogaps_receivers.py); robustness check on `state`.
               "cogaps": "ccc_group_cogaps"}

LIGAND_PANEL = [
    "TGFB1", "TGFB2", "TGFB3", "WNT5A", "NODAL",
    "CXCL1", "CXCL2", "CXCL3", "CXCL6", "CXCL8", "CXCL10", "CCL2", "CCL20",
    "IL1A", "IL1B", "IL6", "MIF", "ICAM1", "VCAM1", "SELE",
]


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("--adata", required=True, type=Path, help="ccc_niche.h5ad")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--schemes", default="main,state,receptor",
                   help="Comma-separated receiver schemes to run")
    p.add_argument("--disease-key", default="disease_group")
    p.add_argument("--expr-prop", type=float, default=0.1)
    p.add_argument("--min-cells", type=int, default=10,
                   help="Min cells per group to keep it in a disease subset")
    p.add_argument("--seed", type=int, default=13)
    return p.parse_args()


def pick_receivers(categories):
    """Pericyte / AT2 receiver groups present in this scheme's grouping."""
    recv = []
    for c in categories:
        if c == "Pericytes" or (c.startswith("Pericyte_") and not c.endswith("unassigned")):
            recv.append(c)
        elif c.startswith("AT2_AGTR2"):
            recv.append(c)
    return recv


def export_expressed_fractions(adata, group_key, path, min_frac=0.01):
    from scipy import sparse
    X = adata.layers.get("logcounts", adata.X)
    detect = (X > 0)
    groups = adata.obs[group_key].astype(str)
    rows = {}
    for g in groups.unique():
        m = (groups == g).to_numpy()
        d = detect[m]
        frac = np.asarray(d.mean(axis=0)).ravel() if sparse.issparse(d) else d.mean(axis=0)
        rows[g] = frac
    frac_df = pd.DataFrame(rows, index=adata.var_names)
    frac_df = frac_df[frac_df.max(axis=1) >= min_frac]
    frac_df.to_csv(path, sep="\t")
    logging.info(f"Saved expressed fractions -> {path.name}: {frac_df.shape}")


def run_one(adata, group_key, expr_prop, min_cells, seed):
    vc = adata.obs[group_key].value_counts()
    keep = vc[vc >= min_cells].index
    sub = adata[adata.obs[group_key].isin(keep)].copy()
    sub.obs[group_key] = sub.obs[group_key].astype("category").cat.remove_unused_categories()
    if sub.obs[group_key].nunique() < 2:
        logging.warning("Fewer than 2 groups; skipping")
        return None
    li.mt.rank_aggregate(
        sub, groupby=group_key, expr_prop=expr_prop, use_raw=False,
        return_all_lrs=False, verbose=False, seed=seed,
    )
    return sub.uns["liana_res"]


def run_scheme(adata, scheme, group_key, args, fig_dir):
    logging.info(f"==== scheme={scheme} (group_key={group_key}) ====")
    # Drop pericytes with no state assignment from the state scheme receivers,
    # but keep them as senders so signaling context is preserved.
    frac_path = args.outdir / f"expressed_fraction_{scheme}.tsv.gz"
    export_expressed_fractions(adata, group_key, frac_path)
    if scheme == "main":  # back-compat name for 02.nichenet.R default
        export_expressed_fractions(adata, group_key,
                                   args.outdir / "expressed_fraction_per_group.tsv.gz")

    diseases = (adata.obs[args.disease_key].cat.categories
                if hasattr(adata.obs[args.disease_key], "cat")
                else adata.obs[args.disease_key].unique())
    all_res = []
    for grp in diseases:
        dsub = adata[adata.obs[args.disease_key] == grp]
        if dsub.n_obs < 200:
            logging.info(f"[{scheme}/{grp}] only {dsub.n_obs} cells; skipping")
            continue
        logging.info(f"[{scheme}/{grp}] running liana on {dsub.n_obs} cells")
        res = run_one(dsub.copy(), group_key, args.expr_prop, args.min_cells, args.seed)
        if res is None:
            continue
        res = res.copy(); res["disease_group"] = grp
        res.to_csv(args.outdir / f"liana_res_{scheme}_{grp}.tsv.gz", sep="\t", index=False)
        all_res.append(res)
    if not all_res:
        logging.error(f"[{scheme}] no liana results produced")
        return
    res = pd.concat(all_res, ignore_index=True)

    receivers = pick_receivers(pd.unique(res["target"]))
    logging.info(f"[{scheme}] receivers: {receivers}")
    into = res[res["target"].isin(receivers)].copy()
    into = into.sort_values(["disease_group", "target", "magnitude_rank"])
    into.to_csv(args.outdir / f"liana_into_receivers_{scheme}.tsv.gz", sep="\t", index=False)
    panel = into[into["ligand_complex"].isin(LIGAND_PANEL)].copy()
    panel.to_csv(args.outdir / f"ligand_panel_into_receivers_{scheme}.tsv.gz",
                 sep="\t", index=False)

    for recv in receivers:
        sub = res[res["target"] == recv]
        if sub.empty:
            continue
        top = (sub.sort_values("magnitude_rank")
                  .drop_duplicates(["source", "ligand_complex", "receptor_complex"])
                  .head(25))
        try:
            p = li.pl.dotplot(
                liana_res=top, colour="magnitude_rank", size="specificity_rank",
                source_labels=list(top["source"].unique()),
                target_labels=[recv], top_n=25, orderby="magnitude_rank",
                orderby_ascending=True, figure_size=(9, 8),
            )
            p.save(fig_dir / f"dotplot_into_{recv}.pdf")
            p.save(fig_dir / f"dotplot_into_{recv}.png", dpi=300)
        except Exception as e:
            logging.warning(f"dotplot for {recv} failed: {e}")


def main():
    args = parse_args()
    configure_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.outdir / "figures"; fig_dir.mkdir(exist_ok=True)

    adata = sc.read_h5ad(args.adata)
    if "logcounts" in adata.layers:
        adata.X = adata.layers["logcounts"]
    logging.info(f"Niche: {adata.shape}")

    for scheme in [s.strip() for s in args.schemes.split(",") if s.strip()]:
        if scheme not in SCHEME_KEYS:
            logging.warning(f"unknown scheme {scheme}; skipping"); continue
        key = SCHEME_KEYS[scheme]
        if key not in adata.obs.columns:
            logging.warning(f"{key} not in obs; skipping scheme {scheme}"); continue
        run_scheme(adata, scheme, key, args, fig_dir)

    session_info.show()


if __name__ == "__main__":
    main()
