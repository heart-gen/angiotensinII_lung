"""
Are the AGTR1-undetected pericytes a distinct population, or just dropout?

AGTR1 is detected in a minority of pericytes. Two readings are possible:

  (i)  BIOLOGY -- the undetected cells are a real AGTR1-negative subpopulation
       that sits somewhere else in the tissue, and
  (ii) TECHNICAL -- AGTR1 is a low-abundance transcript, so at 10x depths most
       cells that express it are still sequenced as zeros.

Reading (i) has a testable consequence: if AGTR1-undetected cells were a distinct
population, their zero rate should be higher than that of equally-abundant genes,
and their tissue position should differ. This script tests both halves.

Outputs (all under --outdir):

  gene_detection_curve.tsv.gz  every gene: pooled abundance vs observed zero rate.
                               Defines the matched-gene reference.
  matched_gene_null.tsv        the K genes matched to AGTR1 on pooled abundance,
                               with their own zero rates -- the null distribution.
  dropout_model_summary.tsv    observed vs expected AGTR1-undetected cells (n and
                               fraction), with the matched-gene interval, an
                               empirical p, and a Poisson reference.
  detection_by_depth.tsv       AGTR1 detection rate by sequencing-depth decile,
                               against the matched-gene expectation in the same
                               cells.
  airspace_effect_ladder.tsv   airspace affinity ~ AGTR1 readout, five model
                               specifications on ONE common complete-case cell set.
  agtr1_cells.tsv.gz           tidy per-cell table backing the figure panels.

Why the matched-gene null and not a parametric one: every parametric dropout model
(Poisson, NB) needs a dispersion assumption that is itself what is in question. The
empirical null asks the narrower and more defensible question -- among genes of the
SAME pooled abundance in the SAME cells, is AGTR1's zero rate unusual? A Poisson
expectation is reported alongside, but as a reference only; real counts are
overdispersed and Poisson under-predicts zeros for essentially every gene, so it is
not a fair test and is labelled as such in the output.

Cell universe: the 11,680 pericytes, identical and identically ordered across
pericyte.hlca_full.dataset.h5ad, pericyte_states.h5ad and the airspace denoising
table (asserted, not assumed).

WHERE THE COUNTS COME FROM. The `counts` layer of these objects is NOT raw: it is
the soupX ambient-corrected matrix (float, minimum 0.00455), and a sampling model
fitted to it is meaningless. True integer counts survive only in `raw/X` of
pericyte.hlca_full.dataset.h5ad, which is what the dropout model uses. The
project's own `AGTR1_detect` is soupX-derived, so both definitions are carried:
panels A and E need raw counts to be internally consistent with the matched genes,
while the airspace ladder keeps the project definition so it stays comparable with
the published estimate. Their agreement is reported in dropout_model_summary.tsv.

Only the pieces actually needed are read (`raw/X`, `raw/var`, `obs`), rather than
whole AnnData objects -- one of these files is 8.4 GB and none of the other layers
are used.
"""
import logging
import argparse
import warnings
import h5py
import numpy as np
import pandas as pd
import session_info
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pathlib import Path
from scipy import sparse
from anndata.io import read_elem

GENE = "AGTR1"

## Matched on log10 pooled abundance. 200 genes is enough for a stable 2.5-97.5
## percentile interval while keeping the abundance window tight (the window is
## reported in the summary, and the caller should check it is narrow).
N_MATCHED = 200

## Genes seen in fewer than this many pericytes have a zero rate that is pinned to
## ~1 by construction and carry no information about dropout; genes seen in all
## cells are equally uninformative at the other end.
MIN_CELLS = 5

DEPTH_BINS = 10


def configure_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--raw-adata", type=Path,
                        default=Path("../../pericyte_analysis/_m/pericyte.hlca_full.dataset.h5ad"),
                        help="Pericyte AnnData whose raw/X holds integer counts")
    parser.add_argument("--state-adata", type=Path,
                        default=Path("../../../pericyte_states/_m/pericyte_states.h5ad"),
                        help="Pericyte AnnData supplying obs (states, AGTR1_detect)")
    parser.add_argument("--qc-adata", type=Path,
                        default=Path("./pericytes_with_airspace_score.h5ad"),
                        help="Airspace AnnData supplying the HLCA QC covariates "
                             "(n_counts/n_genes/pct_mito) used by the published model")
    parser.add_argument("--denoise", type=Path,
                        default=Path("./airspace/pericytes_airspace_denoising.tsv"),
                        help="Long-format scVI denoising table (index x Model)")
    parser.add_argument("--denoise-model", default="Pericyte-only-trained",
                        help="Which scVI model supplies the denoised AGTR1 lens. "
                             "Pericyte-only-trained matches 03.agtr1_lenses.R.")
    parser.add_argument("--outdir", type=Path, default=Path("./agtr1_dropout"))
    parser.add_argument("--n-matched", type=int, default=N_MATCHED)
    return parser.parse_args()


## ---- loading ------------------------------------------------------------
def read_obs(path: Path, columns=None) -> pd.DataFrame:
    """Read only the obs frame of an h5ad."""
    with h5py.File(path, "r") as fh:
        obs = read_elem(fh["obs"])
    return obs if columns is None else obs[[c for c in columns if c in obs]]


def load_raw_counts(path: Path):
    """raw/X as CSC, with symbol var_names.

    Returns (X, var_names, barcodes). var_names are made unique the same way
    scanpy would, so a duplicated symbol cannot silently alias a different gene.
    """
    with h5py.File(path, "r") as fh:
        if "raw" not in fh:
            raise KeyError(f"{path} has no raw/ group -- no integer counts anywhere")
        X = read_elem(fh["raw"]["X"])
        var = read_elem(fh["raw"]["var"])
        barcodes = read_elem(fh["obs"]).index.astype(str)

    X = X.tocsc() if sparse.issparse(X) else sparse.csc_matrix(X)
    sample = X.data[:100000]
    if sample.size and not np.allclose(sample, np.round(sample)):
        raise ValueError(f"raw/X of {path} is not integer-valued")

    names = pd.Index(var["feature_name"].astype(str) if "feature_name" in var
                     else var.index.astype(str))
    dup = names.duplicated(keep=False)
    if dup.any():
        names = names.to_series().reset_index(drop=True)
        counts = names.groupby(names).cumcount()
        names = pd.Index(np.where(counts > 0, names + "-" + counts.astype(str), names))
    logging.info("loaded raw counts: %d cells x %d genes (%d nonzero)",
                 X.shape[0], X.shape[1], X.nnz)
    return X, names, barcodes


## ---- the matched-gene dropout model -------------------------------------
def gene_detection_curve(X, var_names, depth):
    """Pooled abundance and observed zero rate for every gene.

    Abundance is the pooled rate (total counts of the gene / total counts of all
    genes, scaled to counts per 10k), i.e. the maximum-likelihood rate under the
    depth-scaled sampling model. It is NOT the mean of per-cell normalized values,
    which is dominated by shallow cells.
    """
    n_cells = X.shape[0]
    total = np.asarray(X.sum(axis=0)).ravel()
    n_det = np.diff(X.indptr) if sparse.isspmatrix_csc(X) else (X > 0).sum(axis=0)
    n_det = np.asarray(n_det).ravel()
    return pd.DataFrame({
        "gene": np.asarray(var_names),
        "total_counts": total,
        "pooled_cp10k": total / depth.sum() * 1e4,
        "n_detected": n_det,
        "frac_detected": n_det / n_cells,
        "frac_undetected": 1.0 - n_det / n_cells,
    })


def match_genes(curve: pd.DataFrame, gene: str, k: int, min_cells: int):
    """The k genes closest to `gene` in log10 pooled abundance."""
    target = curve.loc[curve["gene"] == gene]
    if len(target) != 1:
        raise KeyError(f"{gene} appears {len(target)} times in the gene curve")
    target = target.iloc[0]
    if target["n_detected"] < min_cells:
        raise ValueError(f"{gene} detected in only {target['n_detected']} cells")

    pool = curve[(curve["gene"] != gene) &
                 (curve["n_detected"] >= min_cells) &
                 (curve["n_detected"] < len(curve))].copy()
    ## Mitochondrial and ribosomal genes are ambient-dominated: their detection
    ## rates reflect soup, not per-cell sampling, so they are not fair matches.
    pool = pool[~pool["gene"].str.match(r"^(MT-|RP[LS]\d|RPLP|RPSA)")]

    lt = np.log10(target["pooled_cp10k"])
    pool["abundance_dist"] = (np.log10(pool["pooled_cp10k"]) - lt).abs()
    matched = pool.nsmallest(k, "abundance_dist").copy()
    matched["target_gene"] = gene
    return target, matched.sort_values("abundance_dist")


def summarise_dropout(target, matched, depth, agtr1_counts):
    """Observed vs expected undetected cells."""
    n_cells = len(depth)
    obs_frac = float(target["frac_undetected"])
    obs_n = int(round(obs_frac * n_cells))

    exp_frac = float(matched["frac_undetected"].median())
    lo, hi = np.percentile(matched["frac_undetected"], [2.5, 97.5])

    ## Two-sided empirical p: how often is a matched gene at least as extreme,
    ## in either direction, as AGTR1? +1 is the usual finite-sample correction.
    d_obs = abs(obs_frac - exp_frac)
    n_extreme = int((np.abs(matched["frac_undetected"] - exp_frac) >= d_obs).sum())
    p_emp = (n_extreme + 1) / (len(matched) + 1)

    sd = float(matched["frac_undetected"].std(ddof=1))
    z = (obs_frac - exp_frac) / sd if sd > 0 else np.nan

    ## Poisson reference: expected zeros if AGTR1 counts were pure depth-scaled
    ## sampling at its own pooled rate. Reported, not tested -- see module docstring.
    lam = float(target["total_counts"]) / float(depth.sum())
    pois_n = float(np.exp(-lam * depth).sum())

    return pd.DataFrame([{
        "gene": target["gene"],
        "n_cells": n_cells,
        "observed_undetected_n": obs_n,
        "observed_undetected_frac": obs_frac,
        "observed_detected_n": n_cells - obs_n,
        "matched_expected_undetected_frac": exp_frac,
        "matched_expected_undetected_n": exp_frac * n_cells,
        "matched_lo_frac": float(lo), "matched_hi_frac": float(hi),
        "matched_lo_n": float(lo) * n_cells, "matched_hi_n": float(hi) * n_cells,
        "matched_sd_frac": sd, "z": z, "p_empirical": p_emp,
        "n_matched_genes": len(matched),
        "pooled_cp10k": float(target["pooled_cp10k"]),
        "matched_cp10k_min": float(matched["pooled_cp10k"].min()),
        "matched_cp10k_max": float(matched["pooled_cp10k"].max()),
        "poisson_reference_undetected_n": pois_n,
        "poisson_reference_undetected_frac": pois_n / n_cells,
        "mean_counts_when_detected": float(
            agtr1_counts[agtr1_counts > 0].mean()) if (agtr1_counts > 0).any() else 0.0,
        "max_counts": float(agtr1_counts.max()),
    }])


def detection_by_depth(X, var_names, matched_genes, gene, depth, n_bins):
    """Observed AGTR1 detection per depth decile, against the matched-gene rate.

    Panel E only earns its place if it says something panel A cannot: A asks
    whether the zero COUNT is unusual, E asks whether the zeros are concentrated
    in the shallow cells (dropout) or spread evenly (biology).
    """
    idx = pd.Index(np.asarray(var_names))
    gcol = np.asarray(X[:, idx.get_loc(gene)].todense()).ravel()
    M = X[:, [idx.get_loc(g) for g in matched_genes]]
    M = (M > 0).astype(np.float32).toarray()

    ## qcut on ranks so ties in depth cannot collapse a bin.
    bins = pd.qcut(pd.Series(depth).rank(method="first"), n_bins, labels=False)
    rows = []
    for b in range(n_bins):
        m = bins.values == b
        n = int(m.sum())
        det = int((gcol[m] > 0).sum())
        rate = det / n
        ## Wilson interval -- the normal approximation is unusable at these rates.
        zc = 1.959964
        cen = (det + zc**2 / 2) / (n + zc**2)
        half = zc / (n + zc**2) * np.sqrt(det * (n - det) / n + zc**2 / 4)
        per_gene = M[m].mean(axis=0)
        rows.append({
            "depth_bin": b + 1, "n_cells": n,
            "depth_median": float(np.median(depth[m])),
            "depth_min": float(depth[m].min()), "depth_max": float(depth[m].max()),
            "n_detected": det, "obs_rate": rate,
            "obs_lo": max(cen - half, 0.0), "obs_hi": min(cen + half, 1.0),
            "matched_rate": float(np.median(per_gene)),
            "matched_lo": float(np.percentile(per_gene, 2.5)),
            "matched_hi": float(np.percentile(per_gene, 97.5)),
        })
    return pd.DataFrame(rows)


## ---- per-cell table and the airspace model ladder -----------------------
def build_cell_table(state_obs, qc_obs, X, var_names, barcodes, denoise,
                     denoise_model, agtr1_counts):
    out = pd.DataFrame({"barcode": barcodes})

    keep = ["donor_id", "dataset", "study", "lung_condition", "sex",
            "pericyte_state", "state_program", "AGTR1_detect", "AGTR1_expr"]
    for k in keep:
        if k in state_obs:
            out[k] = np.asarray(state_obs[k].values)
    out["age"] = pd.to_numeric(state_obs["age_or_mean_of_age_range"],
                               errors="coerce").values

    ## The QC covariates are taken from the HLCA obs rather than recomputed, so the
    ## fully adjusted row of the ladder is fitted on exactly the covariates the
    ## published per-cell model used. Where HLCA does not supply one, fall back to
    ## deriving it from raw counts and say so.
    qc = qc_obs.reindex(pd.Index(barcodes))
    depth_raw = np.asarray(X.sum(axis=1)).ravel()
    if "n_counts" in qc and qc["n_counts"].notna().all():
        out["n_counts"] = qc["n_counts"].to_numpy(dtype=float)
    else:
        logging.warning("n_counts absent from QC obs; using raw library size")
        out["n_counts"] = depth_raw
    if "n_genes" in qc and qc["n_genes"].notna().all():
        out["n_genes"] = qc["n_genes"].to_numpy(dtype=float)
    else:
        logging.warning("n_genes absent from QC obs; deriving from raw counts")
        out["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
    if "pct_mito" in qc and qc["pct_mito"].notna().all():
        pm = qc["pct_mito"].to_numpy(dtype=float)
    else:
        logging.warning("pct_mito absent from QC obs; deriving from raw counts")
        mt = np.array([str(g).startswith("MT-") for g in var_names])
        mito = np.asarray(X[:, np.where(mt)[0]].sum(axis=1)).ravel()
        pm = np.divide(mito, depth_raw, out=np.zeros(len(depth_raw)),
                       where=depth_raw > 0)
    ## HLCA stores pct_mito on 0-100 in some releases; the published model
    ## rescales, so match it.
    out["pct_mito"] = pm / 100.0 if np.nanmax(pm) > 1 else pm

    out["raw_library_size"] = depth_raw
    out["AGTR1_counts"] = agtr1_counts
    ## Raw detection is what the sampling model in panels A/E is about; the
    ## soupX-derived AGTR1_detect that travels with the object is what the rest of
    ## the project has published. Keep both, explicitly named.
    out["AGTR1_detect_raw"] = (agtr1_counts > 0).astype(int)
    if "AGTR1_detect" in out:
        out["AGTR1_detect"] = out["AGTR1_detect"].astype(int)
    else:
        out["AGTR1_detect"] = out["AGTR1_detect_raw"]

    den = denoise[denoise["Model"] == denoise_model][
        ["index", "airspace_score", "AGTR1_scvi"]]
    if den.empty:
        raise KeyError(f"no rows for scVI model {denoise_model!r}")
    merged = out.merge(den, left_on="barcode", right_on="index", how="inner")
    if len(merged) != len(out):
        raise ValueError(f"denoising join lost cells: {len(out)} -> {len(merged)}")
    return merged.drop(columns=["index"])


def _fit_lmm(formula, df, group_col="donor_id"):
    """MixedLM with a donor random intercept and a dataset variance component,
    falling back to donor-clustered OLS if the mixed fit will not converge."""
    groups = df[group_col].astype("category")
    vc = ({"dataset_vc": "0 + C(dataset)"}
          if "dataset" in df and df["dataset"].nunique() > 1 else None)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            model = sm.MixedLM.from_formula(formula, data=df, groups=groups,
                                            re_formula="1", vc_formula=vc)\
                              .fit(method="lbfgs", reml=True, maxiter=3000, disp=False)
        return model, False
    except Exception as exc:                                   # pragma: no cover
        logging.warning("MixedLM failed (%s); clustered OLS fallback", exc)
        model = sm.OLS.from_formula(formula, data=df)\
                      .fit(cov_type="cluster", cov_kwds={"groups": df[group_col]})
        return model, True


def airspace_ladder(cells: pd.DataFrame):
    """Five specifications of `airspace affinity ~ AGTR1 readout`.

    All per-cell specifications are fitted on ONE common complete-case set, so the
    ladder isolates the specification and not the sample. The outcome is z-scored
    within that set and AGTR1_scvi is z-scored too, so every estimate reads as
    "SD of airspace affinity" and the rows are directly comparable.

    The ladder exists because a single adjusted estimate cannot distinguish a real
    but tiny effect from a depth artifact: spec `depth-adjusted` vs `unadjusted`
    shows what the depth covariates do, and `donor-level` shows whether anything
    survives aggregation, where dropout averages out.
    """
    df = cells.dropna(subset=["airspace_score", "AGTR1_detect", "donor_id",
                              "sex", "age"]).copy()
    ## Donors contributing only one AGTR1 class carry no within-donor contrast.
    both = df.groupby("donor_id", observed=True)["AGTR1_detect"].nunique()
    df = df[df["donor_id"].isin(both[both >= 2].index)].copy()
    df["sex"] = df["sex"].astype("category")
    df["dataset"] = df["dataset"].astype("category")

    sd_air = df["airspace_score"].std(ddof=1)
    df["airspace_z"] = (df["airspace_score"] - df["airspace_score"].mean()) / sd_air
    df["AGTR1_scvi_z"] = (df["AGTR1_scvi"] - df["AGTR1_scvi"].mean()) \
        / df["AGTR1_scvi"].std(ddof=1)

    COV = "age + C(sex) + n_counts + n_genes + pct_mito"
    specs = [
        ("binary detectability, unadjusted", "AGTR1_detect",
         "airspace_z ~ AGTR1_detect", "per-cell LMM (1|donor)"),
        ("binary detectability, depth-adjusted", "AGTR1_detect",
         "airspace_z ~ AGTR1_detect + n_counts + n_genes", "per-cell LMM (1|donor)"),
        ("binary detectability, fully adjusted", "AGTR1_detect",
         f"airspace_z ~ AGTR1_detect + {COV}", "per-cell LMM (1|donor)"),
        ("scVI-denoised AGTR1 (per SD)", "AGTR1_scvi_z",
         f"airspace_z ~ AGTR1_scvi_z + {COV}", "per-cell LMM (1|donor)"),
    ]

    rows = []
    for label, term, formula, spec in specs:
        model, fallback = _fit_lmm(formula, df)
        ci = model.conf_int()
        rows.append({
            "label": label, "term": term, "spec": spec, "formula": formula,
            "estimate": float(model.params[term]), "se": float(model.bse[term]),
            "ci_lower": float(ci.loc[term, 0]), "ci_upper": float(ci.loc[term, 1]),
            "pval": float(model.pvalues[term]),
            "n_cells": int(len(df)), "n_donors": int(df["donor_id"].nunique()),
            "fallback_ols": fallback,
        })

    ## Donor level: the same question after aggregation. Its predictor is a
    ## fraction on 0-1, so its estimate is the SD change between a donor with no
    ## detectable AGTR1 and one where every pericyte is detectable -- the largest
    ## contrast on the ladder, which is why a null here is informative.
    donor = df.groupby("donor_id", observed=True).agg(
        airspace_z=("airspace_z", "mean"),
        frac_AGTR1_pos=("AGTR1_detect", "mean"),
        age=("age", "first"), sex=("sex", "first"), n_cells=("AGTR1_detect", "size"),
    ).reset_index()
    ols = smf.ols("airspace_z ~ frac_AGTR1_pos + age + C(sex)", data=donor).fit()
    ci = ols.conf_int()
    rows.append({
        "label": "donor fraction AGTR1-detectable", "term": "frac_AGTR1_pos",
        "spec": "donor-level OLS",
        "formula": "airspace_z ~ frac_AGTR1_pos + age + C(sex)",
        "estimate": float(ols.params["frac_AGTR1_pos"]),
        "se": float(ols.bse["frac_AGTR1_pos"]),
        "ci_lower": float(ci.loc["frac_AGTR1_pos", 0]),
        "ci_upper": float(ci.loc["frac_AGTR1_pos", 1]),
        "pval": float(ols.pvalues["frac_AGTR1_pos"]),
        "n_cells": int(donor["n_cells"].sum()), "n_donors": int(len(donor)),
        "fallback_ols": False,
    })

    out = pd.DataFrame(rows)
    out["airspace_sd"] = sd_air
    return out, df


def main():
    configure_logging()
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    X, var_names, barcodes = load_raw_counts(args.raw_adata)
    state_obs = read_obs(args.state_adata)
    if not np.array_equal(state_obs.index.astype(str).to_numpy(), barcodes):
        raise ValueError("raw counts and state obs are not the same cells in the "
                         "same order -- refusing to align by position")
    qc_obs = read_obs(args.qc_adata, ["n_counts", "n_genes", "pct_mito"])

    depth = np.asarray(X.sum(axis=1)).ravel()
    curve = gene_detection_curve(X, var_names, depth)
    curve.to_csv(args.outdir / "gene_detection_curve.tsv.gz", sep="\t",
                 index=False, compression="gzip")

    target, matched = match_genes(curve, GENE, args.n_matched, MIN_CELLS)
    matched.to_csv(args.outdir / "matched_gene_null.tsv", sep="\t", index=False)

    agtr1_counts = np.asarray(X[:, var_names.get_loc(GENE)].todense()).ravel()
    summary = summarise_dropout(target, matched, depth, agtr1_counts)
    ## How far the soupX-derived project definition of detection sits from raw
    ## detection. The figure legend quotes this, so it is computed here rather
    ## than inferred later.
    if "AGTR1_detect" in state_obs:
        proj = state_obs["AGTR1_detect"].to_numpy(dtype=int)
        summary["detect_agreement_frac"] = float(
            (proj == (agtr1_counts > 0).astype(int)).mean())
        summary["project_detected_n"] = int(proj.sum())
    summary.to_csv(args.outdir / "dropout_model_summary.tsv", sep="\t", index=False)
    logging.info("observed undetected %d/%d (%.3f); matched expectation %.3f "
                 "[%.3f, %.3f], empirical p = %.3f",
                 summary.observed_undetected_n[0], len(depth),
                 summary.observed_undetected_frac[0],
                 summary.matched_expected_undetected_frac[0],
                 summary.matched_lo_frac[0], summary.matched_hi_frac[0],
                 summary.p_empirical[0])

    depth_tab = detection_by_depth(X, var_names, matched["gene"].tolist(),
                                   GENE, depth, DEPTH_BINS)
    depth_tab.to_csv(args.outdir / "detection_by_depth.tsv", sep="\t", index=False)

    denoise = pd.read_csv(args.denoise, sep="\t")
    cells = build_cell_table(state_obs, qc_obs, X, var_names, barcodes, denoise,
                             args.denoise_model, agtr1_counts)
    cells.to_csv(args.outdir / "agtr1_cells.tsv.gz", sep="\t", index=False,
                 compression="gzip")

    ladder, fitted = airspace_ladder(cells)
    ladder.to_csv(args.outdir / "airspace_effect_ladder.tsv", sep="\t", index=False)
    logging.info("airspace ladder on %d cells / %d donors",
                 ladder.n_cells[0], ladder.n_donors[0])
    for _, r in ladder.iterrows():
        logging.info("  %-40s beta = %+.4f  p = %.3g", r.label, r.estimate, r.pval)

    session_info.show()


if __name__ == "__main__":
    main()
