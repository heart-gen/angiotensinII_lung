"""
QC of the airspace-scored HLCA object, and production of the scVI-denoised AGTR1
lens (`AGTR1_scvi`) that the project uses as its dropout-robust readout.

WHICH COUNTS scVI IS TRAINED ON. The `counts` layer of
pericytes_with_airspace_score.h5ad is byte-identical to its `soupX` layer: the
ambient-corrected float matrix, minimum nonzero 1.86e-4, with only 1.9% of the
stored values integer-valued. A negative-binomial likelihood fitted to that
matrix cannot separate a dropout zero from low expression, which is exactly the
distinction the denoised lens exists to make. Training on it produced a variable
with no relationship to observed AGTR1 (donor-level rho = 0.014, p = 0.91). True
integer counts survive in `raw/X` of airspace.hlca_full.dataset.h5ad, which has
the same cells in the same order and the same genes in the same order (keyed on
var['ensembl_id']); that is what --raw-adata supplies and what is used here.

VALIDITY GATE. The denoised lens is not written to its canonical path unless it
passes validate_denoising(): pooling to donors with >= --min-donor-cells
pericytes averages dropout away, so donor-mean denoised AGTR1 must track the
donor AGTR1 detection rate. Thresholds are pre-specified in that function. A
second, independent control uses external biology rather than the same counts --
AGTR1 marks the mural compartment, so the perivascular-trained model must rank
pericytes above capillary endothelium.
"""
import h5py
import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from scipy import sparse
from anndata import AnnData
from anndata.io import read_elem
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score

SEED = 13

## Pre-specified acceptance thresholds for the denoised lens. Fixed before the
## run so the gate cannot be tuned to whatever the model happens to produce.
PASS_RHO = 0.50
WEAK_RHO = 0.25
MIN_DONOR_CELLS = 20

## Studies contributing fewer than this many cells to a training subset cannot
## support their own batch term: seurat_v3 fits a per-batch loess, which is
## singular at n = 2-6, and scVI would spend a free parameter on a batch it can
## never estimate. They are pooled into one residual batch instead of dropped,
## so no cell loses its denoised value.
MIN_BATCH_CELLS = 20
SMALL_BATCH_LABEL = "_small_studies"

sns.set_context("talk")
sns.set_style("whitegrid")

def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )


def parse_args():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--adata", type=Path,
                        default=Path("./pericytes_with_airspace_score.h5ad"),
                        help="HLCA airspace scored AnnData")
    parser.add_argument("--raw-adata", type=Path,
                        default=Path("./airspace.hlca_full.dataset.h5ad"),
                        help="Pre-QC AnnData whose raw/X holds INTEGER counts. The "
                             "scored object only carries the soupX float matrix, "
                             "which an NB likelihood cannot be fitted to.")
    parser.add_argument("--outdir", type=Path, default=Path("./"))
    parser.add_argument("--n-hvg", type=int, default=4000,
                        help="Genes retained for scVI training (AGTR1 always kept)")
    parser.add_argument("--max-epochs", type=int, default=400,
                        help="Upper bound for the pericyte-only model; early "
                             "stopping decides the actual count. Held at 400 so "
                             "the canonical downstream lens stays reproducible.")
    ## The perivascular model ran 400/400 with early stopping NEVER firing, so
    ## it was reported as undertrained. Two things were wrong, and raising the
    ## ceiling alone fixes neither:
    ##   (a) the ceiling was genuinely below where the run would have stopped;
    ##   (b) scvi-tools defaults early_stopping_min_delta to 0.0, so ANY
    ##       improvement in elbo_validation -- including drift in the 4th
    ##       decimal -- resets the patience counter. On 123k cells that counter
    ##       effectively never expires, so "early stopping did not fire" carried
    ##       no information about convergence. Training loss was already flat
    ##       (659 at epoch 300 -> 662 at epoch 400, i.e. no longer improving).
    ## Both are per-model so the pericyte arm is provably untouched.
    parser.add_argument("--max-epochs-perivascular", type=int, default=1500,
                        help="Upper bound for the perivascular model")
    parser.add_argument("--min-delta-perivascular", type=float, default=1.0,
                        help="Minimum elbo_validation improvement that counts as "
                             "progress for the perivascular model. 0 (the "
                             "scvi-tools default, kept for the pericyte model) "
                             "makes early stopping unable to fire on a large set.")
    parser.add_argument("--min-donor-cells", type=int, default=MIN_DONOR_CELLS,
                        help="Donors with at least this many pericytes define the "
                             "dropout-averaged reference used by the validity gate")
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def load_raw_counts(raw_path: Path, adata: AnnData, layer="counts_raw"):
    """Install true integer counts from raw/X as adata.layers[layer].

    Alignment is asserted rather than assumed: cells must match by obs index in
    the same order, and genes by var['ensembl_id'] against the raw var index
    (the raw object is keyed on Ensembl IDs, the scored one on symbols, but the
    gene order is identical).
    """
    logging.info("Loading integer counts from %s", raw_path)
    with h5py.File(raw_path, "r") as fh:
        if "raw" not in fh:
            raise KeyError(f"{raw_path} has no raw/ group; no integer counts available")
        obs_raw = read_elem(fh["obs"])
        var_raw = read_elem(fh["raw"]["var"])
        if not obs_raw.index.equals(adata.obs.index):
            raise ValueError("raw object cells do not match the scored object "
                             "(index or order differs)")
        if "ensembl_id" not in adata.var:
            raise KeyError("scored var lacks 'ensembl_id'; cannot align to raw var")
        if not np.array_equal(adata.var["ensembl_id"].astype(str).values,
                              var_raw.index.astype(str).values):
            raise ValueError("raw object genes do not match the scored object "
                             "(ensembl_id order differs)")
        X = sparse.csr_matrix(read_elem(fh["raw"]["X"]))

    dat = X.data
    if dat.size and not np.allclose(dat, np.round(dat)):
        raise ValueError("raw/X is not integer-valued; refusing to fit an NB "
                         "likelihood to it (this is the original defect)")
    logging.info("integer counts installed: %d nonzero, min %.4g, max %.4g",
                 dat.size, dat.min() if dat.size else 0, dat.max() if dat.size else 0)
    adata.layers[layer] = X
    return adata


def save_figure(fig, base_path: Path):
    fig.savefig(base_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def define_masks(adata: AnnData):
    logging.info("Define masks")

    pericyte_mask = adata.obs["subclusters"].eq("Pericytes")
    ec_types = [
        "EC general capillary", "EC aerocyte capillary",
        "EC arterial", "EC venous",
    ]
    lymphatic_types = [
        "Lymphatic EC differentiating", "Lymphatic EC mature",
        "Lymphatic EC proliferating"
    ]
    vsmc_types = [
        "Smooth muscle", "Smooth muscle FAM83D+", "SM activated stress response"
    ]
    perivascular_mask = (
        pericyte_mask |
        adata.obs["subclusters"].isin(ec_types + lymphatic_types + vsmc_types)
    )
    return pericyte_mask, perivascular_mask


def denoise_agtr1_scvi(adata, batch_key="study", layer_raw="counts_raw", gene="AGTR1",
                       max_epochs=400, min_delta=0.0, seed=SEED, label="model"):
    """Fit scVI and return the posterior mean AGTR1 rate per cell.

    The layer MUST be integer counts -- asserted, not assumed. Empty and
    near-empty batch categories are dropped/flagged first: a subset inherits the
    parent object's full `study` categories, so scVI would otherwise allocate
    batch terms to categories with zero cells.
    """
    logging.info("[%s] training scVI on %d cells x %d genes", label,
                 adata.n_obs, adata.n_vars)
    scvi.settings.seed = seed

    if layer_raw not in adata.layers:
        raise KeyError(f"layer '{layer_raw}' absent; integer counts are required")
    dat = adata.layers[layer_raw].data
    if dat.size and not np.allclose(dat, np.round(dat)):
        raise ValueError(f"layer '{layer_raw}' is not integer-valued; an NB "
                         "likelihood must not be fitted to it")

    use_batch = batch_key if batch_key in adata.obs else None
    if use_batch is not None:
        col = adata.obs[use_batch]
        if str(col.dtype) == "category":
            adata.obs[use_batch] = col.cat.remove_unused_categories()
        counts = adata.obs[use_batch].value_counts()
        tiny = counts[counts < 3]
        if len(tiny):
            logging.warning("[%s] batches with <3 cells carry unstable batch terms: %s",
                            label, tiny.to_dict())
        logging.info("[%s] %d batches, smallest %d cells", label,
                     adata.obs[use_batch].nunique(), int(counts.min()))

    scvi.model.SCVI.setup_anndata(adata, layer=layer_raw, batch_key=use_batch)
    model = scvi.model.SCVI(adata, gene_likelihood="nb")
    ## min_delta is the size of an elbo_validation improvement that still counts
    ## as progress. At the scvi-tools default of 0.0 the patience counter is
    ## reset by arbitrarily small drift, so a run can exhaust max_epochs while
    ## being flat -- "early stopping did not fire" then means nothing.
    logging.info("[%s] early stopping: monitor=elbo_validation patience=20 "
                 "min_delta=%.4g max_epochs=%d", label, min_delta, max_epochs)
    model.train(max_epochs=max_epochs, early_stopping=True,
                early_stopping_patience=20, early_stopping_min_delta=min_delta,
                plan_kwargs={"lr": 1e-3})

    ran = int(model.trainer.current_epoch)
    stopped = ran < max_epochs
    if not stopped:
        logging.warning("[%s] early stopping NEVER FIRED at max_epochs=%d -- "
                        "convergence is NOT established for this model",
                        label, max_epochs)
    else:
        logging.info("[%s] early stopping fired at epoch %d/%d", label, ran, max_epochs)

    ## Record the loss trajectory so convergence is auditable from the outputs
    ## rather than only from a progress bar buried in the job log.
    hist = {}
    for key in ("elbo_validation", "elbo_train"):
        h = model.history.get(key)
        if h is None or not len(h):
            continue
        v = h.iloc[:, 0].to_numpy(dtype=float)
        hist[f"{key}_first"] = float(v[0])
        hist[f"{key}_final"] = float(v[-1])
        hist[f"{key}_best"] = float(np.nanmin(v))
        ## Improvement over the last 20 recorded epochs -- the same window the
        ## patience counter uses. Near zero means the fit has plateaued.
        tail = v[-21:]
        hist[f"{key}_delta_last20"] = float(tail[0] - tail[-1]) if len(tail) > 1 else float("nan")

    ## Only AGTR1 is needed; materialising all genes would be a dense
    ## n_cells x n_genes frame.
    den = model.get_normalized_expression(gene_list=[gene], library_size=1e4,
                                          return_mean=True)
    adata.obs[f"{gene}_scvi"] = den[gene].values

    info = {"model": label, "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
            "epochs_ran": ran, "max_epochs": int(max_epochs),
            "early_stopping_min_delta": float(min_delta),
            "early_stopped": bool(stopped),
            "n_batches": int(adata.obs[use_batch].nunique()) if use_batch else 0}
    info.update(hist)
    return adata, info


def subset_anndata(adata: AnnData, mask=None, label="perivascular", gene="AGTR1",
                   n_hvg=4000, max_epochs=400, min_delta=0.0, seed=SEED):
    logging.info("Subset data for training: %s (%d cells)", label, int(mask.sum()))

    adata_sub = adata[mask].copy()

    ## Genes seen in almost no cell of this subset have no usable variance
    ## estimate and are what makes the HVG loess singular.
    n_before = adata_sub.n_vars
    keep = np.asarray((adata_sub.layers["counts_raw"] > 0).sum(axis=0)).ravel() >= 3
    if gene in adata_sub.var_names:
        keep[adata_sub.var_names.get_loc(gene)] = True
    adata_sub = adata_sub[:, keep].copy()
    logging.info("[%s] genes detected in >=3 cells: %d of %d", label,
                 adata_sub.n_vars, n_before)

    ## Pool batches too small to model (see MIN_BATCH_CELLS).
    counts = adata_sub.obs["study"].value_counts()
    small = counts[(counts > 0) & (counts < MIN_BATCH_CELLS)].index.tolist()
    batch_col = "study_model"
    adata_sub.obs[batch_col] = adata_sub.obs["study"].astype(str)
    if small:
        adata_sub.obs.loc[adata_sub.obs["study"].isin(small), batch_col] = SMALL_BATCH_LABEL
        logging.warning("[%s] pooling %d studies with <%d cells into '%s': %s",
                        label, len(small), MIN_BATCH_CELLS, SMALL_BATCH_LABEL,
                        {s: int(counts[s]) for s in small})
    adata_sub.obs[batch_col] = adata_sub.obs[batch_col].astype("category")

    ## seurat_v3 is the count-based HVG flavour, which is the correct pairing for
    ## the integer layer scVI is trained on.
    try:
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg, batch_key=batch_col,
                                    layer="counts_raw", flavor="seurat_v3")
    except ValueError as err:
        logging.warning("[%s] batch-aware HVG failed (%s); falling back to a "
                        "pooled HVG fit across all cells", label, err)
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=n_hvg,
                                    layer="counts_raw", flavor="seurat_v3")
    if gene in adata_sub.var_names:
        adata_sub.var.loc[gene, "highly_variable"] = True

    ## APPLY the selection. Previously this flag was computed and then never
    ## used, so scVI trained on all 55,329 genes.
    adata_sub = adata_sub[:, adata_sub.var["highly_variable"].values].copy()
    if gene not in adata_sub.var_names:
        raise ValueError(f"{gene} was dropped by HVG selection")
    logging.info("[%s] training gene set: %d genes (%s retained)", label,
                 adata_sub.n_vars, gene)

    adata_sub, info = denoise_agtr1_scvi(
        adata_sub, batch_key=batch_col, layer_raw="counts_raw", gene=gene,
        max_epochs=max_epochs, min_delta=min_delta, seed=seed, label=label
    )
    info["n_studies_pooled"] = len(small)
    adata.obs.loc[mask, f"{gene}_scvi_{label}"] = adata_sub.obs[f"{gene}_scvi"]
    return adata, info


def _verdict(rho, pval):
    if not np.isfinite(rho) or not np.isfinite(pval):
        return "FAIL"
    if rho >= PASS_RHO and pval < 0.05:
        return "PASS"
    if rho >= WEAK_RHO:
        return "WEAK"
    return "FAIL"


def validate_denoising(adata, pericyte_mask, model_cols, min_cells=MIN_DONOR_CELLS,
                       gene="AGTR1"):
    """Pre-specified validity gate for the denoised AGTR1 lens.

    A denoiser is allowed to smooth, shrink and compress, but once dropout has
    been averaged away by pooling it must still agree with the observed signal.
    Donors contributing at least `min_cells` pericytes give a dropout-averaged
    estimate of that donor's true AGTR1 level, so donor-mean denoised AGTR1 has
    to track the donor AGTR1 detection rate. That is the primary criterion, and
    its thresholds (PASS_RHO / WEAK_RHO) are fixed at module level before the
    run so the gate cannot be tuned to the model's output.

        PASS   rho >= 0.50 and p < 0.05
        WEAK   0.25 <= rho < 0.50
        FAIL   rho < 0.25 or p >= 0.05

    Three supporting diagnostics are reported alongside, because each of them
    distinguishes a specific failure mode:

      * cell-level agreement with observed AGTR1 -- near zero means the output
        carries none of the signal it claims to denoise;
      * between-study variance fraction -- a denoised value dominated by batch
        is a study label, not an expression estimate (raw AGTR1 sits near 0.08);
      * AUROC for pericytes vs capillary endothelium (perivascular model only) --
        an independent control that uses external biology rather than the same
        counts, since AGTR1 marks the mural compartment in this tissue.
    """
    logging.info("Validating the denoised lens against pooled detection")
    obs = adata.obs
    per = obs.loc[pericyte_mask, ["donor_id", "study", f"{gene}_detect",
                                  f"{gene}_expr"] + model_cols].copy()
    per[f"{gene}_detect"] = per[f"{gene}_detect"].astype(float)

    rows = []
    for col in model_cols:
        d = per[[col, "donor_id", "study", f"{gene}_detect", f"{gene}_expr"]].dropna()

        ## cell level -- diagnostic, not the criterion
        for ref in (f"{gene}_detect", f"{gene}_expr"):
            r = spearmanr(d[col], d[ref])
            rows.append(dict(model=col, level="cell", reference=ref,
                             statistic="spearman_rho", value=r.correlation,
                             p_value=r.pvalue, n=len(d), verdict=""))

        ## donor level -- the criterion, plus a sensitivity ladder
        g = d.groupby("donor_id", observed=True).agg(
            n=(col, "size"), den=(col, "mean"), det=(f"{gene}_detect", "mean"))
        for cut in (5, min_cells, 50, 100):
            gg = g[g["n"] >= cut]
            if len(gg) < 5:
                rows.append(dict(model=col, level=f"donor>={cut}",
                                 reference=f"{gene}_detect_rate",
                                 statistic="spearman_rho", value=np.nan,
                                 p_value=np.nan, n=len(gg), verdict="n/a"))
                continue
            r = spearmanr(gg["det"], gg["den"])
            rows.append(dict(model=col, level=f"donor>={cut}",
                             reference=f"{gene}_detect_rate",
                             statistic="spearman_rho", value=r.correlation,
                             p_value=r.pvalue, n=len(gg),
                             verdict=_verdict(r.correlation, r.pvalue)
                                     if cut == min_cells else ""))

        ## how much of the denoised value is just the study label
        ss_tot = float(((d[col] - d[col].mean()) ** 2).sum())
        ss_bet = float(d.groupby("study", observed=True)[col].apply(
            lambda s: len(s) * (s.mean() - d[col].mean()) ** 2).sum())
        rows.append(dict(model=col, level="cell", reference="study",
                         statistic="between_study_variance_fraction",
                         value=ss_bet / ss_tot if ss_tot > 0 else np.nan,
                         p_value=np.nan, n=len(d), verdict=""))

    ## independent positive control: does the lens know pericytes from capillary EC?
    ec = ["EC general capillary", "EC aerocyte capillary"]
    ctl = obs.loc[pericyte_mask | obs["subclusters"].isin(ec)].copy()
    ctl["is_pericyte"] = ctl["subclusters"].eq("Pericytes").astype(int)
    for col in model_cols + [f"{gene}_detect"]:
        if col not in ctl or ctl[col].notna().sum() == 0:
            continue
        sub_ = ctl[[col, "is_pericyte"]].dropna()
        if sub_["is_pericyte"].nunique() < 2:
            continue
        auc = roc_auc_score(sub_["is_pericyte"], sub_[col].astype(float))
        rows.append(dict(model=col, level="pericyte_vs_capillary_EC",
                         reference="subclusters", statistic="AUROC", value=auc,
                         p_value=np.nan, n=len(sub_),
                         verdict="PASS" if auc > 0.5 else "FAIL"))

    res = pd.DataFrame(rows)

    ## model-free display: observed detection rate across denoised quintiles
    quint = []
    for col in model_cols:
        d = per[[col, f"{gene}_detect"]].dropna()
        if d[col].nunique() < 5:
            continue
        ## A heavily tied variable collapses bins, so let qcut choose how many
        ## survive and label afterwards rather than fixing the label count.
        binned = pd.qcut(d[col], 5, duplicates="drop")
        codes = binned.cat.codes
        d = d.assign(q=pd.Categorical(
            [f"Q{c + 1}" for c in codes],
            categories=[f"Q{i + 1}" for i in range(len(binned.cat.categories))],
            ordered=True))
        agg = (d.groupby("q", observed=True)
                 .agg(n=(col, "size"), mean_denoised=(col, "mean"),
                      observed_detect_rate=(f"{gene}_detect", "mean"))
                 .reset_index())
        agg.insert(0, "model", col)
        quint.append(agg)
    quint = pd.concat(quint, ignore_index=True) if quint else pd.DataFrame()

    primary = res[(res["level"] == f"donor>={min_cells}") & (res["verdict"] != "")]
    verdicts = set(primary["verdict"])
    overall = "PASS" if verdicts and verdicts <= {"PASS"} else (
        "WEAK" if "PASS" in verdicts or "WEAK" in verdicts else "FAIL")
    logging.info("Validation verdict: %s (%s)", overall,
                 primary.set_index("model")["verdict"].to_dict())
    return res, quint, overall


def write_validation(base: Path, res, quint, overall, min_cells):
    res.to_csv(base.with_suffix(".tsv"), sep="\t", index=False)
    if len(quint):
        quint.to_csv(base.with_name(base.name + "_quintiles").with_suffix(".tsv"),
                     sep="\t", index=False)
    lines = [
        "AGTR1 scVI denoising -- validity gate",
        "=" * 60,
        f"OVERALL VERDICT: {overall}",
        "",
        "Primary criterion: Spearman correlation between donor-mean denoised",
        f"AGTR1 and the donor AGTR1 detection rate, donors with >= {min_cells}",
        f"pericytes. Pre-specified: PASS rho >= {PASS_RHO} and p < 0.05;",
        f"WEAK {WEAK_RHO} <= rho < {PASS_RHO}; FAIL otherwise.",
        "",
        "Pooling to that level averages dropout away, so the detection rate is a",
        "good estimate of the donor's true AGTR1 level. A denoised lens that does",
        "not track it is not measuring AGTR1, and a null result obtained with it",
        "is uninformative rather than evidence of absence.",
        "",
        res.to_string(index=False),
    ]
    if len(quint):
        lines += ["", "Observed detection rate by denoised quintile (model-free):",
                  quint.to_string(index=False)]
    base.with_suffix(".txt").write_text("\n".join(lines) + "\n")
    logging.info("Validation report written to %s", base.with_suffix(".txt"))


def compare_predicted_agtr1(adata, method_col="AGTR1_scvi", pericyte_mask=None):
    logging.info("Compare airspace scores after denoising")

    if pericyte_mask is None:
        pericyte_mask = adata.obs["subclusters"].eq("Pericytes")

    df = adata.obs.loc[pericyte_mask, ["AGTR1_detect", "airspace_score", "donor_id",
                                       method_col]].dropna()
    df["AGTR1_detect"] = df["AGTR1_detect"].astype(int)

    # Calibrate threshold by matching prevalence among pericytes
    prev = df["AGTR1_detect"].mean()
    thr  = df[method_col].quantile(1 - prev)
    df["AGTR1_predpos"] = (df[method_col] >= thr).astype(int)

    # Compare airspace (with clustered SE)
    res = smf.ols("airspace_score ~ AGTR1_predpos", data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df["donor_id"]}
    )

    return {
        "threshold": float(thr),
        "coef": float(res.params["AGTR1_predpos"]),
        "pval": float(res.pvalues["AGTR1_predpos"]),
        "n": int(df.shape[0])
    }


def write_summary_results(
    filepath: Path, correlation_r: float, comparison_results: dict,
    header: str = "AGTR1 scVI Model Evaluation Summary"
):
    """Write correlation and prediction comparison results to a summary text file."""
    logging.info("Writing summary of results to file")

    lines = [
        f"{header}\n",
        "Pericyte AGTR1_scvi correlation (perivascular vs pericyte-only): "
        f"{correlation_r:.4f}\n",
        "\nComparison Results (Perivascular model on Pericytes):",
        f"  Threshold: {comparison_results['threshold']:.4f}",
        f"  Coefficient: {comparison_results['coef']:.4f}",
        f"  p-value: {comparison_results['pval']:.4g}",
        f"  N pericytes: {comparison_results['n']}",
        ""
    ]

    filepath = filepath.with_suffix(".txt")
    with open(filepath, "w") as f:
        f.write("\n".join(lines))

    print(f"Summary results saved to: {filepath}")


def plot_corr(adata: AnnData, pericyte_mask=None, base: Path="./", tsv_path: Path=None):
    logging.info("Plot scatter of correlation faceted by model")

    if pericyte_mask is None:
        pericyte_mask = adata.obs["subclusters"] == "Pericytes"

    df = adata.obs.loc[pericyte_mask, [
        "donor_id", "airspace_score", "AGTR1_scvi_perivascular", "AGTR1_scvi_pericytes"
    ]].dropna()

    # Melt for seaborn
    df_melted = df.reset_index().melt(
        id_vars=["index", "donor_id", "airspace_score"],
        value_vars=["AGTR1_scvi_perivascular", "AGTR1_scvi_pericytes"],
        var_name="Model", value_name="AGTR1_scvi"
    )
    df_melted["Model"] = df_melted["Model"].map({
        "AGTR1_scvi_perivascular": "Perivascular-trained",
        "AGTR1_scvi_pericytes": "Pericyte-only-trained"
    })
    tsv_path = base.with_suffix(".tsv") if tsv_path is None else tsv_path
    df_melted.to_csv(tsv_path, sep="\t", index=False)
    print(f"Saved melted correlation data to: {tsv_path}")

    # Compute correlation for annotations (least 2 observations)
    donor_means = (
        df_melted.groupby(["Model", "donor_id"], observed=False)
        [["AGTR1_scvi", "airspace_score"]].mean().reset_index()
    )
    
    donor_corrs = (
        donor_means.groupby("Model", observed=False)
        .apply(lambda g: pd.Series(pearsonr(g["AGTR1_scvi"], g["airspace_score"]),
                                   index=["r", "pval"])).reset_index()
    )
    
    # Plot
    g = sns.lmplot(
        data=donor_means, x="AGTR1_scvi", y="airspace_score", col="Model",
        col_wrap=2, height=5, aspect=1,  ci=None,
        scatter_kws={"s": 30, "alpha": 0.7},
        line_kws={"linewidth": 1.5}, legend_out=True
    )

    # Annotate correlations
    for ax, (_, row) in zip(g.axes.flat, donor_corrs.iterrows()):
        r_txt = f"$r = {row['r']:.2f}$\n$p = {row['pval']:.1e}$"
        ax.text(0.05, 0.95, r_txt, transform=ax.transAxes,
                verticalalignment='top', fontsize=12)

    g.set_axis_labels("AGTR1 (scVI-denoised)", "Airspace Score")
    plt.tight_layout()
    save_figure(g.fig, base)


def main():
    # Load parser and logging
    args = parse_args()
    configure_logging()

    # Set parameters
    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    # Load data
    adata = sc.read_h5ad(args.adata)

    # Integer counts -- scVI cannot be fitted to the soupX float matrix
    adata = load_raw_counts(args.raw_adata, adata)

    # Define masks
    pericyte_mask, perivascular_mask = define_masks(adata)

    train_info = []

    # Perivascular model
    adata, info = subset_anndata(adata, perivascular_mask, "perivascular",
                                 n_hvg=args.n_hvg,
                                 max_epochs=args.max_epochs_perivascular,
                                 min_delta=args.min_delta_perivascular,
                                 seed=args.seed)
    train_info.append(info)

    # Pericyte-only model
    ## Deliberately unchanged: 400 epochs, min_delta 0. This arm early-stopped
    ## at 152/400 and is the canonical lens every downstream module reads, so
    ## its settings and seed are held fixed to keep it bit-reproducible.
    adata, info = subset_anndata(adata, pericyte_mask, "pericytes",
                                 n_hvg=args.n_hvg, max_epochs=args.max_epochs,
                                 min_delta=0.0, seed=args.seed)
    train_info.append(info)

    # Compute correlation
    r = np.corrcoef(
        adata.obs.loc[pericyte_mask, "AGTR1_scvi_perivascular"],
        adata.obs.loc[pericyte_mask, "AGTR1_scvi_pericytes"],
    )[0, 1]

    # Evalute model
    results = compare_predicted_agtr1(
        adata, method_col="AGTR1_scvi_perivascular",
        pericyte_mask=pericyte_mask
    )

    # Log correlation
    airspace_dir = outdir / "airspace"
    airspace_dir.mkdir(exist_ok=True)

    write_summary_results(
        airspace_dir / "agtr1_model_summary", r, results
    )
    pd.DataFrame(train_info).to_csv(
        airspace_dir / "agtr1_scvi_training_info.tsv", sep="\t", index=False)

    # ---- validity gate -------------------------------------------------
    ## The lens only reaches its canonical path if it demonstrably measures
    ## AGTR1. Otherwise it is written beside it as a CANDIDATE so nothing
    ## downstream silently consumes an unvalidated variable.
    model_cols = ["AGTR1_scvi_perivascular", "AGTR1_scvi_pericytes"]
    res, quint, overall = validate_denoising(
        adata, pericyte_mask, model_cols, min_cells=args.min_donor_cells
    )
    write_validation(airspace_dir / "agtr1_denoising_validation", res, quint,
                     overall, args.min_donor_cells)

    base = airspace_dir / "pericytes_airspace_denoising"
    canonical = base.with_suffix(".tsv")
    if overall == "PASS":
        if canonical.exists():
            backup = canonical.with_suffix(".superseded.tsv")
            canonical.replace(backup)
            logging.info("previous lens preserved at %s", backup)
        tsv_path = canonical
    else:
        tsv_path = base.with_name(base.name + ".CANDIDATE").with_suffix(".tsv")
        logging.error("VALIDATION %s -- writing to %s and leaving the canonical "
                      "table untouched. Do NOT re-run downstream analyses.",
                      overall, tsv_path)

    # Plot scatter
    plot_corr(adata, pericyte_mask, base, tsv_path=tsv_path)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
