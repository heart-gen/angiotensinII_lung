import scvi
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from anndata import AnnData
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score

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
    parser.add_argument("--outdir", type=Path, default=Path("./"))
    return parser.parse_args()


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


def denoise_agtr1_scvi(adata, batch_key="study", layer_raw=None, gene="AGTR1"):
    logging.info("Training scVI model on perivascular cells")

    use_layer = layer_raw if layer_raw in adata.layers else None
    if use_layer is not None:
        X_bak = adata.X
        adata.X = adata.layers[layer_raw]

    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key if batch_key in adata.obs else None)
    model = scvi.model.SCVI(adata, gene_likelihood="nb")
    model.train(max_epochs=60, early_stopping=True, plan_kwargs={"lr": 1e-3})

    den = model.get_normalized_expression(library_size=1e4, return_mean=True)
    adata.layers["scvi_denoised"] = den.values
    adata.obs[f"{gene}_scvi"] = den[gene].values

    if use_layer is not None:
        adata.X = X_bak

    return adata


def subset_anndata(adata: AnnData, mask=None, label="perivascular", gene="AGTR1"):
    logging.info("Subset data for training: %s", label)

    adata_sub = adata[mask].copy()
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=4000, batch_key="study")
    if gene in adata_sub.var_names:
        adata_sub.var.loc[gene, "highly_variable"] = True

    adata_sub = denoise_agtr1_scvi(
        adata_sub, batch_key="study", layer_raw="counts", gene=gene
    )
    adata.obs.loc[mask, f"{gene}_scvi_{label}"] = adata_sub.obs[f"{gene}_scvi"]
    return adata


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


def plot_corr(adata: AnnData, pericyte_mask=None, base: Path="./"):
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
    df_melted.to_csv(base.with_suffix(".tsv"), sep="\t", index=False)
    print(f"Saved melted correlation data to: {base.with_suffix('.tsv')}")

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

    # Define masks
    pericyte_mask, perivascular_mask = define_masks(adata)
    
    # Perivascular model
    adata = subset_anndata(adata, perivascular_mask, "perivascular")
    
    # Pericyte-only model
    adata = subset_anndata(adata, pericyte_mask, "pericytes")

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

    # Plot scatter
    plot_corr(adata, pericyte_mask, airspace_dir / "pericytes_airspace_denoising")

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
