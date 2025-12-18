import pandas as pd
import scanpy as sc
import session_info
import seaborn as sns
import logging, argparse
from pathlib import Path
from scipy.stats import spearmanr

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
                        default=Path("./results/pericytes_with_airspace_score.h5ad"),
                        help="HLCA airspace scored AnnData")
    parser.add_argument("--imputed-agrt1", type=Path,
                        default=Path("./results/airspace/pericytes_airspace_denoising.tsv"),
                        help="Location of denoised AGTR1 data")
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    return parser.parse_args()


def save_figure(fig, base_path: Path):
    fig.savefig(base_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def clean_data(adata, denoised_df, pericyte_label="Pericytes", key="subclusters"):
    meta_df = adata.obs.copy()
    meta_df = meta_df[meta_df[key] == pericyte_label].copy()

    meta_df = meta_df.dropna(subset=[
        "airspace_score", "AGTR1_detect", "donor_id", 
        "sex", "age_or_mean_of_age_range",
    ]).copy()

    meta_df["AGTR1_detect"] = meta_df["AGTR1_detect"].astype(int)
    meta_df["sex"] = meta_df["sex"].astype("category")
    meta_df["age"] = meta_df["age_or_mean_of_age_range"].astype(int)

    df = pd.merg(meta_df, denoised_df, on=[], how="inner")
    return None

    
def main():
    # Load parser and logging
    args = parse_args()
    configure_logging()

    # Set parameters
    outdir = args.outdir
    outdir.mkdir(exist_ok=True, parents=True)

    # Load AnnData
    adata = sc.read_h5ad(args.adata)
    denoised = pd.read_csv(args.imputed_agrt1, sep="\t", index_col=0)

    # Extract, clean, and merge metadata
    keep_cols = []


if __name__ == "__main__":
    main()
