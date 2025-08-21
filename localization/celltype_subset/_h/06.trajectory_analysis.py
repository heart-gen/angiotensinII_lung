## This script performs trajectory analysis
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from os import makedirs, path
import matplotlib.pyplot as plt

def load_adata(filepath):
    '''Load AnnData object from .h5ad file'''
    return sc.read_h5ad(filepath)


def compute_pseudotime(adata, root_cell=None, phate_key='X_phate',
                       groups_key="leiden"):
    '''Compute diffusion pseudotime using PHATE embedding and set root'''
    # Build neighborhood based on PHATE embedding
    sc.pp.neighbors(adata, use_rep=phate_key, random_state=13)
    # Compute diffusion map for pseudotime
    sc.tl.diffmap(adata)    
    # PAGA for trajectory inference
    sc.tl.paga(adata, groups=groups_key)    
    # Root cell: choose lowest PHATE-1 by default
    if root_cell is None:
        root_cell = np.argmin(adata.obsm[phate_key][:, 0])
    adata.uns['iroot'] = root_cell
    # Diffusion pseudotime
    sc.tl.dpt(adata)
    return adata


def plot_pseudotime(adata, phate_key="X_phate", outdir="plots", model="core",
                    save_prefix="pseudotime_scatter"):
    """
    Plot pseudotime on PHATE and Diffusion Map embeddings (side by side).
    Saves as PNG and PDF.
    """
    makedirs(f"{outdir}/{model}", exist_ok=True)
    sc.pl.embedding(
        adata,
        basis=phate_key,
        color="dpt_pseudotime",
        save=None,
        show=False
    )
    plt.savefig(path.join(outdir, model, f"{save_prefix}.png"),
                dpi=300, bbox_inches="tight")
    plt.savefig(path.join(outdir, model, f"{save_prefix}.pdf"),
                bbox_inches="tight")
    plt.close()


def plot_gene_dynamics(adata, genes, outdir="plots", model="core",
                       save_prefix="gene_dynamics"):
    """
    Plot smoothed gene expression along pseudotime for selected marker genes.
    Uses a rolling mean to smooth values.
    """
    makedirs(f"{outdir}/{model}", exist_ok=True)
    pseudotime = adata.obs['dpt_pseudotime']
    df = pd.DataFrame(index=adata.obs_names)
    df['pseudotime'] = pseudotime

    for gene in genes:
        if gene in adata.var_names:
            df[gene] = adata[:, gene].X.toarray().flatten()
        else:
            print(f"Gene {gene} not found, skipping.")

    df = df.sort_values("pseudotime")

    fig, ax = plt.subplots(figsize=(8, 5))
    for gene in genes:
        if gene in df:
            ax.plot(df['pseudotime'], df[gene].rolling(50, min_periods=1).mean(),
                    label=gene, linewidth=2)

    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Smoothed Expression")
    ax.legend()
    plt.tight_layout()
    plt.savefig(path.join(outdir, model, f"{save_prefix}.png"),
                dpi=300, bbox_inches="tight")
    plt.savefig(path.join(outdir, model, f"{save_prefix}.pdf"),
                bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Pericyte/stroma subclustering with model-aware outputs")
    parser.add_argument("--model", type=str, default="core",
                        help="Model type: 'core' or 'full'. Default: core")
    args = parser.parse_args()
    model = args.model

    # Load data
    adata = load_adata(f'pericyte.hlca_{model}.subclustered.h5ad')
    # Trajectory inference & pseudotime
    adata = compute_pseudotime(adata)
    # Save processed data
    adata.write(f'pericyte.hlca_{model}.subclustered.analysis.h5ad')
    # Session information
    session_info.show()


if __name__ == '__main__':
    main()
    
