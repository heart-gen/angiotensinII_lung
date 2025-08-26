"""
Trajectory analysis with PHATE/Diffusion Map pseudotime
- Reproducible seeding
- Auto-compute PHATE if missing
- Model-aware I/O (core/full)
- Side-by-side pseudotime plots (PHATE & Diffusion Map)
- Optional smoothed gene dynamics along pseudotime
- Works with stromal and pericytes subclusters
"""
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import session_info
from os import makedirs, path
import matplotlib.pyplot as plt
from scipy.sparse import issparse

def set_seed(seed=13):
    np.random.seed(seed)
    try:
        import random
        random.seed(seed)
    except Exception:
        pass
    sc.settings.set_figure_params(dpi=300)
    sc.settings.verbosity = 2


def ensure_dir(dirname):
    makedirs(dirname, exist_ok=True)


def load_adata(filepath):
    """Load AnnData object from .h5ad file"""
    return sc.read_h5ad(filepath)


def ensure_phate(adata, use_rep="X_pca_harmony", knn=15, decay=20, seed=13):
    """
    Ensure adata.obsm['X_phate'] exists; compute if missing.
    Requires 'phate' package and a representation to neighbor on.
    """
    if "X_phate" in adata.obsm:
        return adata

    if use_rep in adata.obsm:
        rep = use_rep
    elif "X_pca" in adata.obsm:
        rep = "X_pca"
    elif "PCA" in adata.obsm:
        adata.obsm["X_pca"] = adata.obsm["PCA"]
        rep = "X_pca"
    else:
        rep = None # Fallback

    # Compute neighbors on chosen rep (for PHATE graph)
    sc.pp.neighbors(adata, n_neighbors=50, use_rep=rep, random_state=seed)

    try:
        import phate
        ph = phate.PHATE(knn=knn, decay=decay, n_jobs=-1, random_state=seed)
        X = adata.obsm[rep] if rep is not None else adata.X
        adata.obsm["X_phate"] = ph.fit_transform(X)
    except ImportError:
        raise RuntimeError("phate is not installed but required to compute PHATE embedding.")
    return adata


def compute_pseudotime(adata, root_cell=None, phate_key='X_phate',
                       groups_key="leiden", seed=13):
    """Compute diffusion pseudotime; neighbors on PHATE, then Diffusion Map + DPT."""
    # Check embedding
    if phate_key not in adata.obsm:
        raise KeyError(f"Embedding '{phate_key}' not found. Compute PHATE before pseudotime.")

    # Build neighborhood on PHATE
    sc.pp.neighbors(adata, use_rep=phate_key, random_state=seed)
    
    # Compute diffusion Map
    sc.tl.diffmap(adata)
    
    # PAGA graph (optional but useful for topology)
    if groups_key in adata.obs:
        sc.tl.paga(adata, groups=groups_key)

    # Root cell: lowest PHATE-1 by default
    if root_cell is None:
        root_cell = int(np.argmin(adata.obsm[phate_key][:, 0]))
    adata.uns['iroot'] = root_cell
    
    # Diffusion pseudotime
    sc.tl.dpt(adata)
    return adata


def _save_plot(plot_func, outdir, model, fname, **kwargs):
    formats = ["png", "pdf"]
    plt.figure()
    ax = plot_func(show=False, **kwargs)
    fig = ax.get_figure() if hasattr(ax, "get_figure") else plt.gcf()
    for ext in formats:
        plt.savefig(path.join(outdir, model, f"{fname}.{ext}"),
                    dpi=300, bbox_inches='tight')
    plt.close(fig)


def plot_pseudotime(adata, phate_key="X_phate", outdir="figures", model="core",
                    save_prefix="pseudotime", cmap="viridis"):
    """
    Plot pseudotime on PHATE and Diffusion Map embeddings (side by side).
    Saves as PNG and PDF.
    """
    ensure_dir(path.join(outdir, model))
    # PHATE
    _save_plot(lambda **kwargs: sc.pl.embedding(
        adata, basis=phate_key, color="dpt_pseudotime",
        title="Pseudotime (PHATE)", cmap=cmap, **kwargs),
               outdir, model, f"{save_prefix}_phate")

    # Diffusion map
    _save_plot(lambda **kwargs: sc.pl.embedding(
        adata, basis="X_diffmap", color="dpt_pseudotime",
        title="Pseudotime (Diffusion Map)", cmap=cmap, **kwargs),
               outdir, model, f"{save_prefix}_diffmap")


def _get_layer_or_X(adata, gene, prefer_layer="logcounts"):
    """Return 1D expression vector for a gene from layer if present, else X."""
    if prefer_layer in adata.layers:
        X = adata[:, gene].layers[prefer_layer]
    else:
        X = adata[:, gene].X
    return X.toarray().ravel() if issparse(X) else np.ravel(X)


def plot_gene_dynamics(adata, genes, outdir="figures", model="core",
                       save_prefix="gene_dynamics", window_frac=0.02,
                       prefer_layer="logcounts"):
    """
    Smoothed gene expression (rolling mean) along pseudotime.
    - window size = max(5, floor(window_frac * n_cells))
    """
    ensure_dir(path.join(outdir, model))
    if "dpt_pseudotime" not in adata.obs:
        raise KeyError("dpt_pseudotime not found; run compute_pseudotime first.")
    pseudotime = adata.obs["dpt_pseudotime"].values
    order_pt = np.argsort(pseudotime)
    n = len(order_pt)
    win = max(5, int(np.floor(window_frac * n)))

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for gene in genes:
        if gene not in adata.var["feature_name"].values:
            print(f"Gene {gene} not found; skipping.")
            continue
        g_annot = pd.DataFrame(adata.var)
        new_gene = g_annot[(g_annot["feature_name"] == gene)].index[0]
        vals = _get_layer_or_X(adata, new_gene, prefer_layer=prefer_layer)[order_pt]
        # simple rolling mean
        kernel = np.ones(win) / win
        smooth = np.convolve(vals, kernel, mode="same")
        ax.plot(pseudotime[order_pt], smooth, label=gene, linewidth=2)
        
    ax.set_xlabel("Pseudotime (DPT)")
    ax.set_ylabel("Smoothed Expression")
    ax.legend(frameon=False)
    plt.savefig(path.join(outdir, model, f"{save_prefix}.png"), dpi=300,
                bbox_inches='tight')
    plt.savefig(path.join(outdir, model, f"{save_prefix}.pdf"),
                bbox_inches='tight')
    plt.close(fig)


def plot_paga(adata, outdir="figures", model="core", save_prefix="paga"):
    """Optional: PAGA connectivity & positions (requires sc.tl.paga)."""
    ensure_dir(path.join(outdir, model))
    try:
        _save_plot(lambda **kwargs: sc.pl.paga(adata, **kwargs),
                   outdir, model, f"{save_prefix}")
    except Exception as e:
        print(f"PAGA plot skipped: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Trajectory analysis with model-aware outputs (core/full)")
    parser.add_argument("--model", type=str, default="core",
                        choices=["core", "full"],
                        help="Model type: 'core' or 'full'. Default: core")
    parser.add_argument("--stroma", action="store_true",
                        help="Flag to indicate input data is stroma (default: False)")
    parser.add_argument("--outdir", type=str, default="figures",
                        help="Output directory for figures")
    parser.add_argument("--genes", type=str, nargs="*", default=["ACTA2", "AGTR1"],
                        help="Marker genes for dynamics plot")
    args = parser.parse_args()

    set_seed(13)

    # Load data depending on stroma flag
    if args.stroma:
        cluster_key = "subclusters"; knn=15
        input_file = f"stroma.hlca_{args.model}.clustered.h5ad"
        output_file = f"stroma.hlca_{args.model}.clustered.analysis.h5ad"
    else:
        cluster_key = "leiden"; knn=20
        input_file = f'pericyte.hlca_{args.model}.subclustered.h5ad'
        output_file = f'pericyte.hlca_{args.model}.subclustered.analysis.h5ad'

    # Load data
    adata = load_adata(input_file)

    # Ensure PHATE exists (compute if needed)
    adata = ensure_phate(adata, use_rep="X_pca_harmony",
                         knn=knn, decay=15, seed=13)
    
    # Trajectory inference & pseudotime
    adata = compute_pseudotime(adata, phate_key="X_phate",
                               groups_key=cluster_key, seed=13)

    # Plots
    plot_pseudotime(adata, phate_key="X_phate",
                    outdir=args.outdir, model=args.model)
    plot_gene_dynamics(adata, args.genes,
                       outdir=args.outdir, model=args.model)
    plot_paga(adata, outdir=args.outdir, model=args.model)
    
    # Save processed data
    adata.write(output_file)
    
    # Session information
    session_info.show()


if __name__ == '__main__':
    main()
    
