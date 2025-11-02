import scanpy as sc
import matplotlib.pyplot as plt

def plot_umap_scanvi(adata, out_prefix="umap_X_scANVI", color_keys=None, dpi=300):
    """
    Make UMAP using the SCANVI latent (X_scANVI) and save PNG/PDF.

    Parameters
    ----------
    adata : AnnData
        Must have adata.obsm["X_scANVI"].
    out_prefix : str
        Path/prefix for output files (no extension).
    color_keys : list[str] or None
        obs columns to color by. Defaults to ["predicted_labels", "donor"].
    dpi : int
        Resolution for PNG.
    """
    if color_keys is None:
        color_keys = ["predicted_labels", "donor"]

    # build neighbors/UMAP on X_scANVI
    sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=15)
    sc.tl.umap(adata)

    # 1 row per color key
    n = len(color_keys)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))

    if n == 1:
        axes = [axes]  # force list so we can iterate

    for ax, key in zip(axes, color_keys):
        sc.pl.umap(
            adata, color=key, ax=ax, show=False, title=key,
            legend_loc="on data" if key == "predicted_labels" else "right margin",
        )

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=dpi, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    adata = sc.read_h5ad("lungmap_transferred.h5ad")
    plot_umap_scanvi(adata, out_prefix="qc_plots/umap_scanvi")


if __name__ == "__main__":
    main()
