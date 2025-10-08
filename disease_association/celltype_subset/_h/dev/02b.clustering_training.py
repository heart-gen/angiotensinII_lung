"""
Single-Cell Data Clustering and Visualization

This script clusters the data and extracts marker genes for cell type annotation.

Author: Kynon J Benjamin
Date: August 23, 2024
Revision: September 27, 2025
"""
import torch, scvi
import session_info
import scanpy as sc

torch.set_float32_matmul_precision("high")

def train_model():
    # Load data
    ref_hvg = sc.read_h5ad("ref_hvg.h5ad")
    
    # Setup and train model SCANVI model
    scvi.model.SCVI.setup_anndata(ref_hvg, layer="counts")
    vae = scvi.model.SCVI(ref_hvg, n_latent=30, n_layers=2)
    vae.train()

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="unknown",
        labels_key="celltype"
    )
    scanvi_model.train(
        max_epochs=500,
        early_stopping=True,
        early_stopping_patience=10,
        early_stopping_monitor="elbo_validation",
        plan_kwargs={"weight_decay": 0.0},
        check_val_every_n_epoch=5
    )

    # Save model for downstream
    scanvi_model.save("scanvi_model/", overwrite=True)


def main():
    # Train model
    train_model()

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
