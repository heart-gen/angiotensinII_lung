"""
Single-Cell Data Clustering and Visualization -- Training.
"""
import torch, scvi
import session_info
import scanpy as sc

torch.set_float32_matmul_precision("high")

def train_model(adata, patience: int = 10, max_epochs: int = 500):
    # Setup and train model SCANVI model
    scvi.model.SCVI.setup_anndata(
        adata, layer="counts", labels_key="subcluster",
    )
    vae = scvi.model.SCVI(adata, n_latent=30, n_layers=2)
    vae.train(
        max_epochs=max_epochs, early_stopping=True,
        early_stopping_patience=patience * 2, validation_size=0.1,
        check_val_every_n_epoch=max(1, patience // 2),
        precision="16-mixed" if torch.cuda.is_available() else "32"
    )

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        vae, unlabeled_category="unknown", labels_key="subcluster"
    )
    scanvi_model.train(
        max_epochs=max_epochs, early_stopping=True,
        early_stopping_patience=patience,
        early_stopping_monitor="elbo_validation",
        check_val_every_n_epoch=max(1, patience // 2),
        plan_kwargs={"weight_decay": 0.0}, validation_size=0.1,
        precision="16-mixed" if torch.cuda.is_available() else "32",
    )

    # Save model for downstream
    scanvi_model.save("scanvi_model/", overwrite=True)


def main():
    # Load data
    ref_hvg = sc.read_h5ad("ref_hvg.h5ad")    

    # Train model
    train_model(ref_hvg)

    # Session information
    session_info.show()


if __name__ == "__main__":
    main()
