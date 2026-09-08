"""
Single-Cell Data Clustering and Visualization -- Training.
"""
import pandas as pd
import torch, scvi
import session_info
import scanpy as sc

torch.set_float32_matmul_precision("high")

# REPRODUCIBILITY (P2-23, added 2026-09-07). Only the numpy and scanpy seeds were
# set; scvi-tools draws from its own generator, so no fit here was reproducible.
scvi.settings.seed = 13

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

    # TRAINING CURVES (P2-23, added 2026-09-07). Training stopped at epoch 11 of
    # 500 -- "elbo_validation did not improve in the last 10 records" -- and with
    # check_val_every_n_epoch = 5 and patience = 10 that criterion can fire after
    # very few validation evaluations. No curve or metric table was written, so
    # whether the model converged usefully could not be assessed from any
    # artifact. It can now.
    for name, model in (("scvi", vae), ("scanvi", scanvi_model)):
        hist = getattr(model, "history", None)
        if not hist:
            continue
        df = pd.concat(hist.values(), axis=1)
        df.index.name = "epoch"
        df.to_csv(f"training_history_{name}.tsv", sep="\t")
        print(f"wrote training_history_{name}.tsv ({len(df)} epochs recorded)")

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
