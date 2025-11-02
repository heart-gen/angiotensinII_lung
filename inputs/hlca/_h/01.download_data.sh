#!/bin/bash

echo "Downloading data"
wget https://github.com/LungCellAtlas/HLCA/blob/main/docs/HLCA_metadata_explanation.csv
curl -o hlca_full.h5ad https://datasets.cellxgene.cziscience.com/3ab47484-a3eb-4f6a-beea-670e1a8fc1e8.h5ad
curl -o hlca_core.h5ad https://datasets.cellxgene.cziscience.com/7a3f08f9-5d07-4ddd-a8fe-5967dd34f35f.h5ad
curl -o hlca_core.rds https://datasets.cellxgene.cziscience.com/2aa90e63-9a6d-444d-8343-8fc2a9921797.rds

echo "Downloading models"
wget https://zenodo.org/records/7599104/files/HLCA_full_v1.1_emb.h5ad
wget https://zenodo.org/records/7599104/files/HLCA_reference_model.zip
wget https://zenodo.org/records/7599104/files/HLCA_reference_model_gene_order_ids_and_symbols.csv
wget https://zenodo.org/records/7599104/files/HLCA_surgery_models.zip
