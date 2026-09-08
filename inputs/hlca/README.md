# HLCA version 2

> **Version note.** "Version 2" refers to the **cellxgene HLCA release** (the
> `.h5ad`/`.rds` files below). The Zenodo assets fetched by the same script —
> `HLCA_full_v1.1_emb.h5ad` and `HLCA_reference_model.zip` — are **v1.1**, which
> is the current published embedding/reference-model pair and is *not* a mismatch.
> The two numbering schemes belong to different artifacts.

This directory has the scripts for downloading the
HLCA version 2 data.

This is part of the CZI Single-Cell Biology, 
Human Cell Atlas (HCA).

Data from: [Integrated Human Lung Cell Atlas](https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)

Publication: [Sikkema et al. (2023) Nat Med](https://doi.org/10.1038/s41591-023-02327-2)

> The integrated Human Lung Cell Atlas (HLCA) represents the first large-scale, 
> integrated single-cell reference atlas of the human lung. It consists of over 
> 2 million cells from the respiratory tract of 486 individuals, and includes 49
> different datasets. It is split into the HLCA core, and the extended or full
> HLCA. The HLCA core includes data of healthy lung tissue from 107 individuals,
> and includes manual cell type annotations based on consensus across 6
> independent experts, as well as demographic, biological and technical metadata.
> The datasets in the HLCA core were integrated using scANVI. The HLCA core can 
> be used as a reference to map new datasets onto using scArches. The full HLCA 
> includes 35 further datasets that include donors with various lung diseases. 
>
> These datasets were mapped onto the core with scArches, and include disease 
> annotations as well as cell type annotations transferred from the HLCA core 
> onto the mapped datasets. Note that while the the HLCA includes an integrated,
> batch-correctd low-dimensional embedding, the gene counts themselves were not
> batch-corrected. Both the HLCA core and the full HLCA can be explored below.
> 
> Detailed information about all metadata in the objects, as well as further
> HLCA-related information, can be found on the HLCA landing page: 
> https://github.com/LungCellAtlas/HLCA. Raw counts are available in the 
> downloaded .h5ad at adata.raw.X, in the downloaded .rds at
> seurat_object@assays$RNA@counts, and via CELLxGENE Census.
