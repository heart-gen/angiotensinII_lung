# This module moved

The disease-association analyses left this repository on 2026-09-10. They now
live in **heart-gen/lung-pericyte-analysis**, together with `sensitivity/`,
`niche_index/_h/01.niche_disease_stats.R`, the disease figures and Tables
S13A/S13E/S13F/S14. This repository keeps the localization, pericyte-state,
matrix and signalling work and makes no disease claims.

Start at that repository's `disease_association/README.md`.

## Why this directory still exists on disk

It holds large binaries that were deliberately **not** copied, because they are
re-derivable and moving 27 GB to make a second copy helps nobody:

| File | Size | Rebuild with |
| ---- | ---- | ------------ |
| `_m/stroma.hlca_full.dataset.h5ad` | 5.0 GB | `disease_association/_h/step_0.sh` over there |
| `pericyte_analysis/_m/{query,ref}_hvg.h5ad`, `results/clustered_data.h5ad` | 59 MB | that module is retired; see its `RETIRED.md` |
| `pericyte_analysis/_m/scanvi_model/model.pt` | | as above |

Nothing here is tracked by git and nothing here is read by any script in this
repository. Delete the directory once the disease repository has rebuilt or
copied what it needs.

## Two files did not move -- they were relocated within this repository

Both were shared builders that happened to live under `disease_association/`
despite not being disease analyses:

| Was | Is now | Read by |
| --- | ------ | ------- |
| `_h/00.preprocess_reference.py` | `inputs/hlca/_h/02.preprocess_reference.py` | `cell_communication/_h/step_0.sh` |
| `ipf_analysis/_h/02.generate_ipf_data.R` | `inputs/ipf/_h/01.generate_ipf_data.R` | `basement_membrane/_h/step_3.sh` |

Their outputs moved with them, to `inputs/hlca/_m/hlca_full.dataset.h5ad` and
`inputs/ipf/_m/ipf_dataset.h5ad`. `ipf_analysis/_h/sample_demo.csv` was copied to
`inputs/ipf/_h/sample_demo.csv` for `tables/_h/01.cohort_mouse.R`. The disease
repository carries its own copy of each.
