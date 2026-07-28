#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_ct_dx
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
## Two cheap steps: pull the donor -> dataset/lung_condition map out of the
## stroma h5ad obs (h5py, no X read), then fit the per-cell-type three-group
## AGTR1 model on the donor x cell-type table 01.disease_association.R wrote.
## RM-shared controls memory via cores (2000 MB/core); do NOT set --mem.
#SBATCH --time=00:30:00
#SBATCH --output=logs/agtr1_celltype_disease.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** donor -> dataset map from stroma obs ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/04.donor_dataset_map.py \
       --adata stroma.hlca_full.dataset.h5ad \
       --outfile mean_expr/donor_dataset_map.tsv
if [ $? -ne 0 ]; then log_message "Error: donor map failed"; exit 1; fi
conda deactivate

log_message "**** AGTR1 disease effects across stromal cell types ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/05.agtr1_celltype_disease.R \
        --donor-celltype mean_expr/donor_metadata.tsv \
        --donor-map mean_expr/donor_dataset_map.tsv \
        --outdir mean_expr \
        --min-donors 3
if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
