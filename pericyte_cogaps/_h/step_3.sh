#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_project
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
## pseudobulk read of the 423k-cell niche + a linear projection; modest mem/cpu.
#SBATCH --time=02:00:00
#SBATCH --output=logs/cogaps_project.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

NP=5

## (1) pseudobulk the niche onto the CoGAPS gene space
log_message "**** Niche pseudobulk (scRNA_env) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/03a.niche_pseudobulk.py \
       --niche ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ../_m/cogaps_genes.tsv \
       --outdir ../_m --min-cells 10
if [ $? -ne 0 ]; then log_message "Error: pseudobulk failed"; exit 1; fi
conda deactivate

## (2) projectR transfer of the pericyte patterns onto the niche pseudobulk.
## NOTE: projectR needs internet to install. If absent, install it ONCE on a
## login node:  R_env/bin/Rscript -e '.libPaths("./.Rlib"); BiocManager::install("projectR", lib="./.Rlib", ask=FALSE, update=FALSE)'
## The R script auto-falls-back to an identical OLS projection if projectR is missing.
log_message "**** projectR transfer (R_env) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/03.project_niche.R \
        --pseudobulk ../_m/niche_pseudobulk_logmean.tsv.gz \
        --samples ../_m/niche_pseudobulk_samples.tsv \
        --npatterns ${NP} \
        --loadings ../_m/feature_loadings_np${NP}.tsv.gz \
        --annotation ../../cell_communication/_m/cogaps_receiver_annotation_np${NP}.tsv \
        --outdir ../_m
if [ $? -ne 0 ]; then log_message "Error: projectR transfer failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
