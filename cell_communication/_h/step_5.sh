#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ccc_cogaps
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
## single (cogaps) receiver scheme; ~ one disease-stratified liana pass.
#SBATCH --time=03:00:00
#SBATCH --output=logs/ccc_cogaps.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

## Orthogonal data-driven receiver scheme: assign pericytes to their dominant
## CoGAPS pattern (nP=5, parsimonious) and add it to the CCC niche.
NP=5
log_message "**** Add CoGAPS receiver scheme (nPatterns=${NP}) ****"
python ../_h/00b.cogaps_receivers.py \
       --niche ccc_niche.h5ad \
       --patterns ../../pericyte_cogaps/_m/patterns_np${NP}.tsv.gz \
       --score-spearman ../../pericyte_cogaps/_m/validation_np${NP}/pattern_score_spearman.tsv.gz \
       --npatterns ${NP} --min-rho 0.15 --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: cogaps receivers failed"; exit 1; fi

log_message "**** Disease-stratified liana into CoGAPS receivers ****"
python ../_h/01.run_liana.py --adata ccc_niche.h5ad --outdir "./" --schemes cogaps
if [ $? -ne 0 ]; then log_message "Error: liana cogaps failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
