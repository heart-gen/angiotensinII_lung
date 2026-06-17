#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-small
#SBATCH --job-name=pathway_balance
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=logs/pathway_balance.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** AT1R/AT2R pathway-balance scoring ****"
python ../_h/00.pathway_balance.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: Python failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
