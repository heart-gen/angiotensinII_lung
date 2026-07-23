#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=mouse_comparability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=logs/mouse_comparability.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Cross-species comparability audit + cell-level Agtr1a ****"
python ../_h/04.species_comparability.py \
    --adata mouse_integrated.h5ad \
    --metadata mouse_states_metadata.tsv.gz \
    --outdir stats_data
if [ $? -ne 0 ]; then log_message "Error: comparability audit failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
