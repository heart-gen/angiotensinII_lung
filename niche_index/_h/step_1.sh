#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-small
#SBATCH --job-name=niche_stats
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/niche_stats.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

module purge
module load anaconda3/2024.10-1
module list
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Niche-index disease statistics ****"
Rscript ../_h/01.niche_disease_stats.R

if [ $? -ne 0 ]; then log_message "Error: Rscript failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
