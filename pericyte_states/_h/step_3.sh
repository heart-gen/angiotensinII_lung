#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_agtr1_lenses
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/agtr1_lenses.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** AGTR1 by program: raw / detection / scVI-denoised lenses ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/03.agtr1_lenses.R \
        --meta pericytes_states_metadata.tsv.gz \
        --denoise ../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv \
        --den-model Pericyte-only-trained \
        --outdir stats_data
if [ $? -ne 0 ]; then log_message "Error: R step failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
