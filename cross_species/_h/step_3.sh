#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=mouse_conserv
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/mouse_conserv.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** RETIRED -- 03.conservation_stats.R is no longer run ****"
log_message "Its state-level contrasts are 96% smooth muscle and its continuous"
log_message "test runs on the dense scvi_corrected layer; see the header of"
log_message "../_h/03.conservation_stats.R. Run step_4.sh instead."
conda deactivate
log_message "**** Job ends ****"
