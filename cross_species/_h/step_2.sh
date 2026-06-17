#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=mouse_states
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=logs/mouse_states.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Batch-robust mouse state assignment + receptors ****"
python ../_h/02.conserved_states.py --adata mouse_integrated.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: states failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
