#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=GPU-shared
#SBATCH --gpus=v100-32:1
#SBATCH --job-name=mouse_integrate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=5
#SBATCH --time=02:00:00
#SBATCH --output=logs/mouse_integrate.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** scVI integration of mouse lung niche (batch=dataset_id) ****"
python ../_h/01.integrate_mouse.py --adata mouse_lung.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: integration failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
