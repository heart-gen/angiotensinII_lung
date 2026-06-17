#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ccc_liana
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
## 3 receiver schemes (main + state + receptor) ~ 3x the single-scheme runtime.
#SBATCH --time=08:00:00
#SBATCH --output=logs/ccc_liana.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Disease-stratified liana CCC ****"
python ../_h/01.run_liana.py --adata ccc_niche.h5ad --outdir "./"

if [ $? -ne 0 ]; then log_message "Error: Python failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
