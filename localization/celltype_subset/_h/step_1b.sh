#!/bin/bash
#SBATCH --partition=EM
#SBATCH --job-name=full_subset
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=24
#SBATCH --time=03:00:00
#SBATCH --output=full_subsetting.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
export BASILISK_EXTERNAL_DIR=/ocean/projects/bio250020p/shared/opt/basilisk_cache

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Run subsetting ****"
MODEL="full"

Rscript ../_h/01.subset_ct.py --model "$MODEL"

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
