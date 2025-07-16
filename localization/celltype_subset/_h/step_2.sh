#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=subclustering
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:06:00
#SBATCH --output=visualize.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

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
conda activate /ocean/projects/bio250020p/shared/opt/env/AI_env

log_message "**** Run subclustering ****"

python ../_h/02.sub-clustering.py

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
