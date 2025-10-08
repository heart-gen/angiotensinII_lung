#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=subset_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00
#SBATCH --output=subset_data.log

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

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Run analysis ****"
python ../_h/01.subset_data.py

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
