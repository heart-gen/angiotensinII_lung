#!/bin/bash
#SBATCH --partition=GPU-shared
#SBATCH --job-name=cluster_train
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --gpus=v100-32:1
#SBATCH --ntasks-per-node=5
#SBATCH --time=02:00:00
#SBATCH --output=logs/train_model.log

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
module load cuda
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Run analysis ****"
python ../_h/02.train_model.py

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
