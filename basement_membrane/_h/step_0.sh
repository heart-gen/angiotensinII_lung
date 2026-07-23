#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_score
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/bm_score.log

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

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Score basement-membrane programs on pericytes ****"
python ../_h/00.bm_score.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad \
       --outdir "./"

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
