#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=core_subset
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=40
#SBATCH --cpus-per-task=1
#SBATCH --output=core_subsetting.log

echo "**** Job starts ****"
date

echo "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module load R
module list

echo "**** Run subsetting ****"
export BASILISK_EXTERNAL_DIR=/ocean/projects/bio250020p/shared/opt/basilisk_cache

model="core"
Rscript ../_h/00.subset_ct.R "$model"

echo "**** Job ends ****"
date
