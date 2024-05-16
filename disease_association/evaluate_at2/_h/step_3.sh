#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=at2b_prop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --output=proportion.log
#SBATCH --time=01:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module list

echo "**** Review AGTR2 proportion of AT2b cells ****"

python ../_h/03.proportion_at2b.py

echo "**** Job ends ****"
date
