#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=prep_sce
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --output=prepare.log
#SBATCH --time=08:00:00

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

module load R
module list

## Edit with your job command

Rscript ../_h/00.prepare_sce.R

echo "**** Job ends ****"
date
