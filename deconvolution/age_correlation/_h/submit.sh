#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=age_correlation
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15gb
#SBATCH --output=cell_proportion.log
#SBATCH --time=02:30:00

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

echo "**** Run correlation -- cell proportion ****"

Rscript ../_h/01.gtex_comparison.R

echo "**** Job ends ****"
date
