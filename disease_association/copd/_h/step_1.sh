#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=lung_copd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15gb
#SBATCH --output=summary.log
#SBATCH --time=06:00:00

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

echo "**** Run single cell profile: COPD ****"

Rscript ../_h/01.angiotensinII_analysis.R

echo "**** Job ends ****"
date
