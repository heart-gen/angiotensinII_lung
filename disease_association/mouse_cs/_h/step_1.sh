#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=mc_plotting
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --output=summary.log
#SBATCH --time=01:00:00

echo "**** Job starts ****"
date -Is

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

echo "**** Plotting qPCR results from mouse ****"

Rscript ../_h/01.plotting_cs.R

echo "**** Job ends ****"
date -Is
