#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=run_prism
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --array=1-3%1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3gb
#SBATCH --output=summary.prism.%A_%a.log
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

echo "**** Run BayesPrism ****"

Rscript ../_h/02.run_prism.R --select $SLURM_ARRAY_TASK_ID

echo "**** Job ends ****"
date
