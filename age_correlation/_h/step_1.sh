#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=ontogeny_lungmap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --output=ontogeny_lungmap.log

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

echo "**** Run single cell analysis: lungmap ****"

Rscript ../_h/01.ontogeny_analysis.R

echo "**** Job ends ****"
date -Is
