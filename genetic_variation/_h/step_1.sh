#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=prep_gct
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=prep_gct.log
#SBATCH --time=04:00:00

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

## Edit with your job command

echo "**** Run prepare GCT format ****"

python ../_h/01.prepare_gct.py

echo "**** Job ends ****"
date
