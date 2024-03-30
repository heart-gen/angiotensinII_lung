#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=pericyte_phate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=pericyte_phate.log

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

echo "**** Run PHATE on pericytes ****"

cp -v ../_h/01.pericyte_phate.ipynb .
quarto render 01.pericyte_phate.ipynb --execute
rm 01.pericyte_phate.ipynb

mkdir pericytes
mv -v pericyte_* pericytes/
mv 01.pericyte_phate_files pericytes/

echo "**** Job ends ****"
date
