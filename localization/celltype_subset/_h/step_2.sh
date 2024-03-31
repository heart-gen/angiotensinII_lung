#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=at2_phate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --output=at2_phate.log

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

echo "**** Run PHATE on AT2 ****"

cp -v ../_h/02.at2_phate.ipynb .
quarto render 02.at2_phate.ipynb --execute
rm 02.at2_phate.ipynb

mkdir at2
mv -v at2_* at2/
mv -v 02.at2_phate_files at2/

echo "**** Job ends ****"
date
