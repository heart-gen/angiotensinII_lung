#!/bin/bash
#SBATCH --partition=neuron,caracol,gpu
#SBATCH --job-name=run_tensorQTL
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --output=tensorQTL.log
#SBATCH --time=05:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility

module load tensorqtl
module list

## Edit with your job command

echo "**** Run tensorQTL ****"

python ../_h/01.tensorQTL.py

echo "**** Job ends ****"
date
