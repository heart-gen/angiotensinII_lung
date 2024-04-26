#!/bin/bash
#SBATCH --partition=caracol,gpu
#SBATCH --job-name=tensorqtl_post
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20gb
#SBATCH --output=postprocessing.log
#SBATCH --time=04:00:00

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

echo "**** Run tensorQTL post processing ****"

PREFIX="GTEx_v8"
python ../_h/03.post_processing.py --prefix $PREFIX

echo "**** Job ends ****"
date
