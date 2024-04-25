#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=50G
#$ -N lungmap_dataset
#$ -o ./prepare.log
#$ -e ./prepare.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

## Edit with your job command

echo "**** Run COPD conversion ****"
Rscript ../_h/generate_lungmap_data.R

echo "**** Job ends ****"
date
