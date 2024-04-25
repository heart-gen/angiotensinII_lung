#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=40.0G,h_vmem=40G,h_fsize=50G
#$ -N 'subset_lungmap'
#$ -o ./subset.log
#$ -e ./subset.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
module load gcc/9.1.0
module list

echo "**** Run single cell analysis: lungmap ****"
Rscript ../_h/compartment_subset.R

echo "**** Job ends ****"
date -Is
