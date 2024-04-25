#!/bin/bash
#$ -cwd
#$ -l gpu,mem_free=80G,h_vmem=80G,h_fsize=100G
#$ -N de_markers
#$ -o ./markers.log
#$ -e ./markers.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

echo "**** Run tensorQTL ****"
export CUDA_VISIBLE_DEVICES=3

python3 ../_h/marker_genes.py

echo "**** Job ends ****"
date
