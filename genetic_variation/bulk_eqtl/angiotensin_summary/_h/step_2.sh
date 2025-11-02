#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=ld_comp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=ld_comparison.log
#SBATCH --time=01:00:00

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

module load plink/2.00a4.6
module list

## Edit with your job command

echo "**** Run LD comparison for AGTR1 top variant and rs5186 ****"

plink2 --bfile ../../../plink_format/_m/genotypes \
       --ld chr3_148474928_A_G_b38 chr3_148742201_A_C_b38 \
       --out ld_analysis

echo "**** Job ends ****"
date
