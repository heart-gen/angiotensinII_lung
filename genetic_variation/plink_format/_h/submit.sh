#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=subset_plink
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3G
#SBATCH --output=summary.log
#SBATCH --time=03:00:00
#SBATCH --dependency=afterok:4612642

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

echo "**** Run subsetting of genotypes ****"

GENOTYPE_DIR='../../../../inputs/gtex/genotypes/_m'

plink2 --bfile $GENOTYPE_DIR/gtex_v8 \
       --keep ../../_m/keepFam.txt \
       --make-bed --out genotypes

echo "**** Job ends ****"
date
