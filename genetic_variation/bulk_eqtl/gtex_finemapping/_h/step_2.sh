#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=subset_angiotensin
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --output=angiotensin.log
#SBATCH --time=02:00:00
#SBATCH --dependency=afterok:4617930

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

ANGIOTENSIN="../_h/angiotensin.out"

echo "**** CaVEMaN ****"
CAVEMAN="GTEx_v8_finemapping_CaVEMaN/GTEx_v8_finemapping_CaVEMaN.txt.gz"
cat <(zcat $CAVEMAN | head -1) <(zcat $CAVEMAN | grep -f $ANGIOTENSIN | grep -i lung) \
    > CaVEMaN_angiotensin_lung.txt

echo "**** CAVIAR ****"
CAVIAR="GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt.gz"
cat <(zcat $CAVIAR | head -1) <(zcat $CAVIAR | grep -f $ANGIOTENSIN | grep -i lung) \
    > CAVIAR_angiotensin_lung.txt

echo "**** DAPG ****"
DAPG="GTEx_v8_finemapping_DAPG"
cat <(zcat $DAPG/GTEx_v8_finemapping_DAPG.txt.gz | head -1) \
    <(zcat $DAPG/GTEx_v8_finemapping_DAPG.txt.gz | grep -f $ANGIOTENSIN | grep -i lung) \
    > DAPG_angiotensin_lung.txt
zcat $DAPG/GTEx_v8_finemapping_DAPG.vcf.gz | grep -f $ANGIOTENSIN | \
    grep -i lung > DAPG_angiotensin_lung.vcf
zcat $DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz | grep -f $ANGIOTENSIN | \
    grep -i lung > DAPG_angiotensin_lung.cs95.txt

echo "**** Copy readme files ****"
cp -v GTEx_v8_finemapping_CaVEMaN/README.txt README_CaVEMaN.txt
cp -v GTEx_v8_finemapping_CAVIAR/README.txt README_CAVIAR.txt
cp -v GTEx_v8_finemapping_DAPG/README.md README_DAPG.md

echo "**** Remove old files ****"
rm GTEx_v8_finemapping_*/*
rmdir GTEx_v8_finemapping_*/

echo "**** Job ends ****"
date
