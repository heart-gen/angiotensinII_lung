#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=extract_angiotensin
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --output=summary.log
#SBATCH --time=03:00:00

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
ANGIOTENSIN="../../gtex_finemapping/_h/angiotensin.out"
CONDITIONAL="../../_m/GTEx_v8.conditional.txt.gz"
SIGNIF="../../_m/GTEx_v8.signif_variants.txt.gz"
SUSIE="../../_m/GTEx_v8.susie.txt.gz"

echo "**** Conditional eQTL ****"

cat <(zcat $CONDITIONAL | head -1 | cut -f2,8,12-14,16) \
    <(zcat $CONDITIONAL | grep -f $ANGIOTENSIN | cut -f2,8,12-14,16) \
    > GTEx_v8.conditional.angiotensin.txt

echo "**** Signficant eQTL ****"

cat <(zcat $SIGNIF | head -1) \
    <(zcat $SIGNIF | grep -f $ANGIOTENSIN) \
    > GTEx_v8.signif_variants.angiotensin.txt

echo "**** Fine-mapped eQTL -- SuSIE ****"

cat <(zcat $SUSIE | head -1) \
    <(zcat $SUSIE | grep -f $ANGIOTENSIN) \
    > GTEx_v8.susie.angiotensin.txt

echo "**** Job ends ****"
date
