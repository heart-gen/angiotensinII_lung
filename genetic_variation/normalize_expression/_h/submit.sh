#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=norm_expr
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --output=summary.log
#SBATCH --time=04:00:00
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

module load htslib
module list

## Edit with your job command

echo "**** Run normalize expression ****"

GTF="/dcs05/lieber/hanlab/jbenjami/resources/gtex/gencode.v26.GRCh38.genes.gtf"

ln -sfn ../_h/rnaseqnorm.py .

python ../_h/eqtl_prepare_expression.py \
       --feature gene -o . ../../_m/tpm.gct ../../_m/counts.gct $GTF \
       ../../_m/sample_id_to_subject_id.tsv ../../_m/vcf_chr_list.txt genes

echo "**** Job ends ****"
date
