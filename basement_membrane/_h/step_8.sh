#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_count_null
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --array=1-20
#SBATCH --output=logs/agtr1_count_null_%a.log

## Detection-matched null for the COUNT MODEL -- the arbiter for AGTR1 contrasts.
## glmer.nb on 11,680 cells is slow, so the 175 null genes are split 20 ways.
## 17.agtr1_count_null_summarise.R pools the chunks afterwards.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** chunk ${SLURM_ARRAY_TASK_ID} starts ****"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript ../_h/16.agtr1_count_null.R \
        --input ./agtr1_count_input.tsv.gz \
        --null-counts ./null_gene_count_input.tsv.gz \
        --observed ./stats_data/agtr1_count_models.tsv \
        --outdir ./stats_data \
        --chunk ${SLURM_ARRAY_TASK_ID} --n-chunks 20
rc=$?
conda deactivate
log_message "**** chunk ${SLURM_ARRAY_TASK_ID} ends (rc=$rc) ****"
exit $rc
