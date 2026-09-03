#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=tgfb_specificity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=06:00:00
#SBATCH --output=logs/tgfb_specificity.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

## Is the TGF-beta/BM association TGF-beta signalling, or a generic
## immediate-early / dissociation-stress program? Design and decision rule are
## pre-specified in ../_h/TGFB_SPECIFICITY_PLAN.md.

log_message "**** Detection-matched null panels + leave-one-gene-out ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

python ../_h/11.tgfb_null_panels.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad \
       --outdir "./" \
       --n-null 1000

if [ $? -ne 0 ]; then log_message "Error: null panel generation failed"; exit 1; fi
conda deactivate

log_message "**** Specificity models against the matched null ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript ../_h/12.tgfb_specificity.R \
        --bm-meta ./bm_metadata.tsv.gz \
        --state-meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
        --null-pseudobulk ./tgfb_null_pseudobulk.tsv.gz \
        --outdir ./stats_data \
        --min-cells 5

if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
