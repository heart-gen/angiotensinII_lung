#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_copd_signfix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/agtr1_copd_signfix.log

## Re-export of the R step ONLY. Used for the 2026-07 contrast sign fix and again
## on 2026-09-07 for the P1-19 cross-cohort sign table.
##
## step_1.sh rebuilds gse136831_ras_pseudobulk.tsv.gz from ipf_dataset.h5ad
## first. The sign fix changes no gene, no model and no input, so that rebuild
## would reproduce the same file. This runs the R step against the existing
## pseudobulk and, unlike the manual rerun that produced the shipped tables
## (TODO P2-32), it leaves a log.
##
## Use step_1.sh, not this, whenever the panel or the source object changes.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

if [ ! -f ./gse136831_ras_pseudobulk.tsv.gz ]; then
    log_message "Error: gse136831_ras_pseudobulk.tsv.gz missing -- run step_1.sh instead"
    exit 1
fi

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
eval "$(conda shell.bash hook)"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Independent AGTR1 COPD/IPF evaluation ****"
Rscript ../_h/01.agtr1_copd_stats.R \
        --pseudobulk ./gse136831_ras_pseudobulk.tsv.gz \
        --outdir ./stats_data \
        --min-cells 5 \
        --hlca-effects ../../_m/mean_expr/agtr1_celltype_disease_effects.tsv
if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
