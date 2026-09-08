#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_copd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=logs/agtr1_copd.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

## The pseudobulk builder is REUSED from basement_membrane/_h/05.bm_copd.py
## rather than copied. It is fully parameterised by --genes and already encodes
## the GSE136831 compartment map (Pericyte and SMC kept separate AND pooled as
## Mural) plus the counts-correct pseudobulk (sum counts -> CP10K -> log1p).
## Copying it would have created a second thing to keep in sync with the
## compartment definitions the COPD basement-membrane result already depends on.
log_message "**** Donor x compartment RAS-panel pseudobulk (GSE136831) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../../../basement_membrane/_h/05.bm_copd.py \
       --adata ../../ipf_analysis/_m/ipf_dataset.h5ad \
       --genes ../_h/agtr1_panel_genes.tsv \
       --demo ../../ipf_analysis/_h/sample_demo.csv \
       --outfile ./gse136831_ras_pseudobulk.tsv.gz
if [ $? -ne 0 ]; then log_message "Error: pseudobulk failed"; exit 1; fi
conda deactivate

log_message "**** Independent AGTR1 COPD/IPF evaluation ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/01.agtr1_copd_stats.R \
        --pseudobulk ./gse136831_ras_pseudobulk.tsv.gz \
        --outdir ./stats_data \
        --min-cells 5 \
        --hlca-effects ../../_m/mean_expr/agtr1_celltype_disease_effects.tsv
if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi

## Cell-floor sensitivity (P2-28, added 2026-09-07). These two refits are quoted
## in writings/AGTR1_DISEASE_DIRECTION.md but were previously run by hand with
## `--outdir /tmp/mc3`, so nothing on disk reproduced them. The floor matters
## here: at >=5 cells the Control arm collapses to one donor in some
## compartments, and the >=3 / >=2 rungs are what show whether the direction is
## a floor artifact.
for MC in 3 2; do
    log_message "**** Cell-floor sensitivity: >= ${MC} cells ****"
    Rscript ../_h/01.agtr1_copd_stats.R \
            --pseudobulk ./gse136831_ras_pseudobulk.tsv.gz \
            --outdir ./stats_data/sensitivity_floor${MC} \
            --min-cells ${MC} \
            --hlca-effects ../../_m/mean_expr/agtr1_celltype_disease_effects.tsv
    if [ $? -ne 0 ]; then log_message "Error: floor-${MC} refit failed"; exit 1; fi
done
conda deactivate

## Old-vs-new unit decomposition (P2-28). Explains the fibroblast/myofibroblast
## sign flip between ipf_analysis/_h/01 (per-cell mean, no floor) and this module
## (donor x compartment pseudobulk). Was in no step script at all.
log_message "**** Unit decomposition: per-cell vs pseudobulk ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/02.unit_decomposition.py --outfile ./stats_data/unit_decomp.tsv
if [ $? -ne 0 ]; then log_message "Error: unit decomposition failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
