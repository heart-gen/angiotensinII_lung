#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_signfix_reexport
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/bm_signfix_reexport.log

## Re-export of the two selectivity/ambient R steps ONLY.
##
## step_2.sh rebuilds bm_pseudobulk_celltype.tsv.gz from the 12 GB ccc_niche
## object first. The contrast sign fix changes no gene, no panel and no model --
## the arms added to bm_panels.py are subsets of TGFB_RESPONSE, so the pseudobulk
## input is byte-identical and rebuilding it would burn an hour to reproduce the
## same file. This runs the two R steps against the existing pseudobulk.
##
## Use step_2.sh, not this, whenever the gene panels or the niche object change.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

if [ ! -f ./bm_pseudobulk_celltype.tsv.gz ]; then
    log_message "Error: bm_pseudobulk_celltype.tsv.gz missing -- run step_2.sh instead"
    exit 1
fi

module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Cross-cell-type basement-membrane selectivity ****"
Rscript ../_h/03.bm_selectivity_stats.R \
        --pseudobulk ./bm_pseudobulk_celltype.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: selectivity step failed"; exit 1; fi

log_message "**** Ambient controls + collagen I stoichiometry ****"
Rscript ../_h/08.fibrillar_ambient.R \
        --pseudobulk ./bm_pseudobulk_celltype.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: ambient-control step failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
