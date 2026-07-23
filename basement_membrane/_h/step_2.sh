#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_selectivity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/bm_selectivity.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Donor x cell-type pseudobulk over the CCC niche ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

python ../_h/02.niche_pseudobulk.py \
       --adata ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ./bm_panel_genes.tsv \
       --outfile ./bm_pseudobulk_celltype.tsv.gz

if [ $? -ne 0 ]; then
    log_message "Error: pseudobulk build failed"
    exit 1
fi
conda deactivate

log_message "**** Cross-cell-type basement-membrane selectivity ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript ../_h/03.bm_selectivity_stats.R \
        --pseudobulk ./bm_pseudobulk_celltype.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --outdir ./stats_data \
        --min-cells 5

if [ $? -ne 0 ]; then
    log_message "Error: R execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
