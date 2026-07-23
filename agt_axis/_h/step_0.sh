#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agt_ras_landscape
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/agt_ras_landscape.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

# Reuses the pseudobulk builder from basement_membrane with a different gene
# list, rather than duplicating the aggregation logic across modules.
log_message "**** Donor x cell-type pseudobulk for the local RAS panel ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../../basement_membrane/_h/02.niche_pseudobulk.py \
       --adata ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ../_h/ras_panel.tsv \
       --outfile ./ras_pseudobulk_celltype.tsv.gz
if [ $? -ne 0 ]; then log_message "Error: pseudobulk failed"; exit 1; fi
conda deactivate

log_message "**** RAS source / receiver landscape and circuit completeness ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/01.ras_landscape_stats.R \
        --pseudobulk ./ras_pseudobulk_celltype.tsv.gz \
        --panels ../_h/ras_panel.tsv \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
