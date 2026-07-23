#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_copd
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=logs/bm_copd.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Donor x compartment pseudobulk (GSE136831) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/05.bm_copd.py \
       --adata ../../disease_association/ipf_analysis/_m/ipf_dataset.h5ad \
       --genes ./bm_panel_genes.tsv \
       --demo ../../disease_association/ipf_analysis/_h/sample_demo.csv \
       --outfile ./gse136831_bm_pseudobulk.tsv.gz
if [ $? -ne 0 ]; then log_message "Error: pseudobulk failed"; exit 1; fi
conda deactivate

log_message "**** COPD basement-membrane contrast ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/06.bm_copd_stats.R \
        --pseudobulk ./gse136831_bm_pseudobulk.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
