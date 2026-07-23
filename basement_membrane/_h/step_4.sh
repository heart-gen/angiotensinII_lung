#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_nichenet
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/bm_nichenet.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** NicheNet ligand ranking toward basement-membrane targets ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/07.bm_nichenet_targets.R \
        --priors ../../cell_communication/_m/nichenet_priors \
        --liana-dir ../../cell_communication/_m \
        --frac-file expressed_fraction_main.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --frozen-activities ../../cell_communication/_m/nichenet/ligand_activities_Pericytes.tsv \
        --outdir ./nichenet_bm \
        --nperm 10000
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
