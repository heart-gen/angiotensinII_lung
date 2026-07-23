#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agt_vs_ligands
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/agt_vs_ligands.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** AGT versus the other pericyte ligands ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/02.agt_vs_ligands.R \
        --activities ../../cell_communication/_m/nichenet/ligand_activities_Pericytes.tsv \
        --links ../../cell_communication/_m/nichenet/ligand_target_links_Pericytes.tsv \
        --pseudobulk ./ras_pseudobulk_celltype.tsv.gz \
        --priors ../../cell_communication/_m/nichenet_priors \
        --frac-file ../../cell_communication/_m/expressed_fraction_main.tsv.gz \
        --outdir ./stats_data \
        --nboot 300
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
