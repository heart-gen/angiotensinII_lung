#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_pseudobulk
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
## pseudobulk read of the 423k-cell niche onto the CoGAPS gene space; modest mem.
#SBATCH --time=01:30:00
#SBATCH --output=logs/cogaps_pseudobulk.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

## Pseudobulk the niche onto the CoGAPS gene space (prereq for step_5 projection).
log_message "**** Niche pseudobulk (scRNA_env) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/04.niche_pseudobulk.py \
       --niche ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ../_m/cogaps_genes.tsv \
       --outdir ../_m --min-cells 10
if [ $? -ne 0 ]; then log_message "Error: pseudobulk failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
