#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/cogaps_prep.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Prepare CoGAPS input (HLCA pericyte subset) ****"
python ../_h/00.prepare_cogaps_input.py \
        --adata ../../pericyte_states/_m/pericyte_states.h5ad \
        --outdir ../_m --n-hvg 2000 --seed 13

if [ $? -ne 0 ]; then log_message "Error: prep failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
