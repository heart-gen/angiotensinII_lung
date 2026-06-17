#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=EM
#SBATCH --job-name=ccc_prepare
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=04:00:00
#SBATCH --output=logs/ccc_prepare.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Build CCC niche from full disease object ****"
python ../_h/00.prepare_ccc_input.py \
       --adata ../../disease_association/_m/hlca_full.dataset.h5ad \
       --pericyte-states ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
       --outdir "./"

if [ $? -ne 0 ]; then log_message "Error: Python failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
