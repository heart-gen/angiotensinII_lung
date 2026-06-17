#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-small
#SBATCH --job-name=niche_index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/niche_index.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

module purge
module load anaconda3/2024.10-1
module list
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Build donor-level niche index ****"
python ../_h/00.niche_index.py \
       --states-meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
       --airspace-summary ../../localization/airspace_analysis/_m/airspace/airspace_donor_summary.csv \
       --outdir "./"

if [ $? -ne 0 ]; then log_message "Error: Python failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
