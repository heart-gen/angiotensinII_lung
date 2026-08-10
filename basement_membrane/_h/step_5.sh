#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=bm_state_axis
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --output=logs/bm_state_axis.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** BM vs pericyte clusters, fibrillar axis, AGTR1 and injury continuum ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
## --denoise supplies the scVI-denoised AGTR1 for the direct AGTR1-vs-BM test.
## Same table/model as pericyte_states/_h/03.agtr1_lenses.R: raw AGTR1 and the BM
## panel share a depth confound, so the denoised lens is the discriminating one.
Rscript ../_h/04.bm_state_stats.R \
        --bm-meta ./bm_metadata.tsv.gz \
        --state-meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
        --continuum ../../pericyte_states/_m/continuum_metadata.tsv.gz \
        --denoise ../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv \
        --den-model Pericyte-only-trained \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
