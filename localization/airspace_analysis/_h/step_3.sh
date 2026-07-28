#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_dropout
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=logs/agtr1_dropout.log

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

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

## Reads the pericyte object written by pericyte_states/step_0.sh together with the
## denoising table written by step_2.sh here; it recomputes neither. Everything it
## writes lands in ./agtr1_dropout/ and feeds figures/_h/agtr1_dropout_figure.R.
log_message "**** AGTR1 dropout model + airspace ladder ****"
python ../_h/03.agtr1_dropout.py \
       --raw-adata ../../pericyte_analysis/_m/pericyte.hlca_full.dataset.h5ad \
       --state-adata ../../../pericyte_states/_m/pericyte_states.h5ad \
       --qc-adata ./pericytes_with_airspace_score.h5ad \
       --denoise ./airspace/pericytes_airspace_denoising.tsv \
       --outdir ./agtr1_dropout

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
