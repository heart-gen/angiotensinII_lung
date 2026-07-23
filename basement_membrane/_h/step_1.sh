#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=bm_state_gate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/bm_state_gate.log

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

log_message "**** Pre-specified gate: does a BM panel change the state model? ****"
python ../_h/01.state_gate.py \
       --state-meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
       --bm-meta ./bm_metadata.tsv.gz \
       --reference-map ../../pericyte_states/_m/annotations/state_program_map.tsv \
       --outdir "./"

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
