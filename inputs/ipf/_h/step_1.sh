#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ipf_dataset
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=logs/conversion.log
##
## Builds `_m/ipf_dataset.h5ad` (8.4 GB, 312,928 x 45,947, raw counts) from the
## GEO GSE136831 download in `_m/`. Submit from inputs/ipf/_m:
##     sbatch -D inputs/ipf/_m inputs/ipf/_h/step_1.sh
##
## One module in this repository reads the output: `basement_membrane/_h/step_3.sh`,
## the COPD basement-membrane arm. See the script header for why this builder
## moved here from `disease_association/ipf_analysis/_h/` on 2026-09-10.

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
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell. Without
## these two lines this script cannot run under sbatch at all -- which mattered,
## because it is the ONLY builder of ipf_dataset.h5ad (see the header above).
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Run IPF conversion ****"

Rscript ../_h/01.generate_ipf_data.R

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed (02.generate_ipf_data.R)"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
