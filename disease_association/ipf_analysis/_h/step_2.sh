#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ipf_dataset
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=conversion.log
##
## ** DO NOT RETIRE THIS SCRIPT. ** Everything else in `ipf_analysis/` is retired
## 2024 analysis (see ../RETIRED.md), but `02.generate_ipf_data.R` builds
## `_m/ipf_dataset.h5ad` (8.4 GB, 312,928 x 45,947, raw counts), which four
## CURRENT modules read:
##   * basement_membrane/_h/step_3.sh              (COPD/IPF arm)
##   * disease_association/agtr1_copd_ipf/         (step_1.sh, step_1b.sh, 02)
##   * disease_association/pericyte_analysis/      (01.subset_data.py, 03.transfer_labels.py)
##   * tables/_h/01.cohort_mouse.R                 (Table S1, via ../_h/sample_demo.csv)
## `../_h/sample_demo.csv` is likewise live. This is defect P1-18.

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

Rscript ../_h/02.generate_ipf_data.R

if [ $? -ne 0 ]; then
    log_message "Error: Rscript execution failed (02.generate_ipf_data.R)"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
