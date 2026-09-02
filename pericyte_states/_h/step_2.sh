#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_continuum
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=logs/continuum.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

## The batch shell does not source the login profile, so `module` and `conda`
## must be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Expression-only state continuum (DPT + PAGA) ****"
python ../_h/02.continuum_dpt.py \
       --adata pericyte_states.h5ad --outdir "./" \
       --denoise ../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv \
       --den-model Pericyte-only-trained

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
