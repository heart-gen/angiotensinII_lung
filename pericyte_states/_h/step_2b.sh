#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_continuum_sens
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
## ~28 DPT re-runs over a root/neighbors/n_dcs/subsample grid on 11.7k cells.
#SBATCH --time=03:00:00
#SBATCH --output=logs/continuum_sensitivity.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

## The batch shell does not source the login profile, so `module` and `conda`
## must be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

## --denoise/--den-model are passed explicitly, matching step_2.sh, so the sweep
## covers the same feature list (including AGTR1_scvi) as the headline run it is
## defending.
log_message "**** DPT continuum root/resolution sensitivity ****"
python ../_h/02b.continuum_sensitivity.py \
       --adata pericyte_states.h5ad --outdir "./" \
       --denoise ../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv \
       --den-model Pericyte-only-trained

if [ $? -ne 0 ]; then log_message "Error: Python execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
