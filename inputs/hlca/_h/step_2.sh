#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=EM
#SBATCH --job-name=hlca_prep
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH --output=logs/preprocess_reference.log

## Builds `inputs/hlca/_m/hlca_full.dataset.h5ad` -- the normalized, HVG-selected
## all-cell HLCA object that `cell_communication/_h/00.prepare_ccc_input.py`
## reads. Submit from inputs/hlca/_m:
##     sbatch -D inputs/hlca/_m inputs/hlca/_h/step_2.sh
##
## This step used to live in `disease_association/_h/step_0.sh`. It moved here on
## 2026-09-10 when the disease analyses left this repository: despite its home,
## it was never a disease analysis -- the script's own docstring called it
## "reference data for localization workflow" -- and cell_communication, which
## stays, cannot run without its output.
##
## `--targets full` skips the stromal Harmony correction (75 iterations), which
## only the disease models need. See the script docstring for the second copy of
## this builder and the obligation to keep the two in sync.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Preprocess HLCA reference (all-cell object only) ****"
python ../_h/02.preprocess_reference.py --targets full
if [ $? -ne 0 ]; then log_message "Error: preprocessing failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
