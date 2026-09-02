#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=agtr1_counts
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=logs/agtr1_counts.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "User: ${USER}"; echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

## The batch shell does not source the login profile, so `module` and `conda`
## must be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

log_message "**** Extract integer AGTR1 counts (raw/X, not the soupX matrix) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/09.agtr1_counts_prep.py \
       --raw-adata ../../localization/pericyte_analysis/_m/pericyte.hlca_full.dataset.h5ad \
       --metadata ./bm_metadata.tsv.gz \
       --out ./agtr1_count_input.tsv.gz
if [ $? -ne 0 ]; then log_message "Error: count extraction failed"; exit 1; fi
conda deactivate

log_message "**** AGTR1 count-model lens: NB / binomial GLMM, no imputation ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/10.agtr1_count_models.R \
        --input ./agtr1_count_input.tsv.gz \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: R failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
