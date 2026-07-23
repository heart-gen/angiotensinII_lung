#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=cogaps_validate
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/cogaps_validate.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## Validate the selected rank against the curated state model. Set NPATTERNS to
## the nP chosen de novo by step_2 (cogaps_nP_selection.tsv); keep neighbours if
## you want to show the choice is robust to +/-1.
NPATTERNS="5"   # <- update to the step_2 recommendation before running
for NP in ${NPATTERNS}; do
    log_message "**** Validate CoGAPS nPatterns=${NP} ****"
    Rscript ../_h/03.cogaps_validate.R \
            --indir ../_m \
            --meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
            --hvg-info ../_m/cogaps_hvg_info.tsv.gz \
            --npatterns ${NP} --top-markers 50 \
            --outdir ../_m/validation_np${NP}
    if [ $? -ne 0 ]; then log_message "Error: validate nP=${NP} failed"; exit 1; fi
done
conda deactivate
log_message "**** Job ends ****"
