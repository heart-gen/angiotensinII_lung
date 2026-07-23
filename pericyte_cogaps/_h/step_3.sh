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

## Validate the selected rank(s) against the curated state model. step_2 chose the
## sweep; we carry TWO ranks forward: nP=8 is the MAIN rank (most reproducible band,
## min_r 0.978) and nP=9 is the SENSITIVITY rank (largest with min_r >= 0.80). Both
## get a validation_np${NP}/ dir so the choice is shown robust across the stable band.
NPATTERNS="8 9"   # 8 = main, 9 = sensitivity (see cogaps_nP_selection.tsv)
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
