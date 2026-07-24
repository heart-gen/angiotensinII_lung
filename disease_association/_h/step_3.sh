#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_dx_forest
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
## Pure donor-level TSV job (no h5ad / basilisk): study-adjusted LMM + random-
## effects meta-analysis forest for the Healthy-vs-Fibrotic/ILD pericyte
## injury-program contrast, plus smoking + min-cells sensitivity. RM-shared
## controls memory via cores (2000 MB/core); do NOT set --mem.
#SBATCH --time=00:20:00
#SBATCH --output=logs/disease_forest.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript ../_h/03.disease_forest.R \
        --meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
        --outdir mixed_model_forest \
        --min-cells 10

if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
