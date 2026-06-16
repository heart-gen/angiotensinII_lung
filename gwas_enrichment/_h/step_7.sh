#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gwas_summary
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.%j.log

## Collate MAGMA + S-LDSC + coloc into the enrichment summary and figures.
## Submit:  sbatch -D ../_m ../_h/step_7.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** summary + figures ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
mkdir -p logs
Rscript ../_h/07.summarize_and_figures.R
log "**** done ****"
