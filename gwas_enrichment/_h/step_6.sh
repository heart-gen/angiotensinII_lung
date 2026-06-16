#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=coloc_agtr1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x.%j.log

## AGTR1-locus colocalization (COPD GWAS x lung cis-eQTL).
## coloc is installed in the shared R_env (the script prepends .Rlib but keeps R_env's
## site lib on the path, so coloc resolves there). No per-user install needed.
## Submit:  sbatch -D ../_m ../_h/step_6.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** coloc AGTR1 ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
mkdir -p logs
Rscript ../_h/06.coloc_agtr1.R
log "**** done ****"
