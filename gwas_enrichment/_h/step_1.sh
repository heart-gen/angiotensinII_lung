#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gwas_munge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x.%j.log

## Harmonize staged GWAS to rsID-keyed MAGMA + LDSC inputs.
## Submit:  sbatch -D ../_m ../_h/step_1.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** munge GWAS ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
mkdir -p logs
python ../_h/01.munge_gwas.py --outdir ./
log "**** done ****"
