#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=gwas_specificity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=02:00:00
## 64 cpus x 2000M = 128G: the HLCA core h5ad (5.5G on disk) expands well past 48G when
## read fully into memory; an earlier 24-cpu (48G) run was OOM-killed during the load.
#SBATCH --output=logs/%x.%j.log

## Build autosomal cell-type / pericyte-state specificity + LDSC BEDs.
## Submit:  sbatch -D ../_m ../_h/step_0.sh
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** specificity ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
mkdir -p logs
python ../_h/00.build_specificity.py --outdir ./
log "**** done ****"
