#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_run
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
## De-novo rank selection: sweep nPatterns x seeds. One array task per seed, each
## running the full nP sweep (~7 fits x ~8-12 min = ~1.5 h/task). Seed 13 is the
## canonical (untagged) run reused downstream; the other seeds are tagged so they
## don't clobber it and feed the cross-seed stability estimate in 04.
#SBATCH --array=0-3
#SBATCH --time=06:00:00
#SBATCH --output=logs/cogaps_run_seed%a.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID} (array task ${SLURM_ARRAY_TASK_ID})"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

NP_SWEEP="4,5,6,7,8,9,10"
SEEDS=(13 1 42 2024)
S=${SEEDS[$SLURM_ARRAY_TASK_ID]}
if [ "$S" -eq 13 ]; then TAG=""; else TAG="_seed${S}"; fi   # seed 13 = canonical

log_message "**** CoGAPS sweep nPatterns=${NP_SWEEP}, seed=${S}, tag='${TAG}' ****"
Rscript ../_h/01.run_cogaps.R \
        --indir ../_m --outdir ../_m \
        --npatterns ${NP_SWEEP} --niterations 5000 --nsets 8 --nthreads 1 \
        --seed ${S} --tag "${TAG}"

if [ $? -ne 0 ]; then log_message "Error: CoGAPS run (seed ${S}) failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
