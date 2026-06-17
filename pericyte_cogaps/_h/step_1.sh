#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_run
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
## single-cell distributed CoGAPS over a sweep of nPatterns; give it head room.
#SBATCH --time=12:00:00
#SBATCH --output=logs/cogaps_run.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** CoGAPS pattern sweep (nPatterns 5,7,9) ****"
Rscript ../_h/01.run_cogaps.R \
        --indir ../_m --outdir ../_m \
        --npatterns 5,7,9 --niterations 5000 --nsets 8 --nthreads 1 --seed 13

if [ $? -ne 0 ]; then log_message "Error: CoGAPS run failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
