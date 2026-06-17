#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=cogaps_stability
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
## three extra CoGAPS seeds at the chosen nP=5 (~6-8 min each) + stability calc.
#SBATCH --time=04:00:00
#SBATCH --output=logs/cogaps_stability.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## (1) re-run CoGAPS at nP=5 for three extra seeds, tagged so the canonical
## seed-13 outputs are preserved.
for S in 1 42 2024; do
    log_message "**** CoGAPS nP=5 seed=${S} ****"
    Rscript ../_h/01.run_cogaps.R \
            --indir ../_m --outdir ../_m \
            --npatterns 5 --niterations 5000 --nsets 8 --nthreads 1 \
            --seed ${S} --tag _seed${S}
    if [ $? -ne 0 ]; then log_message "Error: CoGAPS seed ${S} failed"; exit 1; fi
done

## (2) quantify seed stability + nP correspondence.
log_message "**** CoGAPS stability summary ****"
Rscript ../_h/04.cogaps_stability.R \
        --indir ../_m --outdir ../_m \
        --ref-np 5 --seed-tags 1,42,2024 --np-set 7,9
if [ $? -ne 0 ]; then log_message "Error: stability summary failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
