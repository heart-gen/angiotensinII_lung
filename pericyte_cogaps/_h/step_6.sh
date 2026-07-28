#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=cogaps_concord
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --output=logs/cogaps_concordance.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## Main (nP=8) vs sensitivity (nP=9) rank concordance. Needs step_1 (the fits),
## step_3 (validation_np8/ and validation_np9/) to have run.
log_message "**** CoGAPS nP=8 vs nP=9 concordance ****"
Rscript ../_h/06.rank_concordance.R \
        --indir ../_m --outdir ../_m --np-main 8 --np-sens 9
if [ $? -ne 0 ]; then log_message "Error: rank concordance failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
