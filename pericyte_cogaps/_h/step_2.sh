#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=cogaps_select
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/cogaps_select.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## De-novo rank selection from the nP x seed sweep produced by step_1.sh:
## cross-seed pattern robustness + reconstruction error -> recommended nP.
## Inspect cogaps_nP_selection.tsv / figures/cogaps_nP_selection.* and set NP in
## step_3.sh (validate) and step_5.sh (project) accordingly.
log_message "**** CoGAPS de-novo rank selection ****"
Rscript ../_h/02.select_rank.R \
        --indir ../_m --outdir ../_m \
        --np-sweep 4,5,6,7,8,9,10 --seed-tags 1,42,2024 \
        --stability-threshold 0.8
if [ $? -ne 0 ]; then log_message "Error: rank selection failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
