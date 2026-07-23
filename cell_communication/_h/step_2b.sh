#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=nichenet_specificity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
## gene-set permutation null (1000 perms) for the primary pericyte receiver.
#SBATCH --time=02:00:00
#SBATCH --output=logs/nichenet_specificity.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** NicheNet specificity: gene-set permutation null (Pericytes) ****"
Rscript ../_h/02b.nichenet_specificity.R \
        --priors ../_m/nichenet_priors --liana-dir ../_m \
        --frac-file expressed_fraction_main.tsv.gz \
        --receiver Pericytes --n-perm 1000 --top-ligands 15 \
        --outdir ../_m/nichenet
if [ $? -ne 0 ]; then log_message "Error: NicheNet specificity failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
