#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=nichenet_specificity
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
## Gene-set permutation null for the primary pericyte receiver.
## RAISED 1,000 -> 10,000 perms 2026-09-07 (defect P1-6). At 1,000 the empirical
## p is floored at 1/1001 for every top ligand, so the test cannot ORDER them --
## which is precisely the question P1-6 asks (TGFB2 ranks first by AUPR, TGFB1
## first by z). 1,000 perms took ~19 min in the 2026-06-16 run, so 10,000 is ~3.2 h.
#SBATCH --time=05:00:00
#SBATCH --output=logs/nichenet_specificity.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** NicheNet specificity: gene-set permutation null (Pericytes) ****"
Rscript ../_h/02b.nichenet_specificity.R \
        --priors ../_m/nichenet_priors --liana-dir ../_m \
        --frac-file expressed_fraction_main.tsv.gz \
        --receiver Pericytes --n-perm 10000 --top-ligands 15 \
        --outdir ../_m/nichenet
if [ $? -ne 0 ]; then log_message "Error: NicheNet specificity failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
