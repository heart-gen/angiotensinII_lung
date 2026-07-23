#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ccc_donor_val
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=logs/ccc_donor_val.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list

log_message "**** Donor-level NicheNet validation: extract per-donor expression ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/04.donor_validation.py --adata ccc_niche.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: 04 failed"; exit 1; fi
conda deactivate

log_message "**** Donor-level association test ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/05.donor_validation_stats.R .
if [ $? -ne 0 ]; then log_message "Error: 05 failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
