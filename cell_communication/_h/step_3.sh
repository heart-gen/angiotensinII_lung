#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-small
#SBATCH --job-name=ccc_alluvial
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=logs/ccc_alluvial.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list

log_message "**** Marker detection fractions (python) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/03a.marker_fractions.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: 03a failed"; exit 1; fi
conda deactivate

log_message "**** Alluvial figures (R) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/03.figures.R
if [ $? -ne 0 ]; then log_message "Error: 03 failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
