#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=program_category
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
## Program x protein-category enrichment from the pericyte_states h5ad (replaces
## the uninformative state->category alluvial flow with a quantitative readout).
#SBATCH --time=01:00:00
#SBATCH --output=logs/program_category.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** program x category enrichment ****"
python ../_h/03b.program_category_enrichment.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: Python execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
