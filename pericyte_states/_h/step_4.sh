#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_acta2_control
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
## CONTROL: does AGTR1 just track ACTA2+ contractile identity? Recompute the
## leave-ACTA2-out contractile score from the h5ad (Python), then run the
## donor-aware ACTA2 benchmark + AGTR1-vs-contractile correlations (R).
#SBATCH --time=01:00:00
#SBATCH --output=logs/acta2_control.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** (1) leave-ACTA2-out contractile score (Python) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/04.acta2_control.py --adata pericyte_states.h5ad --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: Python step failed"; exit 1; fi
conda deactivate

log_message "**** (2) ACTA2 control stats + correlations (R) ****"
export BASILISK_EXTERNAL_DIR=/ocean/projects/bio260021p/shared/opt/basilisk_cache
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/05.acta2_control.R \
        --meta pericytes_states_metadata.tsv.gz \
        --noacta2 synth_contr_noACTA2.tsv.gz \
        --outdir stats_data
if [ $? -ne 0 ]; then log_message "Error: R step failed"; exit 1; fi
conda deactivate

log_message "**** Job ends ****"
