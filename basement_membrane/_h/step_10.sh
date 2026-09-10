#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=agtr1_lens_cluster
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
## Three lmer fits on 214 pseudobulk rows. Minutes, not hours -- unlike step_6.sh,
## nothing here re-reads an h5ad or fits an NB GLMM on 11,680 cells.
#SBATCH --time=00:20:00
#SBATCH --output=logs/agtr1_lens_by_cluster.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## MUST run after step_6.sh: it consumes that step's agtr1_count_input.tsv.gz and
## cross-checks its pseudobulk unit count against agtr1_count_by_cluster.tsv.
log_message "**** AGTR1 lenses refit at the count model's pseudobulk unit ****"
Rscript ../_h/19.agtr1_lens_by_cluster.R \
        --input ./agtr1_count_input.tsv.gz \
        --state-meta ../../pericyte_states/_m/pericytes_states_metadata.tsv.gz \
        --denoise ../../localization/airspace_analysis/_m/airspace/pericytes_airspace_denoising.tsv \
        --den-model Pericyte-only-trained \
        --count-table ./stats_data/agtr1_count_by_cluster.tsv \
        --outdir ./stats_data \
        --min-cells 5
if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
