#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=ccc_nichenet
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --output=logs/ccc_nichenet.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** NicheNet: main (all pericytes + AGTR2-detectable AT2) ****"
Rscript ../_h/02.nichenet.R \
        --priors ../_m/nichenet_priors --liana-dir ../_m \
        --frac-file expressed_fraction_main.tsv.gz \
        --receivers Pericytes,AT2_AGTR2det \
        --outdir ../_m/nichenet
if [ $? -ne 0 ]; then log_message "Error: NicheNet main failed"; exit 1; fi

log_message "**** NicheNet: state-stratified pericyte receivers ****"
Rscript ../_h/02.nichenet.R \
        --priors ../_m/nichenet_priors --liana-dir ../_m \
        --frac-file expressed_fraction_state.tsv.gz \
        --receivers Pericyte_vascular_stabilizing,Pericyte_inflammatory,Pericyte_synthetic_contractile,Pericyte_activated_migratory,Pericyte_fibroblast_like \
        --outdir ../_m/nichenet_state
if [ $? -ne 0 ]; then log_message "Error: NicheNet state failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
