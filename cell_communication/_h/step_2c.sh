#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=ccc_nichenet_cogaps
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
## NicheNet on the ORTHOGONAL (CoGAPS) receiver scheme + the cross-scheme ligand
## concordance table. step_5.sh added the CoGAPS receivers to the niche and ran
## LIANA on them, but never ran NicheNet, so the "prioritized ligands do not
## depend on how receivers are defined" claim had no data-driven receiver arm.
#SBATCH --time=02:00:00
#SBATCH --output=logs/ccc_nichenet_cogaps.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

## Receiver labels are the dominant-pattern groups written into
## expressed_fraction_cogaps.tsv.gz by 00b.cogaps_receivers.py (nP=8, main rank).
log_message "**** NicheNet: CoGAPS-defined pericyte receivers (nP=8) ****"
Rscript ../_h/02.nichenet.R \
        --priors ../_m/nichenet_priors --liana-dir ../_m \
        --frac-file expressed_fraction_cogaps.tsv.gz \
        --receivers Pericyte_cg_vascular_stabilizing,Pericyte_cg_basement_membrane,Pericyte_cg_fibroblast_like,Pericyte_cg_inflammatory,Pericyte_cg_synthetic_contractile \
        --outdir ../_m/nichenet_cogaps
if [ $? -ne 0 ]; then log_message "Error: NicheNet cogaps failed"; exit 1; fi

log_message "**** Cross-receiver ligand concordance ****"
Rscript ../_h/06.receiver_concordance.R
if [ $? -ne 0 ]; then log_message "Error: receiver concordance failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
