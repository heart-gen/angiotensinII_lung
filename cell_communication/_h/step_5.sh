#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ccc_cogaps
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
## single (cogaps) receiver scheme; ~ one disease-stratified liana pass.
#SBATCH --time=03:00:00
#SBATCH --output=logs/ccc_cogaps.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

## Orthogonal data-driven receiver scheme: assign pericytes to their dominant
## CoGAPS pattern and add it to the CCC niche. nP=8 is the MAIN rank; nP=9 is the
## SENSITIVITY rank (see pericyte_cogaps/_m/cogaps_nP_selection.tsv). Both need their
## validation_np${NP}/ produced first (pericyte_cogaps step_3).
NP_MAIN=8
NP_SENS=9

## Sensitivity rank: write ONLY the annotation table (--write-h5ad false) so the
## dominant-pattern receiver column in ccc_niche.h5ad keeps the MAIN scheme. The
## sensitivity projection (pericyte_cogaps step_5) reads this annotation directly.
log_message "**** CoGAPS receiver annotation, sensitivity nPatterns=${NP_SENS} (annotation only) ****"
python ../_h/00b.cogaps_receivers.py \
       --niche ccc_niche.h5ad \
       --patterns ../../pericyte_cogaps/_m/patterns_np${NP_SENS}.tsv.gz \
       --score-spearman ../../pericyte_cogaps/_m/validation_np${NP_SENS}/pattern_score_spearman.tsv.gz \
       --npatterns ${NP_SENS} --min-rho 0.15 --write-h5ad false --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: cogaps receivers (sens) failed"; exit 1; fi

## Main rank: annotation + write ccc_group_cogaps into ccc_niche.h5ad (last h5ad
## write wins, so the niche carries the MAIN scheme) + the disease-stratified LIANA.
log_message "**** Add CoGAPS receiver scheme, main nPatterns=${NP_MAIN} ****"
python ../_h/00b.cogaps_receivers.py \
       --niche ccc_niche.h5ad \
       --patterns ../../pericyte_cogaps/_m/patterns_np${NP_MAIN}.tsv.gz \
       --score-spearman ../../pericyte_cogaps/_m/validation_np${NP_MAIN}/pattern_score_spearman.tsv.gz \
       --npatterns ${NP_MAIN} --min-rho 0.15 --outdir "./"
if [ $? -ne 0 ]; then log_message "Error: cogaps receivers (main) failed"; exit 1; fi

log_message "**** Disease-stratified liana into CoGAPS receivers (main nP=${NP_MAIN}) ****"
python ../_h/01.run_liana.py --adata ccc_niche.h5ad --outdir "./" --schemes cogaps
if [ $? -ne 0 ]; then log_message "Error: liana cogaps failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
