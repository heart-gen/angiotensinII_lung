#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=assemble_figures
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/assemble_figures.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"
module purge; module load anaconda3/2024.10-1; module list

# Submitted from figures/_m so getwd()/../.. resolves to project root.
log_message "**** Export pericyte UMAP coords (scRNA_env) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/00.export_pericyte_umap.py \
       --adata ../../pericyte_states/_m/pericyte_states.h5ad --outdir ./
if [ $? -ne 0 ]; then log_message "Error: UMAP export failed"; exit 1; fi
conda deactivate

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** Integrated pericyte-layer figure (main + supplement) ****"
Rscript ../_h/pericyte_layer_figure.R
if [ $? -ne 0 ]; then log_message "Error: pericyte-layer figure failed"; exit 1; fi

log_message "**** S3: pericyte-state stability + annotation ****"
Rscript ../_h/state_annotation_figure.R
if [ $? -ne 0 ]; then log_message "Error: state annotation figure failed"; exit 1; fi

log_message "**** S4: AGTR1 dropout vs distinct population ****"
Rscript ../_h/agtr1_dropout_figure.R
if [ $? -ne 0 ]; then log_message "Error: AGTR1 dropout figure failed"; exit 1; fi

log_message "**** S6: unsupervised CoGAPS validation of the programs ****"
Rscript ../_h/cogaps_validation_figure.R
if [ $? -ne 0 ]; then log_message "Error: CoGAPS validation figure failed"; exit 1; fi

log_message "**** S12: sensitivity / robustness ****"
Rscript ../_h/sensitivity_robustness_figure.R
if [ $? -ne 0 ]; then log_message "Error: sensitivity figure failed"; exit 1; fi

log_message "**** Basement-membrane figure + S10, S15 ****"
Rscript ../_h/basement_membrane_figure.R
if [ $? -ne 0 ]; then log_message "Error: basement-membrane figure failed"; exit 1; fi

log_message "**** S7: pericyte continuum stability ****"
Rscript ../_h/continuum_stability_figure.R
if [ $? -ne 0 ]; then log_message "Error: continuum stability figure failed"; exit 1; fi

log_message "**** S8: NicheNet specificity + donor validation ****"
Rscript ../_h/nichenet_specificity_figure.R
if [ $? -ne 0 ]; then log_message "Error: nichenet specificity figure failed"; exit 1; fi

log_message "**** S9: receiver robustness + BM-restricted signaling ****"
Rscript ../_h/receiver_robustness_figure.R
if [ $? -ne 0 ]; then log_message "Error: receiver robustness figure failed"; exit 1; fi

log_message "**** S11: discrete state composition by disease ****"
Rscript ../_h/state_composition_figure.R
if [ $? -ne 0 ]; then log_message "Error: state composition figure failed"; exit 1; fi

log_message "**** Manuscript figures (main + CCC/NicheNet + supplements) ****"
Rscript ../_h/manuscript_mechanism_figure.R
if [ $? -ne 0 ]; then log_message "Error: manuscript figures failed"; exit 1; fi

## Manifest last: it records which supplements actually exist on disk.
log_message "**** Panel manifest ****"
Rscript ../_h/assemble_mechanism_figures.R
if [ $? -ne 0 ]; then log_message "Error: manifest failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
