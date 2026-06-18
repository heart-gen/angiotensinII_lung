#!/bin/bash
#SBATCH --account=bio260021p
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

log_message "**** Supplementary sensitivity/robustness figure ****"
Rscript ../_h/sensitivity_robustness_figure.R
if [ $? -ne 0 ]; then log_message "Error: sensitivity figure failed"; exit 1; fi

log_message "**** Panel manifest ****"
Rscript ../_h/assemble_mechanism_figures.R
if [ $? -ne 0 ]; then log_message "Error: manifest failed"; exit 1; fi

log_message "**** Manuscript figures (main + CCC/NicheNet + supplements) ****"
Rscript ../_h/manuscript_mechanism_figure.R
if [ $? -ne 0 ]; then log_message "Error: manuscript figures failed"; exit 1; fi
conda deactivate
log_message "**** Job ends ****"
