#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=disease_main_fig
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
## Disease main figure (continuous-injury version): graded three-group endpoint,
## component programs, AGTR1 across stromal cell types, and the independent
## GSE136831 COPD/IPF evaluation -- plus the two-panel study-level robustness
## supplement (leave-one-study-out + within-study forest) as
## figureS_disease_robustness. Reads disease_association/_m/
## {mixed_model_forest,mean_expr} and disease_association/agtr1_copd_ipf/_m, so
## run disease_association step_3 and step_5 and agtr1_copd_ipf step_1 first.
## Pure-TSV -> ggplot; RM-shared controls memory via cores (2000 MB/core), do NOT
## set --mem. Submit from figures/_m so getwd()/../.. resolves to the project root.
#SBATCH --time=00:15:00
#SBATCH --output=logs/disease_main_figure.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

module purge
module load anaconda3/2024.10-1
module list

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

Rscript ../_h/disease_main_figure.R
if [ $? -ne 0 ]; then log_message "Error: disease main figure failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
