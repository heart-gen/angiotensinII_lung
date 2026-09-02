#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=supp_tables
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/supp_tables.log

## Builds supplementary Tables S1-S14 into tables/_m/tsv/*.tsv and
## tables/_m/supplementary_tables.xlsx. Submit from tables/_m:
##     sbatch -D tables/_m tables/_h/step_tables.sh
##
## UPSTREAM DEPENDENCIES -- run these first, or parts will be flagged
## `pending_upstream` in the manifest and 08 will report them:
##   sbatch -D pericyte_states/_m    pericyte_states/_h/step_0c.sh  # per-cluster Jaccard
##   sbatch -D pericyte_states/_m    pericyte_states/_h/step_1.sh   # both donor thresholds
##   sbatch -D cell_communication/_m cell_communication/_h/step_4.sh # Fig 4D lm + TGFB2
##   sbatch -D niche_index/_m        niche_index/_h/step_0.sh
##   sbatch -D niche_index/_m        niche_index/_h/step_1.sh
##   sbatch -D sensitivity/_m        sensitivity/_h/step_0.sh
## Plus the full ligand x BM-target matrix for S11C, which runs in seconds and
## needs no SLURM job because --links-only skips the permutation null entirely:
##   (cd basement_membrane/_m && Rscript ../_h/07.bm_nichenet_targets.R --links-only)
##
## 08.assemble_tables.R re-checks a set of manuscript values against the built
## tables and FAILS the job on any mismatch, so a silent disagreement between the
## supplement and the text cannot ship.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

## S3 needs the h5ad to score panel-gene detection, so it runs in the Python env.
log_message "**** S3: curated program gene sets (scRNA_env) ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/02.genesets.py --root ../.. --outdir ./tsv
if [ $? -ne 0 ]; then log_message "Error: gene-set table failed"; exit 1; fi
conda deactivate

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** S1, S2: cohort characteristics and mouse cross-species ****"
Rscript ../_h/01.cohort_mouse.R
if [ $? -ne 0 ]; then log_message "Error: cohort/mouse tables failed"; exit 1; fi

log_message "**** S4, S5: cluster characteristics and AGTR1 models ****"
Rscript ../_h/03.clusters_models.R
if [ $? -ne 0 ]; then log_message "Error: cluster/model tables failed"; exit 1; fi

log_message "**** S6: ACTA2 benchmark and CoGAPS validation ****"
Rscript ../_h/04.validation.R
if [ $? -ne 0 ]; then log_message "Error: validation tables failed"; exit 1; fi

log_message "**** S7, S8, S9: correlations and basement membrane ****"
Rscript ../_h/05.correlations_bm.R
if [ $? -ne 0 ]; then log_message "Error: correlation/BM tables failed"; exit 1; fi

log_message "**** S10, S11: NicheNet injury and BM target sets ****"
Rscript ../_h/06.nichenet.R
if [ $? -ne 0 ]; then log_message "Error: NicheNet tables failed"; exit 1; fi

log_message "**** S12, S13, S14: local RAS, disease, robustness ****"
Rscript ../_h/07.ras_disease.R
if [ $? -ne 0 ]; then log_message "Error: RAS/disease tables failed"; exit 1; fi

## Manifest and workbook last: it verifies every registered part exists and that
## the tables still reproduce the manuscript's quoted values.
log_message "**** Anchor checks + workbook ****"
Rscript ../_h/08.assemble_tables.R
if [ $? -ne 0 ]; then log_message "Error: assembly/anchor checks failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
