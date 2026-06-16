#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ldsc_ldscores
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000M
#SBATCH --array=1-22
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x.%A_%a.log

## ===========================================================================
## S-LDSC step 2/3 — compute custom LD scores for every annotation set, one
## SLURM array task per chromosome.  Submit:  sbatch -D ../_m ../_h/04.compute_ldscores.sh
## ===========================================================================
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
CHR=${SLURM_ARRAY_TASK_ID}
log "**** ldscores chr${CHR} ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics
mkdir -p logs

LDSC_DIR=/ocean/projects/bio250020p/shared/opt/ldsc
RES=/ocean/projects/bio250020p/shared/resources/ldsc
BIM_DIR=${RES}/1000G_EUR_Phase3_plink
HM3=${RES}/hm3_no_MHC.list.txt
WRAP=../_h/ldsc_wrapper.py

mapfile -t ANNOT_DIRS < <(find ./annot -mindepth 1 -maxdepth 1 -type d | sort)
[[ ${#ANNOT_DIRS[@]} -gt 0 ]] || { log "ERROR: no annot dirs; run 03.make_annot.sh"; exit 1; }

for AD in "${ANNOT_DIRS[@]}"; do
    KEY=$(basename "$AD")
    AF="${AD}/${KEY}.${CHR}.annot.gz"
    [[ -f "$AF" ]] || { log "  missing ${AF}"; continue; }
    OUT=./ldscores/${KEY}; mkdir -p "$OUT"
    PFX="${OUT}/${KEY}.${CHR}"
    [[ -f "${PFX}.l2.ldscore.gz" ]] && continue
    python "$WRAP" "$LDSC_DIR" ldsc.py --l2 \
        --bfile "${BIM_DIR}/1000G.EUR.QC.${CHR}" --ld-wind-cm 1 \
        --annot "$AF" --thin-annot --print-snps "$HM3" --out "$PFX"
done
conda deactivate
log "**** ldscores chr${CHR} done ****"
