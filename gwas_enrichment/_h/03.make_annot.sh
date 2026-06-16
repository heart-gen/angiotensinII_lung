#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ldsc_make_annot
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x.%j.log

## ===========================================================================
## S-LDSC step 1/3 — build per-chromosome binary annotations from the cell-type /
## pericyte-state gene BEDs written by 00.build_specificity.py (hg19, ±100kb).
## Submit:  sbatch -D ../_m ../_h/03.make_annot.sh
## ===========================================================================
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** make_annot ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics
## make_annot.py uses pybedtools (.sort/.merge/.intersect), which shells out to the
## legacy sortBed/mergeBed/intersectBed names. The PSC bedtools module ships only the
## singularity-wrapped `bedtools`, so prepend small shims that forward to `bedtools <cmd>`.
module load bedtools/2.30.0
export PATH="$(cd ../_h/bedtools_shims && pwd):${PATH}"
mkdir -p logs

LDSC_DIR=/ocean/projects/bio250020p/shared/opt/ldsc
BIM_DIR=/ocean/projects/bio250020p/shared/resources/ldsc/1000G_EUR_Phase3_plink
WRAP=../_h/ldsc_wrapper.py
BED_DIR=./specificity/beds

mapfile -t BEDS < <(find "$BED_DIR" -name "*_hg19.bed" | sort)
[[ ${#BEDS[@]} -gt 0 ]] || { log "ERROR: no BEDs in ${BED_DIR}; run 00.build_specificity.py"; exit 1; }
log "Found ${#BEDS[@]} annotation BEDs"

for BED in "${BEDS[@]}"; do
    KEY=$(basename "$BED" _hg19.bed)
    OUT=./annot/${KEY}; mkdir -p "$OUT"
    for CHR in $(seq 1 22); do
        PFX="${OUT}/${KEY}.${CHR}"
        [[ -f "${PFX}.annot.gz" ]] && continue
        python "$WRAP" "$LDSC_DIR" make_annot.py \
            --bed-file "$BED" --bimfile "${BIM_DIR}/1000G.EUR.QC.${CHR}.bim" \
            --annot-file "${PFX}.annot.gz"
    done
    n=$(ls "${OUT}/"*.annot.gz 2>/dev/null | wc -l)
    log "  ${KEY}: ${n}/22"
done
conda deactivate
log "**** make_annot done ****"
