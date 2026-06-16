#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=sldsc_h2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000M
#SBATCH --array=0-12
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x.%A_%a.log

## ===========================================================================
## S-LDSC step 3/3 — one trait per array task: munge to .sumstats.gz, then run
## partitioned heritability for that trait against every cell-state annotation,
## conditioning on the baselineLD v2.2 model.  Reports tau* coefficient + p.
## Submit:  sbatch -D ../_m ../_h/05.partitioned_h2.sh   (array auto-skips unstaged traits)
## ===========================================================================
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
module purge; module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics
mkdir -p logs sumstats

LDSC_DIR=/ocean/projects/bio250020p/shared/opt/ldsc
RES=/ocean/projects/bio250020p/shared/resources/ldsc
BASE=${RES}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.
WEIGHTS=${RES}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=${RES}/1000G_Phase3_frq/1000G.EUR.QC.
HM3=${RES}/w_hm3.snplist
WRAP=../_h/ldsc_wrapper.py

# LDSC traits = those with a signed-Z munge input (everything except COPDGene |beta|).
TRAITS=(COPD_UKB IPF_UKB ASTHMA SMKINIT CIGDAY SBP DBP CAD EDU BODYFAT FEV1 FVC FEV1FVC)
idx="${SLURM_ARRAY_TASK_ID:?submit with --array=0-12}"
(( idx < ${#TRAITS[@]} )) || { log "idx ${idx} out of range"; exit 0; }
trait="${TRAITS[idx]}"
log "**** S-LDSC [${idx}]=${trait} ****"; echo "Job: ${SLURM_JOBID:-N/A}"

RAW="./ldsc_raw/${trait}.tsv"
[[ -f "$RAW" ]] || { log "  ${trait}: no LDSC input (not staged / unsigned); skip"; exit 0; }

SS="./sumstats/${trait}.sumstats.gz"
if [[ ! -f "$SS" ]]; then
    log "  munge ${trait}"
    python "$WRAP" "$LDSC_DIR" munge_sumstats.py \
        --sumstats "$RAW" --snp SNP --a1 A1 --a2 A2 --N-col N \
        --signed-sumstats Z,0 --merge-alleles "$HM3" \
        --out "./sumstats/${trait}"
fi

mapfile -t KEYS < <(find ./ldscores -mindepth 1 -maxdepth 1 -type d | xargs -n1 basename | sort)
[[ ${#KEYS[@]} -gt 0 ]] || { log "ERROR: no ldscores; run 04.compute_ldscores.sh"; exit 1; }

# --overlap-annot needs the .annot.gz beside the ldscores; symlink them in.
for k in "${KEYS[@]}"; do
    for chr in $(seq 1 22); do
        dst="./ldscores/${k}/${k}.${chr}.annot.gz"
        [[ -e "$dst" ]] || ln -s "../../annot/${k}/${k}.${chr}.annot.gz" "$dst" 2>/dev/null || true
    done
done

for k in "${KEYS[@]}"; do
    LDPFX="./ldscores/${k}/${k}."
    [[ -f "${LDPFX}1.l2.ldscore.gz" ]] || { log "  incomplete LD scores ${k}"; continue; }
    OUT=./results/${k}; mkdir -p "$OUT"
    [[ -f "${OUT}/${trait}.results" ]] && continue
    log "  h2: ${trait} x ${k}"
    python "$WRAP" "$LDSC_DIR" ldsc.py --h2 "$SS" \
        --ref-ld-chr "${BASE},${LDPFX}" --w-ld-chr "$WEIGHTS" --frqfile-chr "$FRQ" \
        --overlap-annot --thin-annot --print-coefficients \
        --out "${OUT}/${trait}" || log "    failed ${trait} x ${k}"
done
conda deactivate
log "**** S-LDSC [${trait}] done ****"
