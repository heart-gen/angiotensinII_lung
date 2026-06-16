#!/bin/bash
#SBATCH --account=bio260021p
#SBATCH --partition=RM-shared
#SBATCH --job-name=magma_geneprop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x.%j.log

## ===========================================================================
## MAGMA gene-property: is GWAS heritability enriched in genes specific to the
## AGTR1+ injury-niche pericyte state (and lung cell types)?
##   (A) annotate g1000_eur SNPs -> genes (hg19, once)
##   (B) per-trait gene analysis  (SNP-wise mean model)
##   (C) per-trait gene-property: regress gene z on each specificity column
##        - marginal  (one cell type / state per model; Bryois 2020 style)
##        - joint     (5 pericyte states together -> injury-state enrichment
##                     conditional on the other states)
## Run after 01.munge_gwas.py.  Submit:  sbatch -D ../_m ../_h/02.run_magma.sh
## ===========================================================================
set -euo pipefail
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log "**** MAGMA gene-property ****"; echo "Job: ${SLURM_JOBID:-N/A}"
module purge; module load anaconda3/2024.10-1
mkdir -p logs

MAGMA=/ocean/projects/bio250020p/shared/opt/magma-v1.10/magma
REF=/ocean/projects/bio250020p/shared/opt/magma-v1.10/g1000_eur
GENELOC=/ocean/projects/bio250020p/shared/opt/magma-v1.10/NCBI37.3.gene.loc
SPEC=./specificity/specificity_all.tsv
STATE_SPEC=./specificity/specificity_pericyte_state.tsv

mkdir -p magma/anno magma/genes magma/geneprop

## (A) SNP location from the reference bim (SNP CHR BP), then annotate to genes
if [[ ! -f magma/anno/g1000.genes.annot ]]; then
    log "Annotating SNPs to genes (window 10kb up/down)"
    awk '{print $2"\t"$1"\t"$4}' "${REF}.bim" > magma/anno/g1000.snploc
    "$MAGMA" --annotate window=10,10 \
        --snp-loc magma/anno/g1000.snploc --gene-loc "$GENELOC" \
        --out magma/anno/g1000
fi

## Per-column specificity files (marginal models) — split once
COLS=$(head -1 "$SPEC" | tr '\t' '\n' | tail -n +2)
declare -a COLIDX
i=1
for c in $COLS; do i=$((i+1)); COLIDX+=("$i:$c"); done

run_trait() {
    local trait="$1"
    [[ -f "magma/${trait}.pval" ]] || { log "  no pval for ${trait}; skip"; return; }
    # (B) gene analysis
    if [[ ! -f "magma/genes/${trait}.genes.raw" ]]; then
        log "  gene analysis: ${trait}"
        "$MAGMA" --bfile "$REF" \
            --pval "magma/${trait}.pval" use=SNP,P ncol=N \
            --gene-annot magma/anno/g1000.genes.annot \
            --out "magma/genes/${trait}"
    fi
    # (C marginal) one specificity column per model
    for ci in "${COLIDX[@]}"; do
        local idx="${ci%%:*}" name="${ci#*:}"
        local safe; safe=$(echo "$name" | tr ' /' '__')
        local cf="magma/geneprop/_covar.${safe}.tsv"
        cut -f1,"${idx}" "$SPEC" > "$cf"
        "$MAGMA" --gene-results "magma/genes/${trait}.genes.raw" \
            --gene-covar "$cf" \
            --out "magma/geneprop/${trait}.${safe}" >/dev/null 2>&1 || \
            log "    gene-prop failed: ${trait} x ${name}"
    done
    # (C joint) all pericyte states together (conditional enrichment)
    "$MAGMA" --gene-results "magma/genes/${trait}.genes.raw" \
        --gene-covar "$STATE_SPEC" \
        --out "magma/geneprop/${trait}.JOINT_states" >/dev/null 2>&1 || \
        log "    joint gene-prop failed: ${trait}"
    log "  done ${trait}"
}

for pv in magma/*.pval; do
    [[ -e "$pv" ]] || { log "No pval files — run 01.munge_gwas.py first"; break; }
    run_trait "$(basename "$pv" .pval)"
done

log "**** MAGMA done ****"
