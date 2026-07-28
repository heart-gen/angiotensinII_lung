#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=state_jaccard
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/state_bootstrap_jaccard.log

## Per-cluster bootstrap Jaccard for the ALREADY-PUBLISHED clustering solution
## (supplementary Table S4). 00.state_discovery.py computes per-cluster Jaccard
## internally but only ever persisted the median across clusters, so the table had
## no way to say which cluster was least reproducible.
##
## This runs --stability-only, which reads the OUTPUT pericyte_states.h5ad, re-runs
## only the bootstrap at the stored chosen solution, and verifies the recomputed
## labels reproduce `pericyte_state` exactly before writing
## stability/cluster_bootstrap_jaccard.tsv. It does NOT rewrite the h5ad, the
## annotations, or the stability grid -- a full re-run would regenerate the
## published object and could perturb labels every downstream module depends on.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

log_message "**** Per-cluster bootstrap Jaccard (stability-only) ****"
python ../_h/00.state_discovery.py \
       --adata ./pericyte_states.h5ad \
       --outdir "./" \
       --stability-only

if [ $? -ne 0 ]; then
    log_message "Error: stability-only run failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
