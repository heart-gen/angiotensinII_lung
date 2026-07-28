#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=peri_annot_support
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=logs/annotation_support.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading mamba environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

## Reads the h5ad that step_0.sh already wrote and adds the marker dot-plot table
## next to it; it does NOT recluster. The companion per-cluster bootstrap Jaccard
## table that Figure S3 also needs comes from step_0c.sh
## (00.state_discovery.py --stability-only), not from here.
log_message "**** Marker dot-plot table (Figure S3 panel C) ****"
python ../_h/00b.annotation_support.py \
       --adata ./pericyte_states.h5ad \
       --outdir "./" \
       --discovery ../_h/00.state_discovery.py

if [ $? -ne 0 ]; then
    log_message "Error: Python execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
