#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ras_circuit_3
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/ras_circuit_step3.log

## Step 3 of ras_circuit (Figures 5A tables + 5C network). Needs steps 1 AND 2.
## Submit from ras_circuit/_m:
##     sbatch -D ras_circuit/_m ras_circuit/_h/step_3.sh

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"; echo "Node: ${SLURM_NODENAME}"
source /etc/profile.d/modules.sh
module purge; module load anaconda3/2024.10-1; module list
eval "$(conda shell.bash hook)"
mkdir -p logs stats_data
for f in ras_circuit_pseudobulk.tsv.gz receiver_programs.tsv at1r_response_metadata.tsv.gz; do
    [ -f "$f" ] || { log_message "Error: missing $f -- run steps 1 and 2 first"; exit 1; }
done
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/07.ras_network.R --outdir ./stats_data --nboot 1000 \
    || { log_message "Error: 07"; exit 1; }
conda deactivate
log_message "**** Job ends ****"
