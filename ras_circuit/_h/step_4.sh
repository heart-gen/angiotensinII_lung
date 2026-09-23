#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=ras_circuit_4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --output=logs/ras_circuit_step4.log

## Step 4 of ras_circuit (Figure 5F: outgoing pericyte communication).
## Independent of steps 2-3. Submit from ras_circuit/_m:
##     sbatch -D ras_circuit/_m ras_circuit/_h/step_4.sh
##   08  re-filter existing all-pairs LIANA tables (no LIANA re-run)
##   09  NicheNet with Pericytes as SENDER + permutation null + validation design
##   --  pseudobulk of the validation genes over the CCC niche
##   10  donor-level validation against detection-matched null composites

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"; echo "Node: ${SLURM_NODENAME}"
source /etc/profile.d/modules.sh
module purge; module load anaconda3/2024.10-1; module list
eval "$(conda shell.bash hook)"
mkdir -p logs stats_data nichenet_outgoing

conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/receiver_programs.py --out ./receiver_programs.tsv \
       --targets-out ./outgoing_targets.tsv || { log_message "Error: programmes"; exit 1; }
log_message "**** 08 outgoing LIANA ****"
python ../_h/08.liana_outgoing.py --liana-dir ../../cell_communication/_m \
       --outdir ./ || { log_message "Error: 08"; exit 1; }
conda deactivate

log_message "**** 09 NicheNet, pericyte as sender ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/09.nichenet_outgoing.R --outdir ./nichenet_outgoing --n-perm 1000 \
    || { log_message "Error: 09"; exit 1; }
conda deactivate

log_message "**** pseudobulk of validation genes ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../../basement_membrane/_h/02.niche_pseudobulk.py \
       --adata ../../cell_communication/_m/ccc_niche.h5ad \
       --genes ./outgoing_pb_genes.tsv \
       --outfile ./outgoing_pseudobulk.tsv.gz || { log_message "Error: pseudobulk"; exit 1; }
conda deactivate

log_message "**** 10 donor-level validation ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/R_env
Rscript ../_h/10.outgoing_donor_validation.R --outdir ./stats_data --n-null 1000 \
    || { log_message "Error: 10"; exit 1; }

## Last, because it reads the finished output of every other step: the evidence layer
## that annotates the Figure 5A schematic.
log_message "**** 11 DAG annotations ****"
Rscript ../_h/11.dag_annotations.R --outdir ./stats_data \
    || { log_message "Error: 11"; exit 1; }
conda deactivate
log_message "**** Job ends ****"
