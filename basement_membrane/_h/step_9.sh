#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-small
#SBATCH --job-name=bm_smad_repl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
## Independent evaluation of the SMAD -> (BM - fibrillar) association in
## GSE136831 (defect P1-20a). Reads the pseudobulk step_3.sh already built, so
## the 8.4 GB h5ad is NOT re-read; this is a minutes-long job.
#SBATCH --time=00:20:00
#SBATCH --output=logs/bm_smad_replication.log

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Job starts ****"; echo "Job id: ${SLURM_JOBID}"

## The batch shell does not source the login profile, so `module` and `conda` must
## be bootstrapped here rather than inherited from the submitting shell.
source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
module list
eval "$(conda shell.bash hook)"

conda activate /ocean/projects/bio250020p/shared/opt/env/R_env

log_message "**** SMAD arm replication in GSE136831 pericytes ****"
Rscript ../_h/18.smad_replication_gse136831.R \
        --pseudobulk ./gse136831_bm_pseudobulk.tsv.gz \
        --panels ./bm_panel_genes.tsv \
        --outdir ./stats_data \
        --min-cells 5 \
        --compartment Pericyte
if [ $? -ne 0 ]; then log_message "Error: Rscript execution failed"; exit 1; fi

conda deactivate
log_message "**** Job ends ****"
