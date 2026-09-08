#!/bin/bash
#SBATCH --account=bio250020p
#SBATCH --partition=RM-shared
#SBATCH --job-name=mouse_census_dl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/mouse_census_download.log

## CELLxGENE Census download.
##
## CHANGED 2026-09-08. This ran as a bare login-node script, on the stated
## grounds that "CELLxGENE Census requires outbound internet, which PSC compute
## nodes lack". That is not true here -- compute nodes do have outbound access --
## and the login-node run was slow enough to degrade the shared node. It is now
## a normal batch job like every other step in this module.
##
## It also now WRITES A LOG (logs/mouse_census_download.log). There was none
## before, which is why the Census release behind the shipped h5ad was
## unrecoverable -- see the `--census-version` note in 00.download_mouse_lung.py.
## The release is pinned in the script default; do not pass a moving alias here.

log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Mouse lung Census download ****"
echo "Job id: ${SLURM_JOBID}"; echo "Hostname: ${HOSTNAME}"

source /etc/profile.d/modules.sh
module purge
module load anaconda3/2024.10-1
eval "$(conda shell.bash hook)"
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env

python ../_h/00.download_mouse_lung.py --outdir "${1:-./}"
if [ $? -ne 0 ]; then log_message "Error: Census download failed"; exit 1; fi

conda deactivate
log_message "**** Done ****"
