#!/bin/bash
## Run on a login/data node: CELLxGENE Census requires outbound internet, which
## PSC compute nodes lack. Memory is modest (capped cell count). NOT an sbatch job.
log_message() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_message "**** Mouse lung Census download ****"
module purge
module load anaconda3/2024.10-1
conda activate /ocean/projects/bio250020p/shared/opt/env/scRNA_env
python ../_h/00.download_mouse_lung.py --outdir "./"
conda deactivate
log_message "**** Done ****"
