#!/bin/bash
#SBATCH --partition=shared,bluejay
#SBATCH --job-name=phate_plotting
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --output=plotting.log
#SBATCH --time=01:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module load R
module list

echo "**** Run PHATE on AT2 -- plotting ****"

Rscript ../_h/02.plotting_subclusters.R

echo "**** Job ends ****"
date
