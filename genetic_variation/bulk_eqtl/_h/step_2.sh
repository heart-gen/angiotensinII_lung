#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=combine_parquet
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --output=combine.log
#SBATCH --time=05:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility

module load htslib
module list

## Edit with your job command

echo "**** Combine parquet files ****"

python ../_h/02.combine_parquet.py
bgzip -f GTEx_v8.nominal.txt
tabix -f GTEx_v8.nominal.txt.gz

echo "**** Job ends ****"
date
