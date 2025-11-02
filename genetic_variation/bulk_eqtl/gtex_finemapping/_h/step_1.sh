#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=download_finemap
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --output=summary.log
#SBATCH --time=06:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility

module list

## Edit with your job command

echo "**** Download fine mapping ****"
wget https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/fine-mapping-cis-eqtl/GTEx_v8_finemapping_DAPG.tar
wget https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/fine-mapping-cis-eqtl/GTEx_v8_finemapping_CAVIAR.tar
wget https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/fine-mapping-cis-eqtl/GTEx_v8_finemapping_CaVEMaN.tar

echo "**** Untar ****"
tar xvf GTEx_v8_finemapping_DAPG.tar
tar xvf GTEx_v8_finemapping_CAVIAR.tar
tar xvf GTEx_v8_finemapping_CaVEMaN.tar

echo "**** Remove old files ****"
rm *tar

echo "**** Job ends ****"
date
