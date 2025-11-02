#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=50G
#$ -N eQTL_download
#$ -o ./summary.log
#$ -e ./summary.log
#$ -m e

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
echo "**** Download files ****"
gsutil -m cp \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr1.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr10.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr11.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr12.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr13.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr14.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr15.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr16.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr17.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr18.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr19.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr2.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr20.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr21.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr22.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr3.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr4.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr5.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr6.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr7.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr8.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chr9.parquet" \
  "gs://gtex-resources/GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_EUR_eQTL_all_associations/Lung.v8.EUR.allpairs.chrX.parquet" \
  .

echo "**** Job ends ****"
date
