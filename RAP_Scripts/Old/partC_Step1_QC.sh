#!/bin/sh
# This script runs the QC process using PLINK on the merged file generated in 
# partB-merge-files.sh as described in the 

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediary_files"
phenotype="LongCovid_cohort"

run_plink_qc="plink2 --bfile ukb22418_25_merged\
 --keep ${phenotype}.phe --chr 1-23 25\
 --maf 0.001 --mac 10 --geno 0.1 --hwe 1e-15\
 --mind 0.1 --merge-x --write-snplist --write-samples\
 --no-id-header --out  snps_qc_pass_${phenotype}"

dx run swiss-army-knife -iin="/${directory_input}/Intermediary_files/ukb22418_25_merged.bed"\
   -iin="/${directory_input}/Intermediary_files/ukb22418_25_merged.bim"\
   -iin="/${directory_input}/Intermediary_files/ukb22418_25_merged.fam"\
   -iin="/${directory_input}/Initial_input/${phenotype}.phe"\
   -icmd="${run_plink_qc}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory_output}/" --brief --yes --name="StepC_${phenotype}"