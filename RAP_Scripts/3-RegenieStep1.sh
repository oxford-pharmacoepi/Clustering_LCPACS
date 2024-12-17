#!/bin/sh

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediate_files"
phenotype="LongCovid_cohort"
outcome="state"

run_regenie_step1="regenie --step 1\
 --lowmem --out regenie_step_1_${phenotype}_results --bed ukb22418_25_merged\
 --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
 --extract snps_qc_pass_${phenotype}.snplist --phenoCol ${outcome}\
 --covarCol sex\
 --covarCol age\
 --covarCol age2\
 --covarCol agesex\
 --covarCol batch\
 --covarCol pc{1:10}\
 --bsize 1000 --bt --loocv --gz --threads 16"

dx run swiss-army-knife -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.bed" \
   -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.bim" \
   -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.fam"\
   -iin="/${directory_input}/Intermediate_files/snps_qc_pass_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/${phenotype}.phe" \
   --name="Regenie_step_1_"${phenotype}\
   -icmd="${run_regenie_step1}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory_output}" --brief --yes
