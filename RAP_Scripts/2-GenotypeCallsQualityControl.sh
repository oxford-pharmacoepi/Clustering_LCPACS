#!/bin/sh
# This script runs the QC process using PLINK on the merged file generated in 
# partB-merge-files.sh as described in the 

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Genotype_calls_qc"
phenotype="LongCovid_cohort"

run_plink_qc="plink2 --bfile ukb22418_25_merged\
 --chr 1-22\
 --exclude range extended_LD_regions_hg37_GRCh37.txt\
 --maf 0.01\
 --keep ${phenotype}.phe\
 --mac 20\
 --geno 0.1\
 --hwe 1e-6\
 --max-alleles 2\
 --mind 0.1\
 --write-snplist --write-samples\
 --no-id-header --out  snps_qc_pass_${phenotype}"

dx run swiss-army-knife -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.bed"\
   -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.bim"\
   -iin="/${directory_input}/Genotype_calls_merged/ukb22418_25_merged.fam"\
   -iin="/${directory_input}/Initial_input/extended_LD_regions_hg37_GRCh37.txt"\
   -icmd="${run_plink_qc}" --tag="Step1_Regenie" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory_output}/" --brief --yes --name="2-GenotypeCallsQualityControl"