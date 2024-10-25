#!/bin/sh
# This script runs the QC process using PLINK on the merged file generated in 
# partB-merge-files.sh as described in the 

directory_input=""
directory_output=""
phenotype=""

run_plink_qc="plink2 --bfile ukb22418_25_merged\
 --keep Initial_input_${phenotype}.phe --chr 1-22 X Y XY\
 --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15\
 --mind 0.1 --write-snplist --write-samples\
 --no-id-header --out  snps_qc_pass_${phenotype}"

dx run swiss-army-knife -iin="/${directory_input}/Merged_files/ukb22418_25_merged.bed"\
   -iin="/${directory_input}/Merged_files/ukb22418_25_merged.bim"\
   -iin="/${directory_input}/Merged_files/ukb22418_25_merged.fam"\
   -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe"\
   -icmd="${run_plink_qc}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="${project}:/${directory_output}/" --brief --yes --name="StepC_${phenotype}"