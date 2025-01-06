#!/bin/bash

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediate_files"
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
phenotype="LongCovid_cohort"
outcome="state"

for chr in {18..22} X XY; do
run_regenie_cmd="regenie --step 2 --bgen ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen\
    --sample ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_ivqc_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol age2\
    --covarCol agesex\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen"\
-iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample"\
-iin="/${directory_output}/c${chr}_ivqc_${phenotype}.snplist"\
-iin="/${directory_input}/Initial_input/${phenotype}.phe"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_pred.list"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_1.loco.gz"\
-icmd="${run_regenie_cmd}" --tag="RegenieStep2" --instance-type "mem1_ssd1_v2_x16"\
--name="RegenieStep2_chr${chr}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done

for chr in {3..17}; do
run_regenie_cmd="regenie --step 2 --bgen ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen\
    --sample ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_ivqc_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol age2\
    --covarCol agesex\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen"\
-iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample"\
-iin="/${directory_output}/c${chr}_ivqc_${phenotype}.snplist"\
-iin="/${directory_input}/Initial_input/${phenotype}.phe"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_pred.list"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_1.loco.gz"\
-icmd="${run_regenie_cmd}" --tag="RegenieStep2" --instance-type "mem1_ssd1_v2_x36"\
--name="RegenieStep2_chr${chr}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done

for chr in {1..2} X; do
run_regenie_cmd="regenie --step 2 --bgen ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen\
    --sample ${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample\
    --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_ivqc_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol age2\
    --covarCol agesex\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.bgen"\
-iin="${imputed_file_dir}/ukb22828_c${chr}_b0_v3.sample"\
-iin="/${directory_output}/c${chr}_ivqc_${phenotype}.snplist"\
-iin="/${directory_input}/Initial_input/${phenotype}.phe"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_pred.list"\
-iin="/${directory_output}/RegenieStep2-${phenotype}_results_1.loco.gz"\
-icmd="${run_regenie_cmd}" --tag="RegenieStep2" --instance-type "mem1_ssd1_v2_x72"\
--name="RegenieStep2_chr${chr}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done