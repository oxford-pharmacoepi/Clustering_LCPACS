#!/bin/bash

#!/bin/bash


directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediary_files"
imputed_file_dir="GWAS/Plink_files"
phenotype="LongCovid_cohort"
outcome="state"

for chr in {1..22} X XY; do
  run_regenie_cmd="regenie --step 2 --bed plink_files_c${chr} --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/${phenotype}.phe"\
   -iin="/${directory_output}/StepD-${phenotype}_results_pred.list"\
   -iin="/${directory_output}/StepD-${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done

for chr in {3..13}; do
  run_regenie_cmd="regenie --step 2 --bed plink_files_c${chr} --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/${phenotype}.phe"\
   -iin="/${directory_output}/StepD-${phenotype}_results_pred.list"\
   -iin="/${directory_output}/StepD-${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done


for chr in {1..2} X; do
  run_regenie_cmd="regenie --step 2 --bed plink_files_c${chr} --out ${phenotype}_assoc.c${chr}\
    --phenoFile ${phenotype}.phe --covarFile ${phenotype}.phe\
    --bt --approx --firth-se --firth --extract c${chr}_stepE_${phenotype}.snplist\
    --phenoCol ${outcome}\
    --covarCol sex\
    --covarCol age\
    --covarCol batch\
    --covarCol pc{1:10}\
    --pred StepD-${phenotype}_results_pred.list --bsize 200\
    --pThresh 0.05 --minMAC 3 --threads 16 --gz"
    
  dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${chr}.bed"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.bim"\
   -iin="${imputed_file_dir}/plink_files_c${chr}.fam"\
   -iin="/${directory_output}/c${chr}_stepE_${phenotype}.snplist"\
   -iin="/${directory_input}/Initial_input/${phenotype}.phe"\
   -iin="/${directory_output}/StepD-${phenotype}_results_pred.list"\
   -iin="/${directory_output}/StepD-${phenotype}_results_1.loco.gz"\
   -icmd="${run_regenie_cmd}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
   --name="StepF_chr${chr}_${phenotype}"\
   --destination="${project}:/${directory_output}/" --brief --yes
done