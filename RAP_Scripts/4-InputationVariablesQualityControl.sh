#!/bin/bash

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediate_files"
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
phenotype="LongCovid_cohort"
outcome="state"

for i in {18..22} XY; do
run_plink_wes="plink2 --bgen ${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen ref-first\
      --sample ${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --merge-x \
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_ivqc_${phenotype}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
-iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
-iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
-icmd="${run_plink_wes}" --tag="Step2_Regenie" --instance-type "mem1_ssd1_v2_x16"\
--name "4-InputationVariablesQualityControl_chr${i}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done

for i in {3..17}; do
run_plink_wes="plink2 --bgen ${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen ref-first\
      --sample ${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --merge-x \
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_ivqc_${phenotype}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
-iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
-iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
-icmd="${run_plink_wes}" --tag="Step2_Regenie" --instance-type "mem1_ssd1_v2_x36"\
--name "4-InputationVariablesQualityControl_chr${i}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done

for i in {1..2} X; do
run_plink_wes="plink2 --bgen ${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen ref-first\
      --sample ${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --merge-x \
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_ivqc_${phenotype}"

dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
-iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
-iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
-icmd="${run_plink_wes}" --tag="Step2_Regenie" --instance-type "mem1_ssd1_v2_x72"\
--name "4-InputationVariablesQualityControl_chr${i}_${phenotype}"\
--destination="${project}:/${directory_output}/" --brief --yes
done
