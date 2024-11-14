#!/bin/bash

directory_input=""
directory_output=""
imputed_file_dir=""
phenotype=""
outcome=""

for i in {18..22 X}; do
    run_plink_wes="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep Initial_input_${phenotype}.phe --chr 1-22 X Y XY\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {3..17}; do
    run_plink_wes="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep Initial_input_${phenotype}.phe --chr 1-22 X Y XY\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {1..2 Y XY}; do
    run_plink_wes="plink2 --bgen ukb22828_c${i}_b0_v3.bgen ref-first --sample ukb22828_c${i}_b0_v3.sample\
      --no-pheno --keep Initial_input_${phenotype}.phe --chr 1-22 X Y XY\
      --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.bgen" \
     -iin="${imputed_file_dir}/ukb22828_c${i}_b0_v3.sample" \
     -iin="/${directory_input}/Initial_input/Initial_input_${phenotype}.phe" \
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

