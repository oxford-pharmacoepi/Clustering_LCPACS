#!/bin/bash

directory_input="GWAS/Long_covid_gilead"
directory_output="GWAS/Long_covid_gilead/Intermediary_files"
imputed_file_dir="GWAS/Plink_files"
phenotype="LongCovid_cohort"
outcome="state"

for i in {18..22} XY; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --maf 0.001 --mac 10 --merge-x --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim" \
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

for i in {3..17}; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --maf 0.001 --mac 10 --merge-x --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim" \
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x36"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done


for i in {1..2} X; do
    run_plink_wes="plink2 --bfile plink_files_c${i}\
      --no-pheno --keep ${phenotype}.phe --chr 1-23 25\
      --maf 0.001 --mac 10 --merge-x --geno 0.1 --hwe 1e-15 --mind 0.1\
      --write-snplist --write-samples --no-id-header\
      --rm-dup force-first --out c${i}_stepE_${phenotype}"

    dx run swiss-army-knife -iin="${imputed_file_dir}/plink_files_c${i}.bim" \
     -iin="${imputed_file_dir}/plink_files_c${i}.bed" \
     -iin="${imputed_file_dir}/plink_files_c${i}.fam" \
     -iin="${directory_input}/Initial_input/LongCovid_cohort.phe"\
     -icmd="${run_plink_wes}" --tag="Step2" --instance-type "mem1_ssd1_v2_x72"\
     --name "StepE_chr${i}_${phenotype}"\
     --destination="${project}:/${directory_output}/" --brief --yes
done

