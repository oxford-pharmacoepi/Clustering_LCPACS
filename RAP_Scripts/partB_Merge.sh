#!/bin/sh

# How to Run:
# ./partB_Merge.sh on the command line

# What this file does:
# Merge genotype calls

directory_input="GWAS/Long_covid_gilead/Initial_input"
directory_output="GWAS/Long_covid_gilead/Intermediary_files"
phenotype="LongCovid_cohort"

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9XY]* . ;\
        ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt;\
        plink --merge-list files_to_merge.txt --make-bed\
        --chr 1-22 X Y XY --out ukb22418_25_merged;\
        rm files_to_merge.txt;"

dx run swiss-army-knife -iin="/${directory_input}/${phenotype}.phe"\
   -icmd="${run_merge}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/${directory_output}" --brief --yes --name="partB_${phenotype}"
