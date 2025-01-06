#!/bin/sh

merge_cmd='out_file="LongCovid_assoc.regenie.merged.txt"

# Use dxFUSE to copy the regenie files into the container storage
cp /mnt/project/GWAS/Long_covid_gilead/Intermediate_files/LongCovid_assoc.c*_state.regenie.gz .
gunzip LongCovid_assoc.c*_state.regenie.gz

# add the header back to the top of the merged file
echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file

files="./LongCovid_assoc.c*_state.regenie"
for f in $files
do
   tail -n+2 $f | tr " " "\t" >> $out_file
done

# remove regenie files
rm *LongCovid_assoc.c*_state.regenie'

directory="GWAS/Long_covid_gilead"
phenotype="LongCovid_cohort"
outcome="state"

dx run swiss-army-knife -iin="/${directory}/Intermediate_files/${phenotype}_assoc.c1_${outcome}.regenie.gz" \
-icmd="${merge_cmd}" --tag="FinalMerge" --instance-type "mem1_ssd1_v2_x16"\
--destination="${project}:/${directory}/" --brief --yes\
--name="FinalMerge_${phenotype}"