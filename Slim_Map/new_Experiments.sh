cd simulated_reads

CHRID=chr22
output_gam_file=simulated_reads_len_100_errRate_0.006_indelRate_0.004_all_samples.gam
mapped_reads_with_compare_option=simulated_reads_len_100_errRate_0.006_indelRate_0.004_all_samples_with_compare_option.txt


# Mapping the simulated reads to the original graph uisng vg
/usr/bin/time -v .././vg map -t 24 -x ../${CHRID}.xg -g ../${CHRID}.pruned.gcsa --gbwt-name ../${CHRID}.gbwt -G ${output_gam_file} --compare -j > ${mapped_reads_with_compare_option}  # map with compare option 
cat ${mapped_reads_with_compare_option} | jq .correct | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
# 0.808862

##############################
alpha=100
delta=1
CHRID=chr22
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
VCF=retained_snp_indel_variants_${alpha}_${delta}_chr22_plus_other_variants.vg

$bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\n'  $VCF  > snps_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt
$bcftools query -i '(TYPE="indel") && GT="alt"' -f '%POS\n'  $VCF  > indels_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt

Num_Snp=$(wc -l snps_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt | awk '{print $1}')
Num_Indel=$(wc -l indels_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt | awk '{print $1}')
echo "Number of retained snps: ${Num_Snp}"
echo "Number of retained indels: ${Num_Indel}"