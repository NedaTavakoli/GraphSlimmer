
project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
scratch=/storage/hive/scratch/6/ntavakoli6  

cd /storage/hive/project/cse-aluru/ntavakoli6/hged
module load gurobi
module load anaconda3
module load boost
module load cmake

samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux

CHRID=chr22
VCF=data/chr22.vcf.gz
FA=Homo_sapiens.GRCh37.dna.chromosome.22.fa



######## Compute the percentage of INDLs and SNPs at different configs ##############################################################################################
cd /storage/hive/project/cse-aluru/ntavakoli6/hged/

$bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\t%TYPE\n' $VCF  > pos_type_snps_indels_variant_positions_${CHRID}.txt
awk '{count[$2]++} END {print "SNP Count:", count["SNP"], "(", (count["SNP"]/NR)*100, "%)"; print "INDEL Count:", count["INDEL"], "(", (count["INDEL"]/NR)*100, "%)"}' pos_type_snps_indels_variant_positions_${CHRID}.txt


CHRID=chr22
alpha=150
delta=2 
types_file="pos_type_snps_indels_variant_positions_${CHRID}.txt"

awk 'NR==FNR {values[NR]=$1; next} FNR==NR {next} values[FNR] == 1 {print values[FNR], $2}' ILP_sol_${alpha}_${delta}.txt $types_file > type_retained_variants_${alpha}_${delta}_${CHRID}.txt
# awk 'NR==FNR {values[NR]=$1; next} FNR==NR {next} values[FNR] == 1 {print values[FNR], $2}' ILP_sol_100_5_1283792 ../$types_file > ../type_retained_variants_${alpha}_${delta}_${CHRID}.txt
# awk 'NR==FNR {values[NR]=$1; next} FNR==NR {next} values[FNR] == 1 {print values[FNR], $2}' ILP_sol_75_1_1283792 ../$types_file > ../type_retained_variants_${alpha}_${delta}_${CHRID}.txt
# awk 'NR==FNR {values[NR]=$1; next} FNR==NR {next} values[FNR] == 1 {print values[FNR], $2}' ILP_sol_150_1_1463441 ../$types_file > ../type_retained_variants_${alpha}_${delta}_${CHRID}.txt
# awk 'NR==FNR {values[NR]=$1; next} FNR==NR {next} values[FNR] == 1 {print values[FNR], $2}' ILP_sol_150_2_1463441 ../$types_file > ../type_retained_variants_${alpha}_${delta}_${CHRID}.txt


# Calculate and print percentage of SNPs and INDELs
awk '{count[$2]++} END {print "SNP Count:", count["SNP"], "(", (count["SNP"]/NR)*100, "%)"; print "INDEL Count:", count["INDEL"], "(", (count["INDEL"]/NR)*100, "%)"}' type_retained_variants_${alpha}_${delta}_${CHRID}.txt
awk '{count[$2]++} END {print "SNP Count:", count["SNP"], "(", (count["SNP"]/NR)*100, "%)"; print "INDEL Count:", count["INDEL"], "(", (count["INDEL"]/NR)*100, "%)"}' ../type_retained_variants_${alpha}_${delta}_${CHRID}.txt


 ##################################################################################################################################################################

Output:

Totoal: Total numner of variannts in chr22: 1098701, SNP Count: 1056333 ( 96.1438 %), INDEL Count: 42368 ( 3.85619 %)
alpha, delta, reed SNP%, INDEL%
------------------------------
75, 1,  SNP Count: 953586 ( 98.554 %), INDEL Count: 13991 ( 1.44598 %), total_retained: 1098701 - (953586 + 13991) = 133124 ( 12.144 %), removed: 1098701 - 133124 = 965577 ( 87.856 %)
100, 1, SNP Count: 941270 ( 98.628 %), INDEL Count: 13094 ( 1.37201 %), total_retained: 1098701 - (941270 + 13094) = 146337 ( 13.322 %), removed: 1098701 - 146337 = 952364 ( 86.678 %)
150, 1 , SNP Count: 923292 ( 98.7184 %), INDEL Count: 11987 ( 1.28165 %), total_retained: 1098701 - (923292 + 11987) = 155422 ( 14.146 %), removed: 1098701 - 155422 = 943279 ( 85.854 %)
150, 2, SNP Count: 1017688 ( 97.9487 %), INDEL Count: 21313 ( 2.0513 %), total_retained: 1098701 - (1017688 + 21313) =  529700 ( 48.183 %, removed: 1098701 - 529700 = 568001 ( 51.817 %)  
1000, 1, SNP Count: 461500 ( 98.9681 %), INDEL Count: 4812 ( 1.03193 %), total_retained: 1098701 - (461500 + 4812) =  636389 ( 57.905 %), removed: 1098701 - 636389 = 462312 ( 42.095 %)





Report graph sizes: # edges and #nodes



Report mapping accuracy for 100, 1:  Complete_mapping_accuracy: %, reduced_mapping_accuracy: %
Report mapping accuracy for 25, 1:   Complete_mapping_accuracy: %, reduced_mapping_accuracy: %
Report mapping accuracy for 1000, 1: Complete_mapping_accuracy: 80.89%, reduced_mapping_accuracy: 79.66%






# ############################## OLD

# cd simulated_reads

# CHRID=chr22
# output_gam_file=simulated_reads_len_100_errRate_0.006_indelRate_0.004_all_samples.gam
# mapped_reads_with_compare_option=simulated_reads_len_100_errRate_0.006_indelRate_0.004_all_samples_with_compare_option.txt


# # Mapping the simulated reads to the original graph uisng vg
# /usr/bin/time -v .././vg map -t 24 -x ../${CHRID}.xg -g ../${CHRID}.pruned.gcsa --gbwt-name ../${CHRID}.gbwt -G ${output_gam_file} --compare -j > ${mapped_reads_with_compare_option}  # map with compare option 
# cat ${mapped_reads_with_compare_option} | jq .correct | sed s/null/0/ | awk '{i+=$1; n+=1} END {print i/n}'
# # 0.808862

# ##############################
# alpha=100
# delta=1
# CHRID=chr22
# bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
# VCF=retained_snp_indel_variants_${alpha}_${delta}_chr22_plus_other_variants.vg

# $bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\n'  $VCF  > snps_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt
# $bcftools query -i '(TYPE="indel") && GT="alt"' -f '%POS\n'  $VCF  > indels_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt

# Num_Snp=$(wc -l snps_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt | awk '{print $1}')
# Num_Indel=$(wc -l indels_variant_positions_${CHRID}_${alpha}_${delta}_retained.txt | awk '{print $1}')
# echo "Number of retained snps: ${Num_Snp}"
# echo "Number of retained indels: ${Num_Indel}"