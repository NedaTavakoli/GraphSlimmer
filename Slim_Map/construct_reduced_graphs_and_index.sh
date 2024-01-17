
# # On hive
# ssh ntavakoli6@login-hive-slurm.pace.gatech.edu
# salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=1500G -t120:00:00
# salloc -A hive-saluru8 -phive-himem -N1 --ntasks-per-node=1 --mem-per-cpu=3000G -t120:00:00


# project_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged
# scratch=/storage/hive/scratch/6/ntavakoli6  

# cd /storage/hive/project/cse-aluru/ntavakoli6/hged
# module load gurobi
# module load anaconda3
# module load boost
# module load cmake


# samtools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/samtools-1.12/samtools
# bcftools=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/bcftools-1.9/bcftools
# bgzip=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/bgzip
# tabix=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/htslib-1.12/tabix
# GraphAligner=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/GraphAligner/bin/GraphAligner
# k8=/storage/hive/project/cse-aluru/ntavakoli6/hged/software/k8/k8-linux


###################################################################################################
alpha=1000
delta=1
CHRID=chr22
VCF=data/chr22.vcf.gz
FA=Homo_sapiens.GRCh37.dna.chromosome.22.fa


# ***************************************************************************** DO ONLY ONCE  ************************
# Get list of all variant positions using bctfools 
$bcftools query -f '%POS\n' $VCF  > all_variant_positons_${CHRID}.txt

# Get the list of only snps and indels using bcftools
$bcftools query -i '(TYPE="snp" || TYPE="indel") && GT="alt"' -f '%POS\n'  $VCF  > snps_indels_variant_positions_${CHRID}.txt

# Get the positions that are not snps and indels (all other positions):
comm -13 <(sort snps_indels_variant_positions_${CHRID}.txt) <(sort all_variant_positions_${CHRID}.txt) > not_snps_indels_positions_${CHRID}.txt

# ***************************************************************************** Do it for each config ************************
# Get list of all variant positions using bctfools 

# For config ${alpha} and 1 
# Get the retained variant positions from the ILP_sol_${alpha}_${delta}
# check if it zero shows retained variant positions *******
awk 'NR==FNR {values[NR]=$1; next} values[FNR] == 0 {print $1}' ILP_sol_${alpha}_${delta}.txt snps_indels_variant_positions_${CHRID}.txt  > retained_variants_${alpha}_${delta}_${CHRID}.txt

 # get reduced vcf file using retained variants

 #substeps:
 # 1 
 #combine no_snps_indels_positions_${CHRID}.txt and retained_variants_${alpha}_${delta}_${CHRID}.txt
cat retained_variants_${alpha}_${delta}_${CHRID}.txt not_snps_indels_positions_${CHRID}.txt | sort -n > retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.txt


# 2 
#I want to filter a  vcf file such that it only contains the variants in the file retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.txt using bcftools 

positions_of_interest=retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.txt
$bcftools view -h data/chr22.vcf > retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf
grep -Fwf ${positions_of_interest} data/chr22.vcf  >> retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants_copy.vcf
$bgzip -c retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf > retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf.gz
$tabix -p vcf retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf.gz



# *****************************************************************************


# construct the reduced graph
my_new_vcf=retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants.vcf.gz
my_new_vg=retained_snp_indel_variants_${alpha}_${delta}_${CHRID}_plus_other_variants
/usr/bin/time -v ./vg construct -r $FA -v ${my_new_vcf} -a > ${my_new_vg}.vg

# Indexing all variants graph, xg and gbwt indexes
/usr/bin/time -v ./vg index -x ${my_new_vg}.xg -G ${my_new_vg}.gbwt -v ${my_new_vcf}  ${my_new_vg}.vg


#  Prunning all variants graph pruning, Pruning with Haplotypes
/usr/bin/time -v ./vg prune -u -g ${my_new_vg}.gbwt -m node_mapping ${my_new_vg}.vg > ${my_new_vg}.pruned.vg
/usr/bin/time -v ./vg index -g ${my_new_vg}.pruned.gcsa -f node_mapping ${my_new_vg}.pruned.vg   # gcsa index only can obtain if the graph is pruned

# Unfolding the paths with vg prune -u creates duplicated nodes, and if you
#  don't pass the "node mapping" to vg index -g, you get an index that maps to non-existent nodes.

