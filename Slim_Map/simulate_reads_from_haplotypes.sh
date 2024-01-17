
CHRID=chr22

# For short reads
LEN=100
EROR=0.01
INDEL_error=0.0001
COUNT=40 # Total number of reads you want to simulate

For long reads
LEN=1000
EROR=0.025
INDEL_error=0.1
COUNT=40 # Total number of reads you want to simulate
CHRID=chr22

simulated_dir=/storage/hive/project/cse-aluru/ntavakoli6/hged/simulated_reads
output_gam_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.gam
output_fastq_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.fasta
output_sam_file=${simulated_dir}/simulated_reads_len_${LEN}_errRate_${EROR}_indelRate_${INDEL_error}_all_samples.sam

##################################################################################################
# cd /storage/hive/project/cse-aluru/ntavakoli6/hged


# Simulate reads using vg sim and append to the output GAM file
sample_list_file="/storage/hive/project/cse-aluru/ntavakoli6/hged/samples.txt"  # Path to your sample names file

sample_names=()
while IFS= read -r line; do
    sample_names+=("$line")
done < "$sample_list_file"



# Randomly shuffle the sample names array
shuf_sample_names=($(shuf -e "${sample_names[@]}"))


while [[ ${#sample_names[@]} -gt 0 ]]; do
    random_index=$((RANDOM % ${#sample_names[@]}))
    random_sample_name="${sample_names[$random_index]}"
    
    ./vg sim -x ${CHRID}.xg -g ${CHRID}.gbwt -m ${random_sample_name} -e ${EROR} -i ${INDEL_error} -n ${COUNT} -l ${LEN} -a  >> ${output_gam_file}

    # Process the randomly selected sample name here
    
    # Remove the selected sample name from the array
    unset "sample_names[$random_index]"
    sample_names=("${sample_names[@]}")
done


/usr/bin/time -v ./vg view -X ${output_gam_file} > ${output_fastq_file}
/usr/bin/time -v ./vg surject -s ${output_gam_file} -x ${CHRID}.xg > ${output_sam_file} #truth on linear genome






