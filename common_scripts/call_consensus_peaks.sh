## script used to merge all Tn5 shifted bam files and call consensus peaks on these

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input -s suffix -o output -g genome_size -b blacklist -p prefix"
   echo -e "\t-i input dir containing bams"
   echo -e "\t-s suffix input bams"
   echo -e "\t-o output dir"
   echo -e "\t-g genome size"
   echo -e "\t-b encode blacklist"
   echo -e "\t-p prefix output merged files"
   exit 1 # Exit script after printing help
}

while getopts ":i:s:o:g:b:p:" flag; do
    case "${flag}" in
        i) input=${OPTARG};;
        s) suffix=${OPTARG};;
        o) output=${OPTARG};;
        g) genome_size=${OPTARG};;
        b) blacklist=${OPTARG};;
        p) prefix=${OPTARG};;
        ?) helpFunction ;; # Print helpFunction 
    esac
done

module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load bedtools/2.27
module load samtools/1.9
module load ucsc/21072020
module load macs2/2.2.7.1-python-3.7.4

input_bams=( "$input"/*"$suffix" )

echo "merging bam files";
samtools merge ${input}/${prefix}_tn5_shifted.bam ${input_bams[@]}

echo "sorting ${prefix}_tn5_shifted.bam file";
samtools sort ${input}/${prefix}_tn5_shifted.bam -O BAM -o ${input}/${prefix}_tn5_shifted_sorted.bam

echo "indexing ${prefix}_tn5_shifted.bam file";
samtools index -c ${input}/${prefix}_tn5_shifted_sorted.bam ${input}/${prefix}_tn5_shifted_sorted.bam.bai

echo "Calling consensus peaks" 
macs2 callpeak --format BAMPE --treatment ${input}/${prefix}_tn5_shifted_sorted.bam \
    --keep-dup all --outdir ${output} \
    --name ${prefix}_macs2_default --shift 100 --nomodel -B --SPMR \
    --extsize 300 --pvalue 0.01 --call-summits -g ${genome_size}

# echo "Removing peaks within blacklisted regions";
# bedtools intersect -v  \
# -a ${output}/${prefix}_macs2_default_peaks.narrowPeak \
# -b $blacklist | awk 'BEGIN{{OFS="\t"}} {{if ($5>1000) $5=1000; print $0}}' \
#  | awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y|2A|2B)$/'> ${output}/${prefix}_macs2_default_filtered.narrowPeak 

# echo "Sorting final consensus peak set";
# sort -k 8gr,8gr ${output}/${prefix}_macs2_default_filtered.narrowPeak | awk 'BEGIN{{OFS="\t"}}{{$4="Peak_"NR ; print $0}}' \
# |  head -n 300000 | gzip -nc > ${output}/${prefix}_macs2_default_peaks_filtered_sorted.narrowPeak.gz

