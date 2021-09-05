#!/bin/sh
## use this scrip to sample fastq files to test the snakemake pipeline

wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/"

helpFunction()
{
   echo ""
   echo "Usage: $0 -o outdir -i inputdir -n numbreads"
   echo -e "\t-o output dir for the subsampled files"
   echo -e "\t-i input dir relative to base dir containing the fastq files"
   echo -e "\t-n number reads to subsample"
   echo -e "\t base dir: $wd"
   exit 1 # Exit script after printing help
}

while getopts "o:i:n:" flag; do
    case "${flag}" in
        o) outdir=${OPTARG};;
        i) inputdir=${OPTARG};;
        n) numbreads=${OPTARG};;
        ?) helpFunction ;; # Print helpFunction 
    esac
done

## Begin script
for fa in $wd$inputdir/*.fastq.gz; do
  sample="$(echo $(basename $fa)| cut -f 1 -d '/')"
  echo "sampling $numbreads reads from $sample and saving it to $outdir/${sample%.*}.gz";
  gzip -cd "$fa" | head -n $numbreads  | gzip -nc > "$outdir/${sample%.*}.gz" 
done
