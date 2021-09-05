#!/bin/bash

## Script used to create the conda env necessary to run workflow
conda create --name mamba
conda activate mamba

## install mamba in this environment
conda install -c conda-forge mamba

## create an empty env 
conda create --name atac

## install all dependencies you need in this latter env
mamba install -c biocore -c conda-forge -c bioconda -n atac \
snakemake jinja2 networkx graphviz macs2 deepTools fastqc \
bowtie2 ucsc-liftover trimmomatic sra-tools samtools samstats \
multiqc picard bedtools  
coreutils matplotlib 

## export everthing into a yaml file
conda env export > env/atac.yaml

## then activate the env
conda activate atac

# Download files (if they dont exist)
# do it for humans and chimp separately

function check_fastq_files () {
  dir=$1
  expected_no_fastqc=$2
  
  count_files=`ls -1 $dir/*.fastq.gz 2>/dev/null | wc -l`

if [ $count_files == expected_no_fastqc ] 
then
  echo "All files are here"
elif [ $count_files != expected_no_fastqc ]
then
  echo "Downloading files"
  `fastq-dump --split-files --skip-technical --gzip $(<$dir/SraAccList.txt) --outdir $dir `
fi
  
}

echo "Cheking presence human files..."
check_fastq_files "hg38/data/samples" 12

echo "Cheking presence chimp files..."
check_fastq_files "panTro5/data/samples" 10
 
#ENCODE blacklist
echo "Getting hg38 ENCODE blacklisted regions"
encode_dir="data/ENCODE_blacklisted"

if [ ! -d $encode_dir ]; then
  mkdir -p $encode_dir
else
  echo "Directory exists, downloading files"
fi

cd $encode_dir/
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip ./hg38.blacklist.bed.gz

## liftOver
echo "Liftover to create panTro5 blacklisted regions"

liftOver hg38.blacklist.bed  \
../LiftOver_chains/hg38ToPanTro5.over.chain  \
panTro5_blacklist_v2.bed   panTro5_unmapped.beds

echo 'Now you are all set, enjoy' 
