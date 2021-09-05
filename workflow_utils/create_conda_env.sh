#!/bin/bash

# list of modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load anaconda3/2020.11
# module load sra-toolkit/2.10.5-centos_linux64  #for downloading files, otherwise unnecessary
# module load gcc/8.3.0  openmpi/3.1.4
# module load meme/5.1.1-python-3.7.4
# m

## First create and then activate env
conda create --name atac
conda activate atac

## install mamba in this environment
conda install -c conda-forge mamba

## install all dependencies you need in this latter env
mamba install -c biocore -c conda-forge -c bioconda -c anaconda -n atac \
snakemake  macs2 deepTools fastqc \
bowtie2 trimmomatic samtools samstats picard \
python=3.7.4 openssl 

## export everthing into a yaml file
# conda env export > environment.yml


