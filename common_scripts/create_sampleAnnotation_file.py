
## Script used to create the annotation file for human samples to be used by ChrAccR
## launch this script from snakemake so that bamfilesFull reports only the relative path rather than the full path
import os 
from os.path import isfile, join
import re
import csv

## set working directory
# os.chdir(os.getcwd()+'/hg38/')
print('Creating an annotation file from Tn5 shifted BAMs')

Tn5_bam_dir = 'output/Post_alignment/Files/combined/'
outdir = 'output/PeakCalling/QC/'

## get sample IDs
bamfiles = []
for file in os.listdir(Tn5_bam_dir):
    if file.endswith(".bam"):
        bamfiles.append(file)

sampleId=[re.sub('.bam', '', f) for f in bamfiles]

bamfilesFull=[Tn5_bam_dir + f for f in bamfiles]

sample_info=[sampleId,bamfiles,bamfilesFull]

## write annotation file
annotation_file=open(outdir+'SampleAnnotation.txt','w')

annotation_file.write("sampleId\tbamfiles\tbamfilesFull\n")

for row in zip(*sample_info):
    annotation_file.write('\t'.join(row)+'\n')

annotation_file.close()

