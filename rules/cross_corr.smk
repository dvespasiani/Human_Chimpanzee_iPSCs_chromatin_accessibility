##================================================
## 7. Compute ENCODE cross correlation
## NB: use only 1 read even if data are PE 
## trim read to 50 bp using ENCODE script 
## align again but dont filter aligned reads
## convert bam to tagAlign file then
## subsample 25M reads from this file and finally
## compute the cross correlation score
##================================================
rule cc_read1:
  input:
    "../data/input_files/samples/{sample}_1.fastq.gz"
  output:
    temp("output/ENCODE_CC/Files/{sample}_1_trimmed.fastq.gz")
  log:
   "logs/ENCODE_CC/{sample}_fastq.log"
  group:
    "CrossCorrelation"
  run:
    shell("python ../common_scripts/trimfastq.py {input} 50 \
   | gzip -nc > {output} 2> {log} ")

## Align trimmed read to the genome
rule cc_align_read1:
  input:
    rules.cc_read1.output
  output:
    temp("output/ENCODE_CC/Files/{sample}_1_trimmed.bam")
  params:
    index=config['index_genome_dir'],
    threads=8
  group:
    "CrossCorrelation"
  log:
   "logs/ENCODE_CC/{sample}_cc_align_read1.log"
  shell:
    "bowtie2 -x {params.index} \
    --threads {params.threads} \
    -U {input} | \
    samtools view -Su  \
    | samtools sort -o {output} - 2> {log}"

## Filter alignment but don't dedup
rule cc_filter_alignment:
  input:
    rules.cc_align_read1.output
  output:
    temp("output/ENCODE_CC/Files/{sample}_1_trimmed_q30.bam")
  log:
   "logs/ENCODE_CC/{sample}_filter_alignment.log"
  group:
    "CrossCorrelation"
  shell:
    "samtools view -F 1804 -q 30 -b {input} -o {output} 2> {log}"

rule cc_bamt2tagAling:
  input:
    rules.cc_filter_alignment.output
  output:
    temp("output/ENCODE_CC/Files/{sample}_R1_trimmed_q30_SE.tagAlign.gz")
  log:
   "logs/ENCODE_CC/{sample}_bamt2tagAling.log"
  group:
    "CrossCorrelation"
  shell:
    """
    bedtools bamtobed -i {input} | \
    awk 'BEGIN{{OFS="\\t"}}{{$4='N';$5='1000';print $0}}' |\
    gzip -nc > {output} 2> {log}
    """
## Subsample 25 million reads for cross-correlation analysis
## Estimate read length from first 100 reads
rule subsample_aligned_reads:
  input:
   rules.cc_bamt2tagAling.output
  output:
   "output/ENCODE_CC/Files/{sample}.filt.sample.25Mreads.SE.tagAlign.gz"
  log:
   "logs/ENCODE_CC/{sample}_subsample_aligned_reads.log"
  params:
    nreads= 25000000
  group:
    "CrossCorrelation"
  run:
   shell("""
   zcat {input} | grep -v “chrM” | \
   shuf -n {params.nreads} --random-source=<(openssl enc -aes-256-ctr -pass pass:$(zcat -f {input} | \
    wc -c) -nosalt </dev/zero 2>/dev/null) | \
   gzip -nc > {output} 2> {log}
   """) 

rule cross_correlation_SSP:
  input:
    rules.subsample_aligned_reads.output
  output:
    CC_SCORES_FILE="output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.qc",
    CC_PLOT_FILE="output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.plot.pdf"
  group:
    "CrossCorrelation"
  log:
    "logs/ENCODE_CC/{sample}_filt_25Mreads.SE.spp.log"
  shell:
    "Rscript --max-ppsize=500000 $(which run_spp.R) \
    -c={input} -filtchr=chrM \
    -savp={output.CC_PLOT_FILE} \
    -out={output.CC_SCORES_FILE} 2> {log}"

