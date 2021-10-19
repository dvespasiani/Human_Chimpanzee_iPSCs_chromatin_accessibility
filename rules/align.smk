## =================================================
##  1. Align to reference genome using Bowtie2
##  then sort bam file by genomic coordinates
## =================================================
rule align:
  input:
    r1 = reads_input_dir + "{sample}_1.fastq.gz",
    r2 = reads_input_dir + "{sample}_2.fastq.gz"
  output:
    aligment_outdir + "{sample}.bam"
  params:
    index = index_genome_dir
  group:
    main_group
  log:
    aligment_logdir + "{sample}_bowtie_Align_qc.log"
  shell:
    "bowtie2 -q -X 2000 -x {params.index} \
    -1 {input.r1} -2 {input.r2} 2> {log} \
    | samtools view -Su | samtools sort -n -o {output}"

# ## Alignment statistics
# rule alignStat:
#   input:
#     rules.align.output
#   output:
#     "output/Alignment/QCs/{sample}_align_flagstat.qc"
#   group:
#     "Main"
#   log:
#     "logs/Alignment/{sample}_alignStat.log"
#   shell:
#     "samstat {input} > {output} 2> {log}"
 