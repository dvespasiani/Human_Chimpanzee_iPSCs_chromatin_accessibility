## ===========================================================================================
##  Post-alignment QCsto remove unmapped, duplicated, multimapped and low quality reads
##  Removed also multimapping reads (MAPQ <=30) 
## ===========================================================================================
rule rmChrM:
  input:
   rules.align.output
  output:
   "output/Alignment/Files/{sample}_nochrM.bam"
  group:
    "Main"
  log:
    "logs/Alignment/{sample}_rmchrM.log"
  shell:
   "samtools view -h {input} | grep -v 'chrM' | samtools sort -o {output} 2>{log}"

# rule rmchrMStat:
#   input:
#     rules.rmChrM.output
#   output:
#     "output/Alignment/QCs/{sample}_nochrM_flagstat.qc"
#   group:
#     "Main"
#   log:
#     "logs/Alignment/{sample}_nochrM_qc.log"
#   shell:
#     "samstat {input} > {output} 2> {log}"
 
rule encode_filters:
  input:
   rules.rmChrM.output
  output:
    "output/Alignment/Files/{sample}_nochrM_encodefiltered.bam"
  group:
    "Main"
  log:
   "logs/Alignment/{sample}_encodefiltered.log"
  shell:
   "samtools view -b -F 1804 -q 30 -f 2 -u {input} | samtools sort -n -o {output}  2> {log}"

# rule encode_filtersStat:
#   input:
#     rules.encode_filters.output
#   output:
#     "output/Alignment/QCs/{sample}_nochrM_encodefiltered_flagstat.qc"
#   group:
#     "Main"
#   log:
#     "logs/Alignment/{sample}_encodefiltered_qc.log"
#   shell:
#     "samstat {input} > {output} 2> {log}"

rule fixmate:
  input:
   rules.encode_filters.output
  output:
   "output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate.bam"
  group:
    "Main"
  log:
   "logs/Alignment/{sample}_fixmate.log"
  shell:
   "samtools fixmate -r {input} {output} 2> {log}"

# rule fixmateStats:
#   input:
#    rules.encode_filters.output
#   output:
#    "output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_flagstat.qc"
#   group:
#     "Main"
#   log:
#    "logs/Alignment/{sample}_fixmateStats.log"
#   shell:
#    "samstat {input} > {output} 2> {log}"

rule rmOrphanread:
  input:
   rules.fixmate.output
  output:
   "output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread.bam"
  group:
    "Main"
  log:
   "logs/Alignment/{sample}_rmorphanread.log"
  shell:
   "samtools view -F 1804 -f 2 -u {input} | samtools sort -o {output} 2> {log}"

# rule rmOrphanreadStats:
#   input:
#    rules.fixmate.output
#   output:
#    "output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_flagstat.qc"
#   group:
#     "Main"
#   log:
#    "logs/Alignment/{sample}_rmOrphanreadStats.log"
#   shell:
#    "samstat {input} > {output} 2> {log}"

##============================
## Mark and remove duplicates
##============================
## remember to sort samfile by qname for this
## keep parameters USE_JDK_DEFLATER=true USE_JDK_INFLATER=true otherwise there could be java incompatibilities
rule markDups:
  input:
    rules.rmOrphanread.output 
  output:
    bam="output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_dupmark.bam",
    dupQC="output/Alignment/QCs/{sample}_dup.qc"
  group:
    "Main"
  params:
   mem= "-Xmx4g"
  log:
    "logs/Alignment/{sample}_dupmark.log"
  shell:
   "java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={input} O={output.bam} \
   METRICS_FILE={output.dupQC} VALIDATION_STRINGENCY=LENIENT \
   ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log}"

rule dedup:
  input:
    rules.markDups.output.bam
  output:
    "output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam"
  group:
    "Main"
  log:
    "logs/Alignment/{sample}_nodup.log"
  shell:
    "samtools view -h -b  -F 1804 -f 2 {input} > {output} 2> {log}"

# rule dedupStat:
#   input:
#     rules.dedup.output
#   output:
#     "output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_dedup_flagstat.qc"
#   group:
#     "Main"
#   log:
#     "logs/Alignment/{sample}_dedup_qc.log"
#   shell:
#     "samstat {input} > {output} 2> {log}"

##============================
## 6.Index bam file
##============================
rule indexBam:
  input:
    rules.dedup.output
  output:
    "output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bai"
  group:
    "Main"
  log:
    "logs/Alignment/{sample}_indexBam.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"

## Tn5 shift reads
rule Tn5_shift:
  input:
    bam=rules.dedup.output,
    bai=rules.indexBam.output
  output:
    "output/Post_alignment/Files/{sample}_tn5_shifted.bam"
  group:
    "Main"
  log:
    "logs/Post_alignment/{sample}_tn5_shift.log"
  shell:
    "alignmentSieve -b {input.bam} -o {output} --ATACshift 2> {log}"

rule Tn5_shifted_sort:
  input:
    rules.Tn5_shift.output
  output:
    "output/Post_alignment/Files/{sample}_tn5_shifted_sorted.bam"
  group:
    "Main"
  log:
    "logs/Post_alignment/{sample}_tn5_shift_sorted.log"
  shell:
    "samtools sort {input} -O BAM -o {output} 2> {log}"

rule index_Tn5Bams:
  input:
    rules.Tn5_shifted_sort.output
  output:
    "output/Post_alignment/Files/{sample}_tn5_shifted_sorted.bam.bai"
  log:
    "logs/Post_alignment/{sample}_tn5_shifted_sorted_indexBam.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"
