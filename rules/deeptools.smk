
##=========================================
## 8. Run deeTools2 to get: 
## a) BAM files summary 
## b) correlation between samples
## c) fingerprint
## d) BAM coverage
##=========================================
dedup_bam="output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam"
index_bam="output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bai"

# ## First remove blacklisted regions and index resulting bam
# rule blacklist_remove:
#   input:
#     bam = "output/Alignment/Files/{sample}.nochrM.encodefiltered.fixmate.rmorphanread.nodup.bam",
#     bai = "output/Alignment/Files/{sample}.nochrM.encodefiltered.fixmate.rmorphanread.nodup.bai"
#   output:
#     temp("output/DeepTools/Files/{sample}.blacklist_removed.bam")
#   params: 
#    blacklist=config['blacklist']
#   group: 
#    "DeepTools"
#   log:
#    "logs/DeepTools/{sample}_bedtools_intersect.log"
#   shell:
#    "bedtools intersect -nonamecheck -v -a {input.bam} -b {params.blacklist} > {output} 2> {log}"

# rule noblacklist_index:
#   input:
#    rules.blacklist_remove.output
#   output:
#    temp("output/DeepTools/Files/{sample}.blacklist_removed.bai")
#   group: 
#    "DeepTools"
#   log:
#    "logs/DeepTools/{sample}_no_blacklist_index.log"
#   shell:
#    "samtools index -c {input} {output} 2> {log}"

rule deeptools_coverage:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
    "output/DeepTools/Files/{sample}.SeqDepthNorm.bw"
  group: 
   "DeepTools"
  params:
    genome_size = config['genome_size']
  log:
    "logs/DeepTools/{sample}_deepTools_coverage.log"
  shell:
   "bamCoverage \
   --bam {input.bam} \
   --normalizeUsing RPGC \
   --effectiveGenomeSize {params.genome_size} \
   --extendReads \
   -o {output} 2> {log}"

rule deeptools_plot_coverage:
  input:
    bam = expand(dedup_bam,sample=config["samples"]),
    bai = expand(index_bam,sample=config["samples"]),
  output:
   "output/DeepTools/QCs/Samples_plotCoverage.png"
  params:
   minQual=30,
   sample=25000000
  group: 
   "DeepTools"
  log:
   "logs/DeepTools/deepTools_plotCoverage.log"
  shell:
   "plotCoverage \
   --bamfiles {input.bam}\
   --smartLabels \
   --numberOfSamples {params.sample} \
   --minMappingQuality {params.minQual} \
    -o {output} 2> {log}"   

rule deeptools_fingerprint:
  input:
    bam = expand(dedup_bam,sample=config["samples"]),
    bai = expand(index_bam,sample=config["samples"]),
  output:
    fig="output/DeepTools/QCs/Samples_plotfingerprint.png",
    metrics="output/DeepTools/Files/multiBAM_fingerprint_metrics.txt",
    rawcounts="output/DeepTools/Files/multiBAM_fingerprint_rawcounts.txt"
  params:
   minQual=30
  group: 
   "DeepTools"
  log:
   "logs/DeepTools/deepTools_plotFingerprint.log"
  shell:
    "plotFingerprint \
    -b {input.bam} \
    --plotFile {output.fig} \
    --outQualityMetrics {output.metrics} \
    --outRawCounts {output.rawcounts} \
    --smartLabels \
    --minMappingQuality {params.minQual} \
    --numberOfProcessors 'max/2' 2> {log}"

rule computeGCbias:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
   content="output/DeepTools/QCs/{sample}_GC_content.txt",
   plot="output/DeepTools/QCs/{sample}_plot_GC_content.png"
  params:
   genome_size=config['genome_size'],
   genome2bit=config['2bit_index'],
   plot_format='png',
   threads=5
  group: 
   "DeepTools"
  log:
   "logs/DeepTools/{sample}_deepTools_GC_content.log"
  shell:
   "computeGCBias \
   -b {input.bam}\
   --genome {params.genome2bit} \
   --effectiveGenomeSize {params.genome_size}  \
   --plotFileFormat {params.plot_format} \
   --numberOfProcessors {params.threads} \
    -o {output.content} \
    --biasPlot {output.plot} 2> {log}" 

# keep these inputs instead as this rule will need all aligned/filtered bams
rule deeptools_summary:
  input:
    bam = expand(dedup_bam,sample=config["samples"]),
    bai = expand(index_bam,sample=config["samples"]),
  output:
    summary="output/DeepTools/Files/Summary.npz",
    readcounts="output/DeepTools/Files/Readcounts.txt"
  threads: 5
  group: 
   "DeepTools"
  log:
    "logs/DeepTools/deepTools_summary.log"
  shell:
   "multiBamSummary bins \
    -p {threads} \
    -b {input.bam} \
    -out {output.summary} \
    --outRawCounts {output.readcounts} 2> {log}"

rule deeptools_correlation:
  input:
   rules.deeptools_summary.output.summary
  output:
    fig="output/DeepTools/QCs/PearsonCor_multibamsum.png",
    matrix="output/DeepTools/Files/PearsonCor_multibamsum_matrix.txt"
  group: 
   "DeepTools"
  log:
    "logs/DeepTools/deepTools_correlation.log"
  shell:
    "plotCorrelation \
    --corData {input} \
    --plotFile {output.fig} \
    --outFileCorMatrix {output.matrix} \
    --corMethod pearson \
    --whatToPlot heatmap \
    --skipZeros \
    --plotNumbers \
    --colorMap RdYlBu 2> {log}"
