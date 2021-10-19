##=========================================
## 8. Run deeTools2 to get: 
## a) BAM files summary 
## b) correlation between samples
## c) fingerprint
## d) BAM coverage
##=========================================
dedup_bam = aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam"
index_bam = aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bai"

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
#    deeptools_group
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
#    deeptools_group
#   log:
#    "logs/DeepTools/{sample}_no_blacklist_index.log"
#   shell:
#    "samtools index -c {input} {output} 2> {log}"

rule deeptools_coverage:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
    deeptools_outdir + "{sample}.SeqDepthNorm.bw"
  group: 
   deeptools_group
  params:
    genome_size = genome_size
  log:
    deepTools_logdir + "{sample}_coverage.log"
  shell:
   "bamCoverage \
   --bam {input.bam} \
   --normalizeUsing RPGC \
   --effectiveGenomeSize {params.genome_size} \
   --extendReads \
   -o {output} 2> {log}"

rule deeptools_plot_coverage:
  input:
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
   deeptools_outdir + "Samples_plotCoverage.png"
  params:
   minQual = minQual,
   sample = 25000000
  group: 
   deeptools_group
  log:
   deepTools_logdir + "plotCoverage.log"
  shell:
   "plotCoverage \
   --bamfiles {input.bam}\
   --smartLabels \
   --numberOfSamples {params.sample} \
   --minMappingQuality {params.minQual} \
    -o {output} 2> {log}"   

rule deeptools_fingerprint:
  input:
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
    fig = deeptools_outdir + "Samples_plotfingerprint.png",
    metrics = deeptools_outdir + "multiBAM_fingerprint_metrics.txt",
    rawcounts = deeptools_outdir + "multiBAM_fingerprint_rawcounts.txt"
  params:
   minQual = minQual
  group: 
   deeptools_group
  log:
   deepTools_logdir + "plotFingerprint.log"
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
   content = deeptools_outdir + "{sample}_GC_content.txt",
   plot = deeptools_outdir + "{sample}_plot_GC_content.png"
  params:
   genome_size = genome_size,
   genome2bit = genome2bit_index,
   plot_format = 'png',
   threads = 5
  group: 
   deeptools_group
  log:
   deepTools_logdir + "{sample}_GC_content.log"
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
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
    summary = deeptools_outdir + "Summary.npz",
    readcounts = deeptools_outdir + "Readcounts.txt"
  threads: 5
  group: 
   deeptools_group
  log:
    deepTools_logdir + "summary.log"
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
    fig = deeptools_outdir + "PearsonCor_multibamsum.png",
    matrix = deeptools_outdir + "PearsonCor_multibamsum_matrix.txt"
  group: 
   deeptools_group
  log:
    deepTools_logdir + "correlation.log"
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
