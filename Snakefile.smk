import os
import itertools
import glob
##===================
## working directory
##===================
assembly = config['assembly']

os.chdir(assembly)
##================
## config params
##================
sample = config['samples']
index_genome_dir = config['index_genome_dir']
genome_size = config['genome_size']
genome2bit_index = config['genome2bit_index']
merged_sample = config['merged_sample']

augmented_samples = list(itertools.chain(*config['augmented_samples'])) ## flattens list of lists

blacklist = config['blacklist']
chrom_sizes = config['chrom_sizes']

##=====================
## I/O directories
##=====================
reads_input_dir = '../data/input_files/Trimmed_reads/'
base_output_dir = 'output/'
aligment_dir =  {'output': 'output/Alignment/Files/','logs':'output/Alignment/logs/','qc':'output/Alignment/qc/'}
aligment_outdir = aligment_dir['output']
aligment_logdir = aligment_dir['logs']
aligment_qcdir = aligment_dir['qc']

peakcall_dir = {'output': 'output/PeakCalling/Files/','logs':'output/PeakCalling/logs/','qc':'output/PeakCalling/qc/'}
peakcall_outdir = peakcall_dir['output']
peakcall_logdir = peakcall_dir['logs']
peakcall_qcdir = peakcall_dir['qc']

deeptools_dir = {'output': 'output/Deeptools/Files/','logs':'output/Deeptools/logs/','qc':'output/Deeptools/qc/'}
deeptools_outdir = deeptools_dir['output']
deepTools_logdir = deeptools_dir['logs']
deeptools_qcdir = deeptools_dir['qc']

##===========
## params
##===========
minQual = 20

##================
## rule groups
##================
main_group = 'main'
deeptools_group = 'deeptools'
##================
## include rules
##================
include: "rules/align.smk"
include: "rules/filter_bams.smk"
include: "rules/peakcall.smk"
# include: "rules/cross_corr.smk"
include: "rules/deeptools.smk"


rule all:
      input:  
            expand(aligment_outdir + "{sample}.bam",sample = sample), ## align
            expand(aligment_outdir + "{sample}_nochrM.bam",sample = sample), ## rm chrM mapping reads
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered.bam",sample=sample), ## include encode qc filters
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate.bam",sample=sample),## fix read mate
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread.bam",sample=sample), ## rm orphan reads 
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread_dupmark.bam",sample=sample), ## mark dup reads
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam",sample=sample), ## rm dup reads
            expand(aligment_outdir + "{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bai",sample=sample), ## index resulting filtered bams
            expand(aligment_outdir + "{sample}_tn5_shifted.bam",sample=sample),
            expand(aligment_outdir + "{merged_sample}_tn5_shifted.bam",merged_sample = merged_sample),
            expand(aligment_outdir + "{augmented_samples}_tn5_shifted_sorted.bam",augmented_samples = augmented_samples),
            expand(aligment_outdir + "{augmented_samples}_tn5_shifted_sorted.bam", augmented_samples = augmented_samples),
            expand(aligment_outdir + "{augmented_samples}_tn5_shifted_sorted.bam.bai", augmented_samples = augmented_samples),

## Peak call
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks.xls",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_summits.bed",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_treat_pileup.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_control_lambda.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks.narrowPeak",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks_filtered.narrowPeak",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks_filtered_sorted.narrowPeak.gz",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_FE.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal.bedgraph",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal_sorted.bedgraph",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal.bigwig",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_sval",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_ppois.bdg",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal.bedgraph",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal_sorted.bedgraph",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal.bigwig",augmented_samples=augmented_samples),
            expand(peakcall_qcdir + "{augmented_samples}_default.frip.txt",augmented_samples=augmented_samples),
            # "output/PeakCalling/ConsensusPeaks/consensus_peak.bed.gz"

## Cross correlation
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.fastq.gz",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.bam",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed_q30.bam",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_R1_trimmed_q30_SE.tagAlign.gz",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}.filt.sample.25Mreads.SE.tagAlign.gz",sample=sample),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.qc",sample=sample),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.plot.pdf",sample=sample),

## Deeptools 
            # expand(deeptools_outdir + "{sample}.SeqDepthNorm.bw",sample=sample),
            # deeptools_outdir + "Samples_plotCoverage.png",
            # deeptools_outdir + "Samples_plotfingerprint.png",
            # deeptools_outdir + "multiBAM_fingerprint_metrics.txt",
            # deeptools_outdir + "multiBAM_fingerprint_rawcounts.txt",
            # expand(deeptools_outdir + "{sample}_GC_content.txt",sample=sample),
            # expand(deeptools_outdir + "{sample}_plot_GC_content.png",sample=sample),
            # deeptools_outdir + "Summary.npz",
            # deeptools_outdir + "Readcounts.txt",
            # deeptools_outdir + "PearsonCor_multibamsum.png",
            # deeptools_outdir + "PearsonCor_multibamsum_matrix.txt"
