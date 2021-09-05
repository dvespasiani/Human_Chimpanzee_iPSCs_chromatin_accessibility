include: "rules/align.smk"
include: "rules/filter_bams.smk"
include: "rules/peakcall.smk"
# include: "rules/cross_corr.smk"
include: "rules/deeptools.smk"


rule all:
  input:  
## Alignment 
      expand("output/Alignment/Files/{sample}.bam",sample=config["samples"]), ## align
      expand("output/Alignment/QCs/{sample}_align_flagstat.qc",sample=config["samples"]),

## Bam filtering and Tn5 shifting 
      expand("output/Alignment/Files/{sample}_nochrM.bam",sample=config["samples"]), ## rm chrM mapping reads
      # expand("output/Alignment/QCs/{sample}_nochrM_flagstat.qc",sample=config["samples"]), 
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered.bam",sample=config["samples"]), ## include encode qc filters
      # expand("output/Alignment/QCs/{sample}_nochrM_encodefiltered_flagstat.qc",sample=config["samples"]),
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate.bam",sample=config["samples"]),## fix read mate
      # expand("output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_flagstat.qc",sample=config["samples"]),
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread.bam",sample=config["samples"]), ## rm orphan reads 
      # expand("output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_flagstat.qc",sample=config["samples"]), 
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_dupmark.bam",sample=config["samples"]), ## mark dup reads
      # expand("output/Alignment/QCs/{sample}_dup.qc",sample=config["samples"]),
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam",sample=config["samples"]), ## rm dup reads
      # expand("output/Alignment/QCs/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_dedup_flagstat.qc",sample=config["samples"]),
      expand("output/Alignment/Files/{sample}_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bai",sample=config["samples"]), ## index resulting filtered bams
      expand("output/Post_alignment/Files/{sample}_tn5_shifted.bam",sample=config["samples"]),
      expand("output/Post_alignment/Files/{sample}_tn5_shifted_sorted.bam",sample=config["samples"]),
      expand("output/Post_alignment/Files/{sample}_tn5_shifted_sorted.bam.bai", sample=config["samples"]),

## Peak call
      expand("output/PeakCalling/Files/{sample}_macs2_default_peaks.xls",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_summits.bed",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_treat_pileup.bdg",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_control_lambda.bdg",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_peaks.narrowPeak",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_peaks_filtered.narrowPeak",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_default_peaks_filtered_sorted.narrowPeak.gz",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_FE.bdg",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_fc_signal.bedgraph",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_fc_signal_sorted.bedgraph",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_fc_signal.bigwig",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_ppois_sval",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_macs2_ppois.bdg",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_ppois_signal.bedgraph",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_ppois_signal_sorted.bedgraph",sample=config["samples"]),
      expand("output/PeakCalling/Files/{sample}_ppois_signal.bigwig",sample=config["samples"]),
      expand("output/PeakCalling/QC/{sample}_default.frip.txt",sample=config["samples"]),
    # "output/PeakCalling/ConsensusPeaks/consensus_peak.bed.gz"

## Cross correlation
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.fastq.gz",sample=config["samples"]),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.bam",sample=config["samples"]),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed_q30.bam",sample=config["samples"]),
      # expand("output/ENCODE_CC/Files/{sample}_R1_trimmed_q30_SE.tagAlign.gz",sample=config["samples"]),
      # expand("output/ENCODE_CC/Files/{sample}.filt.sample.25Mreads.SE.tagAlign.gz",sample=config["samples"]),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.qc",sample=config["samples"]),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.plot.pdf",sample=config["samples"]),

## Deeptools 
      ## expand("output/DeepTools/Files/{sample}.blacklist_removed.bam",sample=config["samples"]),
      ## expand("output/DeepTools/Files/{sample}.blacklist_removed.bai",sample=config["samples"]),
      expand("output/DeepTools/Files/{sample}.SeqDepthNorm.bw",sample=config["samples"]),
      "output/DeepTools/QCs/Samples_plotCoverage.png",
      "output/DeepTools/QCs/Samples_plotfingerprint.png",
      "output/DeepTools/Files/multiBAM_fingerprint_metrics.txt",
      "output/DeepTools/Files/multiBAM_fingerprint_rawcounts.txt",
      expand("output/DeepTools/QCs/{sample}_GC_content.txt",sample=config["samples"]),
      expand("output/DeepTools/QCs/{sample}_plot_GC_content.png",sample=config["samples"]),
      "output/DeepTools/Files/Summary.npz",
      "output/DeepTools/Files/Readcounts.txt",
      "output/DeepTools/QCs/PearsonCor_multibamsum.png",
      "output/DeepTools/Files/PearsonCor_multibamsum_matrix.txt",