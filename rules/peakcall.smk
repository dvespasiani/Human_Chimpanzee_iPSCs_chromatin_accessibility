##=========================================
## 10. Peak calling (MACS2)
##=========================================
rule peakCalling:
  input:
    bam = rules.Tn5_shifted_sort.output,
    bai = rules.index_Tn5Bams.output
    # bai = "output/Post_alignment/Files/{augmented_samples}_tn5_shifted_sorted.bam.bai"
  output:
    peaks_xls = peakcall_outdir + "{augmented_samples}_macs2_default_peaks.xls",
    pileup = temp(peakcall_outdir + "{augmented_samples}_macs2_default_treat_pileup.bdg"),
    lamb = temp(peakcall_outdir + "{augmented_samples}_macs2_default_control_lambda.bdg"),
    summit = temp(peakcall_outdir + "{augmented_samples}_macs2_default_summits.bed"),
    narrowPeak = peakcall_outdir + "{augmented_samples}_macs2_default_peaks.narrowPeak"
  group: 
    main_group
  params:
    name = "{augmented_samples}_macs2_default",
    outdir = peakcall_outdir,
    fragment_size = 300,
    shift = 100,
    genome_size = genome_size,
    pval_thresh = 0.01
  log:
    peakcall_logdir + "{augmented_samples}_peak_calling.log"
  shell:
    "macs2 callpeak --format BAMPE --treatment {input.bam} \
    --keep-dup all \
    --outdir {params.outdir} --name {params.name} --shift {params.shift} \
    --nomodel -B --SPMR --extsize {params.fragment_size} \
    --pvalue {params.pval_thresh} --call-summits -g {params.genome_size} 2> {log}"

##--------------------------------------------
## Remove peaks within black listed regions
##--------------------------------------------
# rule rm_black_peaks:
#   input:
#    rules.peakCalling.output.narrowPeak
#   output:
#    peak_dir + "{augmented_samples}_macs2_default_peaks_filtered.narrowPeak"
#   params:
#    blacklist = blacklist
#   group: 
#    main_group
#   log:
#    "logs/PeakCalling/{augmented_samples}_blacklist_removed.log"
#   shell:
#    """
#    bedtools intersect -v -a {input} -b {params.blacklist} \
#    | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' \
#    | awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y|2A|2B)$/' > {output} 2> {log}
#    """
## As suggested by ENCODE, sort NarrowPeaks by column 8 in descending order 
## and replace long peak names in column 4 with Peak_<peakRank>
## keep only first NPEAKS
# rule peakSort:
#   input:
#    rules.rm_black_peaks.output
#   output:
#    peak_dir + "{augmented_samples}_macs2_default_peaks_filtered_sorted.narrowPeak.gz"
#   log:
#    "logs/PeakCalling/{augmented_samples}_sorted_narrowpeaks.log"
#   params:
#     NPEAKS=300000
#   group: 
#     main_group
#   shell:
#     """
#     sort -k 8gr,8gr {input} | \
#     awk 'BEGIN{{OFS="\\t"}}{{$4="Peak_"NR ; print $0}}' \
#     | head -n {params.NPEAKS} | \
#     gzip -nc > {output} 2> {log}
#     """
##------------------------------------------------
## Get peak Fold Enrichments (FE) and 
## then convert it to bigWig format for display
##------------------------------------------------
rule FE_peak_signal_tracks:
  input:
   treatment = rules.peakCalling.output.pileup,
   control = rules.peakCalling.output.lamb
  output:
   temp(peakcall_outdir + "{augmented_samples}_macs2_FE.bdg")
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_Signaltrack.log"
  shell:
   """
   macs2 bdgcmp -t {input.treatment} \
   -c {input.control} \
   --ofile {output} \
   -m FE \
   2> {log}
   """

rule clean_FE_signal:
  input:
   rules.FE_peak_signal_tracks.output
  output:
   temp(peakcall_outdir + "{augmented_samples}_fc_signal.bedgraph")
  log:
   peakcall_logdir + "{augmented_samples}_cleanFE_peaks.log"
  params:
    chrom_sizes = chrom_sizes
  group: 
    main_group
  shell:
   """
   slopBed -i {input} -g {params.chrom_sizes} -b 0 \
   | bedClip stdin {params.chrom_sizes} {output}  2> {log}
   """

rule sort_FE_bedGraph:
  input:
   rules.clean_FE_signal.output
  output:
   temp(peakcall_outdir + "{augmented_samples}_fc_signal_sorted.bedgraph")
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_sort_bedGraph.log"
  shell:
   "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule FE_bedGraph2bigWig:
  input:
   rules.sort_FE_bedGraph.output
  output:
   peakcall_outdir + "{augmented_samples}_fc_signal.bigwig"
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_bedGraph2bigWig.log"
  params:
    chrom_sizes = chrom_sizes
  shell:
   """
   bedGraphToBigWig {input} {params.chrom_sizes} {output} 2> {log}
   """

##------------------------------------
## now repeat these latter steps 
## but to calculate peak p-values
##-----------------------------------

## sval counts the number of tags per million in the (compressed) BED file
rule sval:
  input:
  #  rules.peakSort.output # change this with MACS output
   rules.peakCalling.output.narrowPeak
  output:
   peakcall_outdir + "{augmented_samples}_ppois_sval"
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_sval_calculation.log"
  shell:
   """
   (wc -l <(gunzip -nc {input}) \
   | awk '{{printf "%f", $1/1000000}}' >{output} ) 2> {log}
   """

rule ppois_peak_signal_tracks:
  input:
  #  sval = peakcall_outdir + "{augmented_samples}_ppois_sval",
  #  treatment = peakcall_outdir + "{augmented_samples}_macs2_default_treat_pileup.bdg",
  #  control = peakcall_outdir + "{augmented_samples}_macs2_default_control_lambda.bdg"
   sval = rules.sval.output,
   treatment = rules.peakCalling.output.pileup,
   control = rules.peakCalling.output.lamb
  output:
   temp(peakcall_outdir + "{augmented_samples}_macs2_ppois.bdg")
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_ppois_signaltrack.log"
  run:
   for filename in glob.glob(peakcall_outdir + '*_ppois_sval'):
     sval = open(filename).read()
   shell("""
   macs2 bdgcmp -t {input.treatment} -c {input.control} --ofile {output} \
   -m ppois -S {sval} 2> {log}
   """)

rule clean_ppois_signal:
  input:
   rules.ppois_peak_signal_tracks.output
  output:
   temp(peakcall_outdir + "{augmented_samples}_ppois_signal.bedgraph")
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_clean_ppois_peaks.log"
  params:
    chrom_sizes = chrom_sizes
  shell:
   """
   slopBed -i {input} -g {params.chrom_sizes} -b 0 \
   | bedClip stdin {params.chrom_sizes} {output} 2> {log}
   """

rule sort_ppois_bedGraph:
  input:
   rules.clean_ppois_signal.output
  output:
   temp(peakcall_outdir + "{augmented_samples}_ppois_signal_sorted.bedgraph")
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_sort_ppois_bedGraph.log"
  shell:
   "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule ppois_bedGraph2bigWig:
  input:
   rules.sort_ppois_bedGraph.output
  output:
   peakcall_outdir + "{augmented_samples}_ppois_signal.bigwig"
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_ppois_bedGraph2bigWig.log"
  params:
    chrom_sizes = chrom_sizes 
  shell:
   """
   bedGraphToBigWig {input} {params.chrom_sizes} {output} 2> {log}
   """

##----------------------------------------------------
## FRiP on pooled bam files and filtered narrowPeaks
##----------------------------------------------------
rule frip:
  input:
   pooledbam = rules.Tn5_shifted_sort.output,
   pooledpeaks = rules.peakCalling.output.narrowPeak
  #  pooledpeaks = rules.peakSort.output # change this with MACS output
  output:
   peakcall_qcdir + "{augmented_samples}_default.frip.txt"
  group: 
    main_group
  log:
   peakcall_logdir + "{augmented_samples}_calculate_FRiP.log"
  run:
   shell("python3 ../common_scripts/encode_frip.py \
    {input.pooledbam} {input.pooledpeaks} > {output} 2> {log}")

# ##------------------------
# ## Call consensus peak
# ##------------------------
# rule consensus_peak:
#   input:
#    peak_merged_sample = peak_dir + "merged_macs2_default_peaks_filtered_sorted.narrowPeak.gz",
#    peak_other_samples = expand(peak_dir + "{augmented_samples}_macs2_default_peaks_filtered.narrowPeak.gz", sample=config["samples"])
#   output:
#    "output/PeakCalling/ConsensusPeaks/consensus_peak.bed.gz"
#   params:
#    min_overlap = 0.5,
#    sample_names = config["samples"]
#   group:
#     main_group
#   log:
#    "logs/PeakCalling/Consensus_peak.log"
#   shell:
#    """
#    bedtools intersect -f {params.min_overlap} -r \
#    -a {input.peak_merged_sample} -b {input.peak_other_samples} \
#    | gzip -nc > {output} 2> {log}
#    """
