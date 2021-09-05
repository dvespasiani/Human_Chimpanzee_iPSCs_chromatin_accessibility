# Human vs Chimp iPSCs ATAC-seq 
This repo contains all scripts used to analyse ATAC-seq files from <a href='https://www.biorxiv.org/content/10.1101/466631v1.full.pdf'>Gallego Romero *et al.,* 2018<a>. Files accession number: [GSE122319](https://www.ncbi.nlm.nih.gov/gds/?term=GSE122319[Accession]).

## Table of Contents
1. [Set up environment & download files](#Set-up-environment-&-download-files)
2. [Snakemake pipeline](#Snakemake-pipeline)
   - [FastQC](#FastQC)
   - [Adaptor trimming](#Adaptor-trimming)
   - [Genome alignment](#Genome-alignment)
   - [Post-alignment QCs](#post_align_qc)
   - [Strand cross-correlation](#encode_cc)
   - [Deeptools QCs](#Deeptools-QCs)
   - [Peak call](#Peak-call)
   - [Peak call QCs](#Peak-call-QCs)
____
## Set up environment & download files
To set up the environment and get ready to go look into the `workflow_utils` dir, which contains some useful scripts, including:
1. `download_files.sh` to download both SRA and ENCODE blacklisted files. The list of SRA files can be found in `data/SRaAccList.txt`. ENCODE blacklisted regions for hg38 are downloaded from [here](http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz). Genomic coordinates are then converted from hg38 to panTro5 using liftover 
2. `create_conda_env.sh` to...create the conda env for running the snakemake pipeline
3. `subsample_files.sh` to.... subsample the SRA files and creating a new directory where you can test the snakemake pipeline without waiting for ages
_____

# Snakemake pipeline
If not otherwise specified, all these steps below have been taken from the ENCODE ATAC-seq standard processing [pipeline](https://www.encodeproject.org/atac-seq/). For full protocol specifications check [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit).

## FastQC
When running FastQC, you can expect 3 modules returning a warining/failure signal:
1. Per base sequence content because Tn5 has a strong sequence bias at the insertion site. 
2. Sequence Duplication Levels caused by PCR duplicates
3. Overrepresented sequences

## Adaptor trimming
Contrarily to what ENCODE did, for this step I used Trimmomatic to remove Illumina Nextera adapter sequeces. Fasta files containing the actual adapter sequences were obtained from  [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa).<br/>
To identify adapters I've allowed for 2 max mismatches, a threshold of Q = 30 for PE palindrome read alignment (this can control for short adapter retentions at the 3' end of each read) and a threshold of Q = 10 for a simple alignment match between adapters and read. <br />
PS: In trimmomatic each match increase Q score by 0.6 whereas mismatches reduce Q by Q/10. Thus for a Q of 30 there must be around 50 matches and for a Q of 10 there are around 16 matches. <br />
Finally, I've removed all initial/terminal sequences having a phred score <20 and all reads that, after these quality steps were <20bp long. <br/>
Trimmomatic then splits fastq files in `*.paired.fastq` and `*.unpaired.fastq`. All downstream analyses are carried out only on paired trimmed reads.

## Genome alignment 
I used Bowtie2 to align reads to hg38 and panTro5 genomes. Differences with ENCODE alignment:
* I have specified a max fragment length of 5kb (encode sets it to 2kb)
* I have allowed dovetail (i.e. mates extending past each other)
* The maximum number of distinct, valid alignments was set to 1, so no multimapping was allowed (encode offers both solutions). Anyhow, Bowtie2 by default reports the best alignment, and in case of multimapping with same score a random selection is performed.
* I have ran Bowtie2 on a very sensitive mode (i.e., `--very-sensitive`) which by default means:
      * No mismatches are allowed 
      * The length of the seed substrings to align during multiseed alignment is 20 (ps: smaller values make alignment more sensitive)
      * Up to 20 consecutive seed extension attempts that can fail before Bowtie2 moves on, using the best alignment found (nb: a seed extension fails if it does not yield a new best or a new second-best alignment)
      * A maximum of 3 re-seeds of reads with repetitive seeds is allowed (i.e. Bowtie2 chooses a new set of reads with same length and same number of mismatches at different offsets and searches for more alignments).
See [Bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-i) for a more detailed description. <br/>
**NB**: Concordant alignment occurs when reads aligned with their expected mate orientation and with the expected range of distance between mates. Discordant alignment instead occurs when both mates have unique alignments, but the mates aren't in their expected relative orientations or aren't within their expected range ditance (or both).

## Post-alignment QCs <a name="post_align_qc"></a>
Following alignment BAM files are processed using samtools and/or picard tools to: 
1. remove reads mapping to chrM
2. keep only proper paired reads (`-F 1804 -f 2`) and reads with MAPQ score >= 30 (`-q 30`). For flag explanation check [here](https://broadinstitute.github.io/picard/explain-flags.html).
3. remove reads with mate unmapped (`samtools fixmate`)
4. remove (again) read pairs mapping to different chromosomes
5. remove duplicated reads 

PS: step 3 and 4 seem redundant also by looking at the number of reads retained after each step. However, since encode does it and for the moment it is not killing anythng I'll keep it. <br/>

For every step a flagstaf QC file is generated reporting number and type of reads filtered out. However, to get a spreadsheet with the number of reads for each SRA sample after every step execute `common_scripts/make_fileQC_spreadsheet.sh`. <br/>
Finally, for each sample I've indexed the resulting BAM file.  <br/>

**NB:** library complexity and fragment length statistics are calculated in later steps.
 
## Strand cross-correlation <a name="encode_cc"></a>
Cross-correlation (CC) doesnt seem to be vital for ATAC-seq (I'll find again Anshul's comment on this), but I will still include it. For info regarding the rationale of this step read [Landt *et al.,* 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/) or [this webpage](https://genome.ucsc.edu/ENCODE/qualityMetrics.htm). For a detailed description of the CC plot look [here](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/CC_metrics_extra.html) instead. <br/>

In summary, high-quality ChIP-seq experiments would return significant clustering of enriched DNA sequence tags/reads at locations bound by the protein of interest. Such enrichmnet is to be expected both on the + and - strands.
The strand CC is a measure of the enrichment derived without dependence on prior determination of the enriched regions. It is computed as the Pearson’s linear correlation between the coverage for each complementary base (on the - strand and the + strand), by systematically shifting the - strand by *k* bp at a time (i.e. reads are incrementally shifted away from each other). The resulting CC profile will represent the correlation between the number of mapping reds on the 2 strands for a given genomic region. Finally, for any given sample,the correlation values are calculated across every peak for each chromosome, multiplied by a scaling factor and then summed across all chromosomes. The highest CC value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of the clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor or feature. <br/> 

The resulting plot shows the CC scores against the shift value. This plot generally returns 2 enrichment peaks corresponding to:
* the predominant fragment length (highest correlation value)
* the read length (or phantom peak) <br/>
For short-read datasets (< 100 bp reads) and large genomes with a significant number of non-uniquely mappable positions (e.g., human and mouse), the cross-correlation phantom-peak is observed at a strand-shift equal to the read length. This read-length phantom peak is an effect of the variable and dispersed mappability of positions across the genome. For a significantly enriched dataset, the fragment length cross-correlation peak (red line), representing clustering of fragments around target sites, should be larger than the mappability-based read-length phantom peak (blue line)

Metrics that are calculated and that assess the signal-to-noise ratios in a ChIP-seq experiments are:
1. Normalized strand coefficient (NSC) <br/>
This is the ratio of the maximal cross-correlation value (which occurs at strand shift equal to fragment length) divided by the background cross-correlation (i.e. the minimum cross-correlation value over all possible strand shifts. NSC values range from a minimum of 1 (no enrichment) to larger positive numbers with 1.1 considered as the critical threshold. Datasets with NSC values < 1.1 (e.g. < 1.05) tend to have low signal to noise or few peaks. This score is sensitive to both technical effects, e.g. high-quality antibodies such as H3K4me3 and CTCF score well for all cell types and biological effects, e.g. narrow marks score higher than broad marks such as H3K4me3 vs H3K36me3, H3K27me3.

2. Relative strand correlation (RSC) <br/>
This is the ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. RSC values range from a minimum of 0 (no signal) to larger positive values (highly enriched experiments) with 1 considered as the critical threshold. RSC values significantly lower than 1 (e.g. < 0.8) tend to have low signal to noise thus low quality experiments.

3. Qtag (a thresholded version of RSC) 

## DeepTools QCs
DeepTools were used for:
1. Get bam summary and coverage
2. Plot the correlation between bam coverages
3. Calculate the GC bias
3. Calculate the cumulative enrichment for each sample (i.e. fingerprint)
For a complete description of the tools see [deepTools manual](https://deeptools.readthedocs.io/en/develop/index.html). <br/>
Importantluy, some of the DeepTools command require you to specify the effective the genome size. For hg38 this was easily found [here](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html). However for panTro5 I've used `pantro5/scripts/unique-kmers.py` from [MR Crusoe *et al.*, 2015](http://dx.doi.org/10.12688/f1000research.6924.1) as follows:
```
python unique-kmers.py -k 50 /data/projects/punim0586/shared/genomes/panTro5/panTro5.fa
```
which returned
```
Total estimated number of unique 50-mers: 2792339170
```
Importantly, from the filtered BAM files the series of I/O files necessary for deepTools analyses constitute a specific running group that do not affect downstream processes such as peak calling. They represent a separate branch that allow to perform post-alignment QCs similarly to the ENCODE CC analysis.

### BAM Summary and Coverage (FROM HERE IT NEEDS TO BE UPDATED)
Following BAM filtering I removed reads overlapping ENCODE blacklisted regions and then calculated the bam coverage with `bamCoverage` for all samples. Coverage is calculated as the number of reads per bin (i.e. short consecutive counting windows of defined size) and it can be scaled in RPKM,CPM. <br/>
Resulting summary was plotted using `plotCoverage` to visually assess the sequencing depth of each sample after sampling 25*10^6 reads. This command counts the number of overlapping reads and returns 2 plots indicating:
1. The frequencies of the observed read coverages per sample
2. The fraction of the genome covered by >= a given number of reads 

To investigate the correlation between the samples read coverages I've used `multiBamSummary`, which computes the read coverages over the entire genome and/or specific genomic regions for >=2 BAM files. Resulting pairwise Pearson correlation values are plotted into a heatmap  with `plotCorrelation`. 

### GC bias
This QC step tests the assumption of an expected uniform distribution of sequenced reads across the genome, regardless of their base-pair composition. However, PCR steps performed during library preparation can enrich for GC-rich fragments and this will influence the sequencing results.<br/>
The deepTools command `computeGCbias` first calculates the expected GC profile from the given genome by calculating the GC % for DNA fragments of a fixed size. It then compares the expectations with the observed GC content within the sequenced reads. <br/> 
The resulting plots report the number of reads overlapping genomic regions of increasing GC content and the log2(observed/expected) ratio of reads per GC content region. If there is no GC bias then the observed and expected GC profiles will be similar, i.e. ratio of observed/expected GC content per fragment length would = 1. <br/>
In this experiment there is a linear increase in the exp/obs read coverage across all samples. Possible solutions:
1. Run a simple, non-conservative peak calling on the uncorrected BAM file first to obtain a BED file of peak regions that are then supplied to `computeGCbias` to limit the calculation
2. Run the `correctGCbias` deepTool, which removes reads from regions of too high-coverage and add reads to regions of low-coverage
2. Leave it as it is since the bias is consistent across all samples and it might be a biological effect 

### Cumulative enrichment (BAM fingerprint)
<p> The deepTools `plotFingerprint` samples indexed BAM files and calculates the samples cumulative enrichment by counting all reads overlapping a window (bin) of a given length. The resulting plot is useful to assess how well the signal of the sequenced reads can be differentiated from the background distribution of reads in the control input sample. <br/>

An ideal input with perfect uniform distribution of reads along the genome (i.e. without enrichments in specific open chromatin regions) and infinite sequencing coverage would generate a straight diagonal line. On the other hand, a strong enrichment in very specific (i.e. fragment size) chromatin regions will result in a prominent and steep rise of the cumulative sum towards the highest rank. This means that a big chunk of reads from the test (ChIP/ATAC) sample is located in few bins which corresponds to high, narrow enrichments typically seen for transcription factors. <br/>

**PS**: check where cumulative curve starts on the plots as this will give you an indication of the percentage of the genome that was not sequenced at all (i.e. bins containing 0 reads). However, such percentage can be quite high, especially for extremely high enrichment which will result in the vast majority of reads occurring within few peaks, [see examples](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html).  </p>

## Peak call
Peak calling was performed using MACS2 on Tn5-shifted BAM files (from `ATACseqQC.R` script). Because coverage is very low for each sample I first merged all the technical replicates and then called the peaks on each sample (i.e. biological replicate).  <br/>
MACS2 flags set according to ENCODE guidelines.
Output files:
* $_default_peaks.narrowPeaks
BED6+4 file format. 
   * 5th column (i.e. score) is calculated as `int(-10*log10pvalue)`<br/>
   * 7th colum reports fold-change at peak summit <br/>
   * 8th and 9th colums report -log10pvalue and -log10qvalue respectively at peak summit <br/>
   * 10th column represents the relative summit position to peak start <br/>
* $_default_peaks.xlsx
Spreadsheet containing info about called peaks (i.e. chr; start; end; peak summit position; p/q-vals etc...)
* $_default_summits.bed
BED file containing the peak summits locations for every peak
* $_default_treat_pileup.bdg
bedGraph file containing the pileup signals of the treatment sample
* $_control_lambda.bdg
bedGraph file containing the local biases estimated for each genomic location from the control sample, or from treatment sample when the control sample is absent
<p>
These latter 2 files are then compared using `bdgcmp` subcommand to generate a signal track bedGraph file containing scores: p-value, q-value, log-likelihood, and log fold changes for each peak.

Chimp `chrom.sizes` file for bedtools was generated using:
```
samtools faidx panTro5.fa
cut -f 1,2 panTro5.fa.fai > chrom.sizes
```
whereas to get chrom.sizes for hg38 chromosomes simply run:
` wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes`

## Peak call QCs
### Fraction Reads in Peaks (FRiP)
It represents the proportion of all mapped reads that fall into the called peak regions. FRiP scores positively correlates with the number of regions. According to ENCODE, FRiP should be >0.3, though values greater than 0.2 are acceptable.

### TSS enrichment
This metric is used as another signal to noise indication. Reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp). This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome)there should be an increase in signal up to a peak in the middle. We take the signal value at the center of the distribution after this normalization as our TSS enrichment metric.

## Testing for differential accessible (DA) chromatin  (need to re-write this nicely)
To assess whether the accessible chromatin landscapes differs between human and chimp samples I followed the indications of [Reske et al. 2020](https://github.com/reskejak/ATAC-seq/blob/master/csaw_workflow.R), check script `common_scripts/csaw_DA_analysis.R`. Also see [csaw manual](https://bioconductor.org/packages/3.12/workflows/vignettes/csawUsersGuide/inst/doc/csaw.pdf)<br/>

Importantly, before normalising count data and performing the differential accessible analysis I've used liftover to convert MACS2 peak files from one species to the other and then I've used liftover to align back the peaks to the original (species) genome build. I've then discarded those peaks whose widths was < 20% of their original one.

### Count normalisation
`common_scripts/csaw_DA_analysis.R` evaluates both TMM and LOESS normalisation methods. I will pick the most suitable one. 

### Retrieving DA regions
the above script returns multiple plots and the list of DA peaks between the species. Reads from tn5-shifted and black-list removed bam files are counted if overlapping peaks.
First constructed the DGEList object (edgeR) and computed the average log2 CPM for each row of counts.
I've removed regions with low abundance of reads (i.e., logCPM > -3) as this migth reduce the severity of the multiple testing correction.
For each bam file, using a sliding window approach, I then count the fragments overlapping a fixed-width genomic interval. I then computed the TMM and/or LOESS normalization factors based on binned counts.
I then designed the matrix, estimated the dispersion parameters, fitted the model to the data and finally tested each region for DA between human and chimp (all through edgeR) 


