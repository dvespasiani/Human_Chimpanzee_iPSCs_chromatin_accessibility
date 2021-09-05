##==========================================================================================================================
## use this script to perform some post-alignment QCs after retaining only proper paired reads,
## e.g. following removal of chrM, duplicated and low quality (see snakefile) 
## plot fragment size distribution, return spreadsheet with NRF, PCR bottlenecks and ENCODE cc results for each sample
## TN5 read shifting (start + 4 if read maps to + strand otherwise end -5 if read maps to - strang)
##==========================================================================================================================

library(data.table);library(magrittr);library(dplyr);library(openxlsx)
library(ATACseqQC);library(ChIPpeakAnno); library(GenomicAlignments)
library(ggplot2);library(ggpubr)

# setwd("/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/")

options(width = 100,scipen=999)

alignment_dir <- "./output/Alignment/Files"
output_dir <- "./output/Post_alignment/"

## Load filtered bam file (after duplicates removal)
bamfile <- list.files(path = alignment_dir,pattern = "*nodup.bam", full.names = T, recursive = T)
# index <- list.files(path = alignment_dir,pattern = "*.bai", full.names = T, recursive = T)

samples = bamfile %>% 
lapply(
  function(x) x = sub('.*\\/', '', x)
  )%>%lapply(
    function(x)
    x = gsub("\\..*", "", x)
  )
## get bamQC metrics for each sample
bamQCs <- lapply(
  bamfile,function(x)
  x = bamQC(bamfile = x,outPath = NULL, doubleCheckDup = FALSE)
    )

## get the result of the ENCODE cross correlation analysis
phantompeak_colnames=c('Filename',
'numberReads','estimatedFragmLength','correlation_estimatedFragmLength','phantomPeak',
'correlation_phantomPeak','minCorr_strandShift','minCorr','NSC','RSC','QualityTag')

cross_corr=list.files('output/Post_alignment/ENCODE_crosscorrelation/QCs',pattern = "*cc.qc", full.names = T, recursive = T)%>%
lapply(function(x)x=fread(x,sep='\t',header=F,col.names=phantompeak_colnames))%>%rbindlist()
cross_corr=cross_corr[
  ,Filename:=stringr::str_remove(Filename,'.tagAlign.gz')
]
## make data.table with bamQC metrics you want to report
bamQC_metrics_toreport<- copy(bamQCs)%>% 
lapply(
  function(x)
  x=x[c(1:10)]%>%rbind()%>%as.data.table()
  )
names(bamQC_metrics_toreport)=samples
bamQC_metrics_toreport <- Map(mutate,bamQC_metrics_toreport,'Sample'=samples)%>%rbindlist()
bamQC_metrics_toreport <- bamQC_metrics_toreport %>% dplyr::select(c('Sample',everything())) %>% mutate_if(is.list,as.numeric)

## here combine the QC files and report a single excel spreadsheet
qc_spreadsheet=list(bamQC_metrics_toreport,cross_corr)
names(qc_spreadsheet)=c('PostAlignmentQC','ENCODE_CrossCorr')
write.xlsx(qc_spreadsheet,paste(output_dir,'ATACseqQC/Samples_qc_metrics.xlsx',sep=''))

## make a spreadsheet with all MapQ values and the relative number of reads for every sample
mapq_qc=copy(bamQCs)%>%lapply(function(x)x=x[[11]])
mapq_qc=lapply(
  mapq_qc,
  function(x)x=x%>%as.data.table()%>%setnames(c('MAPQ','Number_Reads')) %>%
  mutate_if(is.list,as.numeric)
)
names(mapq_qc)=samples
write.xlsx(mapq_qc,paste(output_dir,'ATACseqQC/Samples_MAPQ_scores.xlsx',sep=''))

##==================================================================================
## For each sample report a summary of the reads mapping per chromosome 
## PS: group all *_random and chrUn* seqnames into a single category
##==================================================================================

bam_idxstats=copy(bamQCs)%>% lapply(function(x) x=x$idxstats%>%as.data.table())
bam_idxstats=lapply(bam_idxstats,function(x) x=x[,total_mapped_reads:=sum(mapped)])

chrnames=paste('chr',c(seq(1:22),'X','Y'),sep='')

other_chrs=copy(bam_idxstats)%>%
lapply(function(x)
x=x[
  !seqnames%in%chrnames,mapped:=sum(mapped)
  ][
  !seqnames%in%chrnames,unmapped:=sum(unmapped)
  ][
    !seqnames%in%chrnames
    ][,c('mapped','unmapped','total_mapped_reads')]%>%unique()
)
other_chrs=lapply(other_chrs,function(x)x=x[,seqnames:='other_chrs'][,seqlength:=0]%>%
dplyr::select(c('seqnames','seqlength',everything())))

bam_idxstats <- lapply(bam_idxstats,function(x) x=x[seqnames%in%chrnames])
final_bam_idxstats <-  purrr::map2(bam_idxstats,other_chrs,rbind)

final_bam_idxstats <- lapply(final_bam_idxstats,function(x)x=x[,fraction_mapped_reads_per_chr:=mapped/total_mapped_reads])
names(final_bam_idxstats)=samples
final_bam_idxstats <- Map(mutate,final_bam_idxstats,'Sample'=names(final_bam_idxstats))%>%rbindlist()

##==================================================================================
## For each sample plot the raw number of reads mapping per every chromosome 
## and their fraction across all mapped reads
##==================================================================================
plot_mapped_reads=function(x){
fraction_reads=ggplot(x,aes(x=seqnames,y=fraction_mapped_reads_per_chr,fill=Sample))+
geom_bar(stat='identity',position ='dodge')+
xlab('')+ylab('Fraction of mapped reads')+
theme(axis.text.x = element_text(angle=50,vjust = 0.9,hjust = 1))

raw_reads=ggplot(x,aes(x=seqnames,y=mapped,fill=Sample))+
geom_bar(stat='identity',position ='dodge')+
xlab('')+ylab('Raw number of mapped reads')+
theme(axis.text.x = element_text(angle=50,vjust = 0.9,hjust = 1))

final=ggarrange(raw_reads, fraction_reads,ncol = 1)
return(final)
}

pdf(paste(output_dir,'ATACseqQC/Plot_Mapped_reads.pdf',sep=''),width=7,height = 7)
plot_mapped_reads(final_bam_idxstats)
dev.off()

## Now plot the fragment size distribution for every sample
pdf(paste(output_dir,'ATACseqQC/Plot_fragment_size_distribution.pdf',sep=''),width=7,height = 7)
fragSizeDist(bamfile,samples)
dev.off()

##=====================================
## TN5 read shifting and write files
##=====================================
## read bam files, shift the reads and then return the shifted files
shifted_reads_dir <- paste(output_dir,"Files/",sep='')

tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
original_bamfiles <- lapply(bamfile, function(x)readBamFile(x, tag=tags, asMates=TRUE, bigFile=T))

shifted=copy(original_bamfiles)
names(shifted)=paste(shifted_reads_dir,samples,sep='')
shifted_filenames <- file.path(names(shifted), "_tn5_shifted.bam",fsep='')

shifted_bamfiles=lapply(shifted,function(x)shiftGAlignmentsList(x))
mapply(rtracklayer::export,shifted_bamfiles, con = shifted_filenames,format='bam')

