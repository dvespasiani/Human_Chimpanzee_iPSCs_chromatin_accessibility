## input: the human and chimpanzee peaks called on the merged bams in the respective coordinates 
## output: the set of peaks that have orthologous regions in the other species
## qc: check how many individuals support each consensus peak
## use liftover back and forth

library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(UpSetR)
library(GenomicRanges)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility')

scripts_dir <- './post_processing_analyses/scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

peakDir = "output/PeakCalling/Files/"
outplot_dir <- create_dir(plot_dir, 'consensusPeaks')
outfile_dir <- create_dir(outdir,'files/consensusPeaks')

##--------------
## Read peaks 
##--------------
## consensus peaks
read_consensus_peaks <- function(genome,peakfile){
    peaks <- fread(paste(genome,peakDir,peakfile,sep=''),sep='\t',header=F,select=c(1:3))%>%
        setnames(c(range_keys))%>%makeGRangesFromDataFrame()%>%reduce()%>%as.data.table()
    peaks <- peaks[,peakID:=paste('peak_',1:nrow(peaks),sep='')][,c('width','strand'):=NULL]
    setkeyv(peaks,range_keys)
    return(peaks)
}
human_consensus_peaks <- read_consensus_peaks('./hg38/','human_hg38_macs2_default_peaks.narrowPeak')
chimp_consensus_peaks <- read_consensus_peaks('./panTro5/','chimp_pantro5_macs2_default_peaks.narrowPeak')


## individual peaks
read_individual_peaks = function(species,patterns){
    file_path=dir(paste(species,'/output/PeakCalling/Files',sep=''),pattern=patterns,full.names=T)
    file_names=gsub('\\_.*','',basename(file_path))
    files=lapply(
        file_path,function(x)
        fread(x,sep='\t',header=F,select=c(1:3),col.names=range_keys)[
            !seqnames %like% '_'
            ]%>%setorderv(range_keys)%>%unique())
    files = Map(mutate,files,sample=file_names)
    lapply(files,function(x)setkeyv(x,range_keys))   
    names(files)=file_names
    return(files)
}

chimp_individual_peaks = read_individual_peaks('./panTro5/',"^C.*macs2_default_peaks.narrowPeak$")
human_individual_peaks = read_individual_peaks('./hg38/',"^H.*macs2_default_peaks.narrowPeak$")

peak_intersection <- function(peaks,consensus){
    peak_intersection=lapply(peaks,function(x)x=foverlaps(x,consensus,type='any')%>%na.omit())
    peakIDs=copy(peak_intersection)%>%lapply(function(x)x$peakID)
    output=list(peakIDs,peak_intersection)
    names(output)=c('IDs','intersection')
    return(output)
}

human_support_peaks <- peak_intersection(human_individual_peaks,human_consensus_peaks)
chimp_support_peaks <- peak_intersection(chimp_individual_peaks,chimp_consensus_peaks)


pdf(paste(outplot_dir,'human-support-consensus-peak.pdf',sep=''),width=10,height=7)
upset(fromList(human_support_peaks$IDs),nsets = 6, order.by = "freq")
dev.off()

pdf(paste(outplot_dir,'chimp-support-consensus-peak.pdf',sep=''),width=10,height=7)
upset(fromList(chimp_support_peaks$IDs),nsets = 6,order.by = "freq")
dev.off()

## export new set of peaks
fwrite(human_consensus_peaks,paste(outfile_dir,'human_consensus_peaks.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
fwrite(chimp_consensus_peaks,paste(outfile_dir,'chimp_consensus_peaks.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

