## use this script to qc the consensus peak call
## in particular see how many individuals support each peak call
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(dplyr)
library(liftOver)
library(rtracklayer)
library(UpSetR)

options(width=200)

## directories
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

peakDir = "output/PeakCalling/ConsensusPeaks/"
plot_dir = './post_processing_analyses/output/plots/orth_peaks_liftover/'

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

##----------------------
## read consensus peak
##----------------------
read_consensus_peaks = function(organimsDir){
    consensus_peak = fread(
        paste0(organimsDir,paste(peakDir,'consensus_peak_05.bed.gz',sep='')),
        sep='\t',header=T,select = c(1:3),col.names = range_keys
        )[
            ,width:=end-start
        ]%>%setorderv(c(range_keys),1)%>%unique()%>%
        makeGRangesFromDataFrame(keep.extra.columns=T)%>%
        reduce(min.gapwidth=20L)%>%as.data.table()
    consensus_peak = consensus_peak[
        ,peakID:=paste('peak_',1:nrow(consensus_peak),sep='')
        ][
            !seqnames %like% '_'
        ]## remove peaks mapping non standard chrs
    setkeyv(consensus_peak,range_keys)
    return(consensus_peak)
}

human_consensus_peaks = read_consensus_peaks('./hg38/')[,species:='human']
chimp_consensus_peaks = read_consensus_peaks('./pantro5/')[,species:='chimp']

## read all peaks for all individuals
read_peaks = function(species){
    file_path=dir(paste(species,'/output/PeakCalling/Files/',sep=''),pattern='*default_peaks_filtered.narrowPeak*',full.names=T)
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

chimp_peaks = read_peaks('./pantro5')
human_peaks = read_peaks('./hg38')

## get list of peaks
get_peak_list=function(peaks,consensus){
    peak_intersection=lapply(peaks,function(x)x=foverlaps(x,consensus,type='any')%>%na.omit())
    peakIDs=copy(peak_intersection)%>%lapply(function(x)x$peakID)
    output=list(peakIDs,peak_intersection)
    names(output)=c('IDs','intersection')
    return(output)
}

human_support_peaks=get_peak_list(human_peaks[c(1:6)],human_consensus_peaks)
chimp_support_peaks=get_peak_list(chimp_peaks[c(1:6)],chimp_consensus_peaks)

pdf(paste(plot_dir,'human_support_consensus_peak.pdf',sep=''),width=10,height=7)
upset(fromList(human_support_peaks$IDs),nsets = 6, order.by = "freq")
dev.off()

pdf(paste(plot_dir,'chimp_support_consensus_peak.pdf',sep=''),width=10,height=7)
upset(fromList(chimp_support_peaks$IDs),nsets = 6,order.by = "freq")
dev.off()

## remove peaks found in just one individual
## take the list of peaks in common with the merged file of all but i-th sample

remove_unique_peaks=function(peaks,sample){

    peaks_all_samples_but_ith=copy(peaks[[1]][c(-sample)])%>%unlist()
    peaks_ith_sample=copy(peaks[[1]][[sample]])
    ids_to_remove=peaks_ith_sample[!peaks_ith_sample %in% peaks_all_samples_but_ith]
    
    consensus_peaks_no_ith_sample = copy(peaks[[2]][[sample]])[!peakID %in% ids_to_remove]
    
    return(consensus_peaks_no_ith_sample)
}

sample_numbers=as.list(1:6)

human_filtered_consensus_peaks=copy(human_support_peaks)
human_filtered_consensus_peaks = lapply(sample_numbers,function(x)remove_unique_peaks(human_filtered_consensus_peaks,x))

chimp_filtered_consensus_peaks=copy(chimp_support_peaks)
chimp_filtered_consensus_peaks = lapply(sample_numbers,function(x)remove_unique_peaks(chimp_filtered_consensus_peaks,x))

## final set of consensus peaks
chimp_final_consensus_peaks=copy(chimp_filtered_consensus_peaks)%>%rbindlist()
chimp_final_consensus_peaks = chimp_final_consensus_peaks[,c(..range_keys)]%>%unique()

human_final_consensus_peaks=copy(human_filtered_consensus_peaks)%>%rbindlist()
human_final_consensus_peaks = human_final_consensus_peaks[,c(..range_keys)]%>%unique()

write.table(chimp_final_consensus_peaks,paste('./pantro5/',peakDir,'filtered_consensus_peaks.txt',sep=''),sep='\t',quote=F,row.names=F,col.names=T)
write.table(human_final_consensus_peaks,paste('./hg38/',peakDir,'filtered_consensus_peaks.txt',sep=''),sep='\t',quote=F,row.names=F,col.names=T)

