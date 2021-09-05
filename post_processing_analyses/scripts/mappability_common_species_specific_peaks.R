## use this script to plot mappability scores for regions associated with common and species-specific peaks
## remember that a mappability value of 1 at position i indicates that the k-mer in the sequence at position i 
## occurs only once in the sequence with up to e errors. 
## Vice versa, a low mappability value indicates that this k-mer belongs to a repetitive region

library(data.table)
library(magrittr)
library(dplyr)
library(edgeR)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(liftOver)
library(rtracklayer)

options(width=200)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

peakDir = "output/PeakCalling/ConsensusPeaks/orthologous_peaks/"
plot_dir = './post_processing_analyses/output/plots/orth_peaks_liftover/'

##-------------------------
##read mappability results
##-------------------------
read_map = function(species){
    file=dir(paste(species,'/output/Genome_mappability/mappability',sep=''),pattern='.bed',full.names=T)
    map = fread(file,sep='\t',header=T)[
        ,seqnames:=ifelse(seqnames%like% 'chr',seqnames,paste('chr',seqnames,sep=''))
        ]
    return(map)
}
 
hg38_map = read_map('./hg38')
pantro5_map = read_map('./pantro5')

## read peaks
read_peaks = function(genome,file){
    peak = fread(paste(genome,peakDir,file,sep=''),sep='\t',header=T)%>%unique()
    peak = peak[,peakID:=paste(peakID,species,sep='.')]
    setkeyv(peak,range_keys)

    return(peak)
}

human_and_common_peaks = read_peaks('./hg38/','common_and_human_specific_peaks_hg38.bed')
chimp_and_common_peaks = read_peaks('./pantro5/','common_and_chimp_specific_peaks_pantro5.bed')

## get map score for peak regions 
human_peakmap = foverlaps(hg38_map,human_and_common_peaks,type='any')%>%na.omit()
chimp_peakmap = foverlaps(pantro5_map,chimp_and_common_peaks,type='any')%>%na.omit()

## calculate and plot the mean mappability score for each peak
human_meanmap = copy(human_peakmap)[ , .(mean_map = mean(map_score)), by = peakID][,file:='human'][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]
chimp_meanmap = copy(chimp_peakmap)[ , .(mean_map = mean(map_score)), by = peakID][,file:='chimp'][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]

pdf(paste0(plot_dir,'peak_mean_mappability_score.pdf',sep=''),width = 10, height = 5)
ggplot(rbind(human_meanmap,chimp_meanmap), aes(x =factor(round(mean_map,1)),fill = species))+
    geom_bar(position='dodge')+xlab(' ')+ylab('peak mean mappability score')+
    facet_wrap(file~.,ncol=2)+
    theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
dev.off()

