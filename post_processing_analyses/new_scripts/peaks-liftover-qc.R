## Script used to get orthologuous regions for consensus peaks 
## from the consensus peaks use liftover to convert the genomic coordinates
## then keep liftOvered peaks based on width thresholds 
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(dplyr)

options(width=200)

## directories
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

indir <- paste(outdir,'files/consensusPeaks/',sep='')
outplot_dir <- create_dir(plot_dir,'orthologous_consensus_peaks')

species = c('chimp','human')

consensusPeaks <- list.files(indir,recursive=F,full.names=T,pattern='*consensus_peaks.bed')%>%lapply(function(x)fread(x,sep='\t'))
names(consensusPeaks) = species

liftBackForthConsensusPeaks <- list.files(indir,recursive=F,full.names=T,pattern='*liftBack.bed')%>%lapply(function(x)fread(x,sep='\t'))
names(liftBackForthConsensusPeaks) = species

## calculate proportion mapped peaks per species and per support
consensusPeaks <- lapply(consensusPeaks,function(x){
    x <- x[,numbPeaks := .N,by=.(support)][,width := end-start]
})

liftBackForthConsensusPeaks <- lapply(liftBackForthConsensusPeaks,function(x){
    x <- x[,numbPeaks := .N,by=.(support)][,width := end-start]
})


propMappedSupportPeaks <- purrr::map2(consensusPeaks,liftBackForthConsensusPeaks,function(c,m){
    c <- copy(c)[,c('support','numbPeaks','peakID')]%>%unique()
    m <- copy(m)[,c('support','numbPeaks','peakID')]%>%unique()
    return <- c[m,on=c('support','peakID'),nomatch=0][,c('support','peakID','numbPeaks','i.numbPeaks')][
        ,prop := i.numbPeaks/numbPeaks
    ]
    return(return)
})

propMappedSupportPeaks <- Map(mutate,propMappedSupportPeaks,species=names(propMappedSupportPeaks))%>%rbindlist()
propMappedSupportPeaks <- propMappedSupportPeaks[,peakID := NULL]%>%unique()


## plot histogram number of peaks and proportion of retained mapped peaks by support
p1 <- ggplot(propMappedSupportPeaks[support>0],aes(x=reorder(support,support),y=i.numbPeaks))+
geom_histogram(stat='identity')+xlab('individuals supporting the peak call') + ylab('number of peaks')+
facet_wrap(species~.,ncol=2)+
theme_classic()+
theme()

p2 <- ggplot(propMappedSupportPeaks[support>0],aes(x=reorder(support,support),y=prop))+
geom_histogram(stat='identity')+xlab('individuals supporting the peak call') + ylab('proportion mapped peaks')+
facet_wrap(species~.,ncol=2)+
theme_classic()+
theme()

pdf(paste(outplot_dir,'liftoverQC-prop-peaks-retained.pdf',sep=''),width=7,height = 7)
list(p1,p2)
dev.off()


## plot distribution differences peak sizes before/after liftover
## and remove peaks that after double liftover whose width is > x% diverse where x is visually determined
peakSizes <- purrr::map2(consensusPeaks,liftBackForthConsensusPeaks,function(c,m){
    c <- copy(c)[,c('support','width','peakID')]%>%unique()
    m <- copy(m)[,c('support','width','peakID')]%>%unique()
    return <- c[m,on=c('support','peakID'),nomatch=0][
        ,diffPeakSize := i.width/width
    ]
    return(return)
})

peakSizes <- Map(mutate,peakSizes,species=names(peakSizes))%>%rbindlist()

p3 <- ggplot(peakSizes[support>0],aes(x=diffPeakSize,fill=species))+
geom_density(alpha=0.5)+xlab('ratio liftovered peak size/original peak size') + ylab('Density')+
scale_fill_manual(values = da_species_palette[-1])+
theme_classic()+
theme()

pdf(paste(outplot_dir,'liftoverQC-distribution-peak-sizes.pdf',sep=''),width=7,height = 7)
p3
dev.off()
