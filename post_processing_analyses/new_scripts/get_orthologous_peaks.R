## Script used to get orthologuous regions for consensus peaks 
## from the consensus peaks use liftover to convert the genomic coordinates
## then keep liftOvered peaks based on width thresholds 
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(dplyr)
library(liftOver)
library(rtracklayer)

options(width=200)

## directories
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

peak_outdir <- 'output/final_peak_set/'
outplot_dir <- create_dir(plot_dir,'orthologous_peaks')
peak_outdir <- create_dir(outdir,'final_peak_set')

## read species peaks 
peak_file <- list.files(path=paste('../',genome,sep=''),recursive=T,full.names= T,pattern="^merged.*narrowPeak$")
peaks <- fread(peak_file,sep='\t',header=F,select=c(1:3),col.names=c(range_keys))%>%makeGRangesFromDataFrame()%>%reduce(min.gapwidth=50L)%>%as.data.table()
peaks <- peaks[,peakID:=paste('peak_',1:nrow(peaks),sep='')][,c('width','strand'):=NULL]

## plot distribution peak sizes
plot_peak_sizes <- function(peaks){
    df <- copy(peaks)[,width:=end-start]
    plot <- ggplot(df,aes(x=width))+
    geom_histogram(binwidth=150,color='black')
    return(plot)
}

pdf(paste(outplot_dir,"distribution_peak_sizes_before_liftover.pdf",sep=''),width=10,height=7)
plot_peak_sizes(peaks)
dev.off()


## liftover back and forth peaks to only keep orthologous ones
chain_files <- list('hg38ToPanTro5.over.chain','panTro5ToHg38.over.chain')

lifted_forth <- convert_coord(peaks=peaks,chain_file=chain_files[[1]])

pdf(paste(outplot_dir,"distribution_peak_sizes_after_liftover.pdf",sep=''),width=10,height=7)
plot_peak_sizes(lifted_forth)
dev.off()

## liftover back to original coords
liftover_back <- convert_coord(lifted_forth,chain_files[[2]])
liftover_back <- liftover_back[,width:=end-start]

##-----------------------------
## keep orthologous peaks 
##-----------------------------
## keep peaks that have orthologous elements 
## and that overlap for 80% width with original peaks
peaks <- peaks[,width:=end-start]

orthologous_peaks <- peaks[
        liftover_back,on=c('peakID'),nomatch=0
        ][
            ,keep:=ifelse(round(i.width/width,2)>=0.5,'keep','no')
            ][ 
                keep=='keep'
                ][
                    ,c(..range_keys,'peakID','width')
]%>%setorderv(range_keys,1)

##-----------------
## export results
##-----------------
filename <- paste(peak_outdir,genome,'_all_orthologous_peaks.txt',sep='')
write.table(orthologous_peaks, file = filename,col.names = T, row.names = F, sep = "\t", quote = F)
