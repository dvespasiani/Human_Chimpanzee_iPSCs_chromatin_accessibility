## input: the human and chimpanzee peaks called on the merged bams in the respective coordinates 
## output: the set of peaks that have orthologous regions in the other species
## qc: check how many individuals support each consensus peak
## use liftover back and forth

library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(liftOver)
library(rtracklayer)
library(UpSetR)
library(ComplexHeatmap)
library(viridis)
library(GenomicRanges)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility')

scripts_dir <- './post_processing_analyses/scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

peakDir = "output/PeakCalling/Files/"
outplot_dir <- create_dir(paste('./post_processing_analyses/',plot_dir,sep=''),'consensus_orthologous_peaks')

outtable_dir <- create_dir('post_processing_analyses/output/tables/','consensusPeaks')

human_chimp_col <- c('#14213d','#fca311')
names(human_chimp_col) = c('chimp','human')

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

fwrite(human_consensus_peaks,paste(outtable_dir,'human_consensus_peaks.bed.gz',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
fwrite(chimp_consensus_peaks,paste(outtable_dir,'chimp_consensus_peaks.bed.gz',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

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

# human_support_peaks <- peak_intersection(human_individual_peaks,human_consensus_peaks)
# chimp_support_peaks <- peak_intersection(chimp_individual_peaks,chimp_consensus_peaks)


# pdf(paste(outplot_dir,'human_support_consensus_peak.pdf',sep=''),width=10,height=7)
# upset(fromList(human_support_peaks$IDs),nsets = 6, order.by = "freq")
# dev.off()

# pdf(paste(outplot_dir,'chimp_support_consensus_peak.pdf',sep=''),width=10,height=7)
# upset(fromList(chimp_support_peaks$IDs),nsets = 6,order.by = "freq")
# dev.off()

# ##------------------------------------
# ## Filter peaks with low mappability  
# ##------------------------------------
# read_map = function(genome){
#     file=dir(paste(genome,'/genome_mappability/mappability',sep=''),pattern='*genmap.bedgraph',full.names=T)
#     map = fread(file,sep='\t',header=F)%>%setnames(c(range_keys,'mappability_score'))
#     map <- map[seqnames %in% standard_chr]

#     return(map)
# }
 
# # hg38_map = read_map('./hg38')
# pantro5_mappability = read_map('./panTro5')




# ## read peaks
# read_peaks = function(genome,file){
#     peak = fread(paste(genome,peakDir,file,sep=''),sep='\t',header=T)%>%unique()
#     peak = peak[,peakID:=paste(peakID,species,sep='.')]
#     setkeyv(peak,range_keys)

#     return(peak)
# }

# human_and_common_peaks = read_peaks('./hg38/','common_and_human_specific_peaks_hg38.bed')
# chimp_and_common_peaks = read_peaks('./pantro5/','common_and_chimp_specific_peaks_pantro5.bed')

# ## get map score for peak regions 
# human_peakmap = foverlaps(hg38_map,human_and_common_peaks,type='any')%>%na.omit()
# chimp_peakmap = foverlaps(pantro5_map,chimp_and_common_peaks,type='any')%>%na.omit()

# ## calculate and plot the mean mappability score for each peak
# human_meanmap = copy(human_peakmap)[ , .(mean_map = mean(map_score)), by = peakID][,file:='human'][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]
# chimp_meanmap = copy(chimp_peakmap)[ , .(mean_map = mean(map_score)), by = peakID][,file:='chimp'][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]

# pdf(paste0(plot_dir,'peak_mean_mappability_score.pdf',sep=''),width = 10, height = 5)
# ggplot(rbind(human_meanmap,chimp_meanmap), aes(x =factor(round(mean_map,1)),fill = species))+
#     geom_bar(position='dodge')+xlab(' ')+ylab('peak mean mappability score')+
#     facet_wrap(file~.,ncol=2)+
#     theme(
#     axis.line = element_blank(),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
#     )
# dev.off()


## Use liftover to get orthologous regions
hg38_to_pantro5 <- convert_coord(human_consensus_peaks,'hg38ToPanTro5.over.chain')
hg38_to_pantro5_back <- convert_coord(hg38_to_pantro5,'panTro5ToHg38.over.chain')

pantro5_to_hg38 <- convert_coord(chimp_consensus_peaks,'panTro5ToHg38.over.chain')
pantro5_to_hg38_back <- convert_coord(pantro5_to_hg38,'hg38ToPanTro5.over.chain')

## keep peaks that have orthologous elements 
## and that have exactly same start/end coord of original peaks
filter_orth_peaks =function(consensus,lifted){
    consensus_peak = copy(consensus)
    final_peak = copy(lifted)[
        consensus_peak,on=c(range_keys,'peakID'),nomatch=0
        ][
            ,width:=end-start
            ][
                width>=150 & width<=2000
                ]
    setkeyv(final_peak,range_keys)
    return(final_peak)
}

chimp_filt_orth_peaks <- filter_orth_peaks(chimp_consensus_peaks,pantro5_to_hg38_back)
human_filt_orth_peaks <- filter_orth_peaks(human_consensus_peaks,hg38_to_pantro5_back)

human_support_orthologous_peaks <- peak_intersection(human_individual_peaks,human_filt_orth_peaks)
chimp_support_orthologous_peaks <- peak_intersection(chimp_individual_peaks,chimp_filt_orth_peaks)


pdf(paste(outplot_dir,'human_support_orthologous_consensus_peak.pdf',sep=''),width=10,height=7)
upset(fromList(human_support_orthologous_peaks$IDs),nsets = 6, order.by = "freq")
dev.off()

pdf(paste(outplot_dir,'chimp_support_orthologous_consensus_peak.pdf',sep=''),width=10,height=7)
upset(fromList(chimp_support_orthologous_peaks$IDs),nsets = 6,order.by = "freq")
dev.off()

## keep peaks that are found in at least 2 individuals
filter_peaks <- function(df,consensus,numb.samples){
    keep <- copy(df)%>%rbindlist()
    keep <- keep[,c('peakID','sample')]%>%unique()
    keep <- keep[,numbsamples:=.N,by=.(peakID)][numbsamples>=numb.samples]%>%dplyr::pull('peakID')%>%unique()

    return_peaks <- copy(consensus)[peakID %in% keep]
    return(return_peaks)
}

human_filtered_peaks <- filter_peaks(human_support_orthologous_peaks[[2]],human_consensus_peaks,2)
chimp_filtered_peaks <- filter_peaks(chimp_support_orthologous_peaks[[2]],chimp_consensus_peaks,2)

## export new set of peaks
outfile_dir <- './post_processing_analyses/output/files/orthologous_consensus_filtered_peaks/'

write.table(human_filtered_peaks,paste(outfile_dir,'human_filtered_orth_consensus_peaks.txt',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
write.table(chimp_filtered_peaks,paste(outfile_dir,'chimp_filtered_orth_consensus_peaks.txt',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
