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
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

peakDir = "output/PeakCalling/ConsensusPeaks/"
plot_dir = './post_processing_analyses/output/plots/orth_peaks_liftover/'
out_dir = './post_processing_analyses/output/final_peak_set/'

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

##-------------------------------------
## Read filtered consensus peaks 
##-------------------------------------
read_filtered_consensus_peaks = function(organimsDir){
    peaks = fread(
        paste0(organimsDir,paste(peakDir,'filtered_consensus_peaks.txt',sep='')),sep='\t',header=T)%>%
        setorderv(c(range_keys),1)%>%unique()%>%
        makeGRangesFromDataFrame(keep.extra.columns=T)%>%as.data.table()
    peaks =peaks[,peakID:=paste('peak_',1:nrow(peaks),sep='')]
    return(peaks)
}

human_consensus_peaks = read_filtered_consensus_peaks('./hg38/')[,species:='human']
chimp_consensus_peaks = read_filtered_consensus_peaks('./pantro5/')[,species:='chimp']

## plot distrbution of consensus peak sizes
all_consensus_peaks = rbind(human_consensus_peaks,chimp_consensus_peaks)
all_consensus_peaks = bin_distance(all_consensus_peaks,all_consensus_peaks$width)

binned_width_levels = c(
    '0-49',
    '50-150',
    '151-300',
    '301-450',
    '451-600',
    '601-750',
    '751-900',
    '901-1000',
    '1001-1150',
    '1151-1301',
    '1301-1500',
    '1501-2000',
    '2001-3000',
    '3001-4000',
    '4001-5000',
    '5001-6000',
    '>6000'
)

plot_distr_binned_width = function(x){
    p <- ggplot(x,aes(x=factor(binned_column,levels=binned_width_levels),fill=species))+
    geom_bar()+
    facet_wrap(species~.,scale='free')+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
    return(p)
}

## plot distributon binned width peaks
pdf(paste(plot_dir,'consensus_peak_width_distribution.pdf',sep=''),width=7,height=5)
plot_distr_binned_width(all_consensus_peaks)
dev.off()

##---------------------------------------------------------
## now liftover and then plot again the size distribution
##---------------------------------------------------------

hg38_to_pantro5 <- convert_coord(human_consensus_peaks,'hg38ToPanTro5.over.chain')%>%unique()
hg38_to_pantro5_back <- convert_coord(hg38_to_pantro5,'panTro5ToHg38.over.chain')%>%unique() 

pantro5_to_hg38 <- convert_coord(chimp_consensus_peaks,'panTro5ToHg38.over.chain')[seqnames!='chrM']%>%unique()
pantro5_to_hg38_back <- convert_coord(pantro5_to_hg38,'hg38ToPanTro5.over.chain')[seqnames!='chrM']%>%unique()

all_liftovered_peaks = rbind(hg38_to_pantro5_back,pantro5_to_hg38_back)[,width:=end-start]
all_liftovered_peaks = bin_distance(all_liftovered_peaks,all_liftovered_peaks$width)

## plot the distribution of the new lifted peaks to see if there are significant differences
pdf(paste(plot_dir,'consensus_liftovered_peak_width_distribution.pdf',sep=''),width=7,height=5)
plot_distr_binned_width(all_liftovered_peaks)
dev.off()

##-----------------------------
## filter consensus peaks
##-----------------------------
## keep peaks that have orthologous elements 
## and that have exactly same start/end coord of original peaks
filter_orth_peaks =function(consensus,lifted){
    consensus_peak = copy(consensus)[,width:=end-start]
    lifted_peak = copy(lifted)[,width:=end-start]

    filtered_peaks = consensus_peak[
        lifted_peak,on=c(range_keys,'peakID','species','width'),nomatch=0
        ][
            width>=150 & width<=2000 ## finally keep peaks between 150-2kb
            ][
                ,c('width','strand','i.strand'):=NULL
                ]
    filtered_peaks = filtered_peaks[,peakID:=paste('peak_',1:nrow(filtered_peaks),sep='')] 
    setkeyv(filtered_peaks,range_keys)
    return(filtered_peaks)
}

chimp_filt_orth_peaks = filter_orth_peaks(chimp_consensus_peaks,pantro5_to_hg38_back)
human_filt_orth_peaks = filter_orth_peaks(human_consensus_peaks,hg38_to_pantro5_back)

##-----------------------
## mappability scores
##-----------------------
## read map scores
read_map = function(species){
    file=dir(paste(species,'/output/Genome_mappability/mappability',sep=''),pattern='.bed',full.names=T)
    map = fread(file,sep='\t',header=T)[
        ,seqnames:=ifelse(seqnames%like% 'chr',seqnames,paste('chr',seqnames,sep=''))
        ]
    return(map)
}
 
hg38_map = read_map('./hg38')
pantro5_map = read_map('./pantro5')

## get map score for peak regions 
human_peakmap = foverlaps(hg38_map,human_filt_orth_peaks,type='any')%>%na.omit()
chimp_peakmap = foverlaps(pantro5_map,chimp_filt_orth_peaks,type='any')%>%na.omit()

## calculate and plot the mean mappability score for each peak
human_meanmap = copy(human_peakmap)[,peakID:=paste(peakID,species,sep='.')][ , .(mean_map = mean(map_score)), by = peakID][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]
chimp_meanmap = copy(chimp_peakmap)[,peakID:=paste(peakID,species,sep='.')][ , .(mean_map = mean(map_score)), by = peakID][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]

mean_map_scores=rbind(human_meanmap,chimp_meanmap)[,mean_map:=round(mean_map,1)]

pdf(paste0(plot_dir,'peak_mean_mappability_score.pdf',sep=''),width = 10, height = 5)
ggplot(mean_map_scores, aes(x=factor(mean_map),fill = species))+
    geom_bar(position='dodge')+xlab(' ')+ylab('peak mean mappability score')+
    facet_wrap(species~.,ncol=2)+
    theme(
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
dev.off()

## only retain peaks with high mappability scores (i.e. 1)
peaks_to_keep =copy(mean_map_scores)[mean_map==1]%>%split(by='species')

human_peaks_to_keep=copy(peaks_to_keep[[1]])
chimp_peaks_to_keep=copy(peaks_to_keep[[2]])

human_peak_set = human_filt_orth_peaks[peakID%in%human_peaks_to_keep$peakID]
chimp_peak_set = chimp_filt_orth_peaks[peakID%in%chimp_peaks_to_keep$peakID]

## write peaks
## for the moment avoid wrting headers
write.table(human_peak_set,paste(out_dir,'human_peaks_hg38.bed',sep=''),sep='\t',col.names=F,quote=F,row.names=F)
write.table(chimp_peak_set,paste(out_dir,'chimp_peaks_pantro5.bed',sep=''),sep='\t',col.names=F,quote=F,row.names=F)

## convert chimp orth peaks to hg38 coords
chimp_peak_set_hg38 = convert_coord(chimp_peak_set,'panTro5ToHg38.over.chain')
## tmp file 
write.table(chimp_peak_set_hg38,paste(out_dir,'chimp_peaks_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)


## check peaks that are in common btwn human/chimp 
## use bedtools
get_common_peaks = "bedtools intersect -f 0.5 -r -a post_processing_analyses/output/final_peak_set/human_peaks_hg38.bed -b post_processing_analyses/output/final_peak_set/chimp_peaks_hg38.bed"
common_peaks <- fread(cmd=get_common_peaks)

common_peaks_hg38coord = foverlaps(human_peak_set,chimp_peak_set_hg38,type='any')%>%na.omit()%>%
    setnames(old=c('start','end','i.start','i.end','peakID','i.peakID'),
            new=c('chimp_start','chimp_end','human_start','human_end','chimp_peakID','human_peakID'))%>%
            dplyr::select(-c(contains('species'))
)

##-------------------
## common peaks QCs
##-------------------
## check how many peaks are in common between species 
## check how far start common peaks are from each other
common_peaks_qc = copy(common_peaks_hg38coord)[
    ,start_dist:=abs(chimp_start-human_start)
]
common_peaks_qc = bin_distance(common_peaks_qc,common_peaks_qc$start_dist)

## plot distribution start distances between common peaks
pdf(paste(plot_dir,'distance_btwn_common_peaks.pdf',sep=''),width=7,height=7)
ggplot(common_peaks_qc,aes(x=factor(binned_column,levels=binned_width_levels)))+
    geom_bar()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()

## keep as common peaks all those peaks within 300bp of distance from each other
final_set_common_peaks = copy(common_peaks_qc)[start_dist<=300][,c('start_dist','binned_column'):=NULL]

## label common and species-specific peaks
chimp_peaks_pantro5 = copy(chimp_peak_set)[
    ,species:=ifelse(peakID %in% final_set_common_peaks$chimp_peakID,'common','chimp')
]
human_peaks_hg38 = copy(human_peak_set)[
    ,species:=ifelse(peakID %in% final_set_common_peaks$human_peakID,'common','human')
]

## keep common peaks where both peaks have mappability scores >0.5
common_peaks_hg38coord = common_peaks_hg38coord[
    ,chimp_map_score:=ifelse(chimp_peakID%in%chimp_peaks_to_keep$peakID,'high','low')
    ][
        ,human_map_score:=ifelse(human_peakID%in%human_peaks_to_keep$peakID,'high','low')
        ][
           (chimp_map_score=='high' & human_map_score=='high')
           ][
               ,c('chimp_map_score','human_map_score'):=NULL
]


## write files (these will be input of the DA analysis)
write.table(human_peak_set,paste(out_dir,'common_and_human_specific_peaks_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
write.table(chimp_peak_set,paste(out_dir,'common_and_chimp_specific_peaks_pantro5.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
write.table(common_peaks_hg38coord,paste(out_dir,'common_peaks_hg38coord.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

## For all other downstream analyses in which u compare the common vs species-specific peaks save also a file with the common regions (not peaks)
## this file is the same of the common_peaks_hg38coord file but with ranges being regions not the peaks. 
## this spares me from repeating the code below over and over again across scripts

## common regions
common_regions_hg38 = copy(common_peaks_hg38coord)[
    , max_chimp_end:= max(chimp_end), by=.(chimp_peakID)
    ][
        , min_chimp_start:= min(chimp_start), by=.(chimp_peakID)
        ][
            ,region_start:=ifelse(human_start<min_chimp_start,human_start,min_chimp_start)
            ][
                ,region_end:=ifelse(human_end>max_chimp_end,human_end,max_chimp_end)
                ][
                    ,c('seqnames','region_start','region_end','chimp_peakID','human_peakID')
]
colnames(common_regions_hg38)[1:3] = range_keys

## plot width distribution of these regions
common_regions_hg38_qc = copy(common_regions_hg38)[,width:=end-start]
common_regions_hg38_qc = bin_distance(common_regions_hg38_qc,common_regions_hg38_qc$width)

## plot distribution start distances between common peaks
pdf(paste(plot_dir,'common_regions_size_distribution.pdf',sep=''),width=7,height=7)
ggplot(common_regions_hg38_qc,aes(x=factor(binned_column,levels=binned_width_levels)))+
    geom_bar()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()

## write file
write.table(common_regions_hg38,paste(out_dir,'common_regions_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

