## use this script to run some QC on the DA peaks obtained from csaw analysis
## particularly, check the numebr of DA peaks in both species
## after liftovering chimp, check how many DA peaks are in common between the species
## etc...
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(liftOver)

options(width=200)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

plot_dir = './post_processing_analyses/output/plots/DA/'
peakDir = './post_processing_analyses/output/DA/peaks/'

## read peaks
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_peaks = read_da_peaks('common_da_results.txt',c('significant','non_significant'))

## function to count number peaks per peaktype (DA/nonDA)
prop_peaks = function(x){
    x=x[
        ,numb_peaktype :=.N,by=.(significant)
        ][
            ,tot_peaks :=.N
            ][
                ,prop_peaktype:=numb_peaktype/tot_peaks
                ]
    return(x)
}


## first plot the proportion of signif/n.s. peaks across all 3 groups
all_species = list(common_peaks,chimp_peaks,human_peaks)
names(all_species) = c('common','chimp','human')

all_species = lapply(all_species,function(x)x=prop_peaks(x))
all_species = Map(mutate,all_species,species=names(all_species))

proportion_peaks = lapply(all_species,function(x)x=x[,c('species','prop_peaktype','significant')]%>%unique())%>%rbindlist()

pdf(paste(plot_dir,'prop_da_nonda_peaks_species.pdf',sep=''),width=7,height=7)
ggplot(proportion_peaks,aes(x=species,y=prop_peaktype,fill=significant))+
    geom_bar(position="stack", stat="identity")
dev.off()

## now convert pantro5 common peaks to hg38 coordinates 
## then get the common region between those peaks

chimp_common_peaks_pantro5 = copy(common_peaks)%>%dplyr::select(c(contains('chimp'),'significant'))
colnames(chimp_common_peaks_pantro5)[1:4] = c('peakID',range_keys)
setkeyv(chimp_common_peaks_pantro5,range_keys)

human_common_peaks_hg38 = copy(common_peaks)%>%dplyr::select(c(contains('human'),'significant'))
colnames(human_common_peaks_hg38)[1:4] = c('peakID',range_keys)
setkeyv(human_common_peaks_hg38,range_keys)


chimp_common_peaks_hg38 = convert_coord(chimp_common_peaks_pantro5,'panTro5ToHg38.over.chain')
setkeyv(chimp_common_peaks_hg38,range_keys)


## calculate and plot the distance between the starts of each common peak
common_regions = foverlaps(human_common_peaks_hg38,chimp_common_peaks_hg38,type='any')[
    ,same_significance:=ifelse(significant==i.significant,'y','n')
    ][
        ,significant:=ifelse(same_significance=='y',significant,'significant')
        ][
            ,start_dist:=abs(start-i.start)
]
common_regions = bin_distance(common_regions,common_regions$start_dist)

binned_width_levels = c(
        '0-49',
        '50-150',
        '151-300',
        '301-450',
        '451-600',
        '601-750',
        '751-900',
        '901-1000',
        '1001-1500',
        '1501-2000',
        '2001-3000',
        '3001-5000'
)

pdf(paste(plot_dir,'distance_btwn_common_species_common_peaks.pdf',sep=''),width=7,height=7)
ggplot(common_regions,aes(x=factor(binned_column,levels=binned_width_levels)))+
    geom_bar()
dev.off()

## now proceed with extraxting the common regions between the common peaks
common_regions = common_regions[
    ,start:=ifelse(start < i.start,i.start,start)
    ][
        ,end:=ifelse(end > i.end,end,i.end)
        ][
            ,peakID:=paste(paste('C',peakID,sep='_'),paste('H',i.peakID,sep='_'),sep='.')
        ][
            ,c(..range_keys,'peakID','significant')
]%>%unique()

##------------------------
## export common regions
##------------------------

write.table(common_regions,paste(peakDir,'common_regions.txt',sep=''),quote=F,col.names = TRUE,row.names=F,sep='\t')
