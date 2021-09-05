library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(liftOver)
library(rtracklayer)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

peakDir = "post_processing_analyses/output/DA/peaks/"
output_dir = './post_processing_analyses/output/homer/homer_input/'

## using peaks after DA analysis (both da and non DA)
arrange_columns_for_homer = function(x,othercols){x=x[
    ,strand:= '+'][
            ,c(..range_keys,'peakID','strand',..othercols)]
}

human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))%>%arrange_columns_for_homer(othercols=c("significant"))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))%>%arrange_columns_for_homer(othercols=c("significant"))
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))%>%arrange_columns_for_homer(othercols=c("significant"))
# common_peaks= read_da_peaks('common_da_results.txt',c('significant','non_significant'))
# chimp_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('human')))%>%setnames(old=c(1:4),new=c('peakID',range_keys)) %>%arrange_columns_for_homer(othercols=c("significant"))
# human_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('chimp')))%>%setnames(old=c(1:4),new=c('peakID',range_keys)) %>%arrange_columns_for_homer(othercols=c("significant"))

all_peaks = list(chimp_peaks,common_regions,human_peaks)
names(all_peaks) = c('chimp_specific_pantro5','common_regions_hg38','human_specific_hg38')
mapply(write.table, all_peaks, file = paste0(output_dir,names(all_peaks),sep='.txt'),sep='\t',col.names=T,row.names=F,quote=F)



