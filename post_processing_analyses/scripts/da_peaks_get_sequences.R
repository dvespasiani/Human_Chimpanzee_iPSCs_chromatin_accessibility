## use this script to generate a fasta file containing the sequence
## in the respective species coordinates for each the peak 

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)


options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

peakDir = "post_processing_analyses/output/DA/peaks/"
output_dir = './post_processing_analyses/output/sequences/'

## get DNA sequences 
get_sequence = function(x,genome){
    peaks = copy(x)%>%as.data.table()
    peaks = peaks[
        ,sequence:=as.character(getSeq(genome, seqnames,start, end))
        ]
    return(peaks)
}

## create fasta file
extract_dna_seq = function(x){
    dna_seq=DNAStringSet(x$sequence)
    names(dna_seq) = x$peakID
    return(dna_seq)
}

human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))%>%get_sequence(Hsapiens)%>%extract_dna_seq()
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))%>%get_sequence(Ptroglodytes)%>%extract_dna_seq()
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))%>%get_sequence(Hsapiens)%>%extract_dna_seq()

# common_peaks = read_da_peaks('common_da_results.txt',c('significant','non_significant'))
# chimp_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('human')))%>%setnames(old=c(1:4),new=c('peakID',range_keys))
# human_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('chimp')))%>%setnames(old=c(1:4),new=c('peakID',range_keys))

# all_human_seq = rbind(human_peaks,human_common_peaks)%>%get_sequence(Hsapiens)%>%extract_dna_seq()
# all_chimp_seq = rbind(chimp_peaks,chimp_common_peaks)%>%get_sequence(Ptroglodytes)%>%extract_dna_seq()

all_seq =list(chimp_peaks,common_regions,human_peaks)
names(all_seq) = c("chimp_peaks_pantro5",'common_regions_hg38','human_peaks_hg38')

mapply(writeXStringSet, all_seq, file = paste0(output_dir,names(all_seq),sep='.fa'),append=F,compress=F, compression_level=NA, format="fasta")

