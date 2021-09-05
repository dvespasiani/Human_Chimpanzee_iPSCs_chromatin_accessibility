## script used to calculate the fraction of CTCF binding sites 
## occurring in our dataset vs the total amount of genome-wide CTCF binding sites
## these info both come from HOMER. for the 1st one I used findMotifs.pl
## for the second one I instead used scanMotifGenomeWide.pl (see scanGenome_for_CTCF.sh)

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

ctcf_dir = './post_processing_analyses/output/homer/homer_output/ctcf/'
ctcf_sequence_dir = paste(ctcf_dir,'ctcf_peak_sequences/',sep='')
peakDir = './post_processing_analyses/output/DA/peaks/'
plot_dir = './post_processing_analyses/output/plots/ctcf'

standard_chr = paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes

##--------------------------------------------------
## read the genomic location of all CTCF motifs
##--------------------------------------------------
genome_wide_ctcf = list.files(paste(ctcf_dir,'genome_wide_ctcf/',sep=''),recursive=F,full.names=T)%>%
    lapply(function(x)fread(x,sep='\t',header=F,col.names=range_keys)[seqnames%in%standard_chr])

names(genome_wide_ctcf) = gsub("\\..*","",
list.files(paste(ctcf_dir,'genome_wide_ctcf/',sep=''),recursive=F,full.names=F)
)

## read peaks that are predicted to have ctcf motifs 
ctcf_matches = list.files(paste(ctcf_dir,'target_sequences/',sep=''),recursive=F,full.names=T)%>%
    lapply(function(x)fread(x,sep='\t',header=T)[,c(1:3)]%>%setnames(old=c('FASTA ID','Sequence'),new=c('peakID','sequence')))

names(ctcf_matches) = gsub("\\..*","",
list.files(paste(ctcf_dir,'target_sequences/',sep=''),recursive=F,full.names=F)
)
 
## get genomic locations of these ctcf peaks
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))

# common_peaks = read_da_peaks('common_da_results.txt',c('significant','non_significant'))
# chimp_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('human')))%>%setnames(old=c(1:4),new=c('peakID',range_keys))
# human_common_peaks = copy(common_peaks)%>%dplyr::select(-c(contains('chimp')))%>%setnames(old=c(1:4),new=c('peakID',range_keys))

# human_all_peaks_hg38 = rbind(human_peaks,human_common_peaks)
# chimp_all_peaks_pantro5 = rbind(chimp_peaks,chimp_common_peaks)

all_peaks = list(chimp_peaks,common_regions,human_peaks)%>%lapply(function(x)x=x[,c(..range_keys,'peakID','significant')])
names(all_peaks) = species_names

# all_peaks = list(chimp_all_peaks_pantro5,human_all_peaks_hg38)%>%lapply(function(x)x=x[,c(..range_keys,'peakID','significant')])
# names(all_peaks) = species_names[-2]

species_specific_peaks = list(chimp_peaks,human_peaks) ## u'll need this below

peaks_with_ctcf = purrr::map2(ctcf_matches,all_peaks,inner_join,by='peakID')

# ## these lines below simply merge all the ctcf peaks into 2 files 1 for chimp and 1 for human 
# ## to compare the set of ctcf peaks against all the genome-wide locations of the CTCF binding sites
# all_peaks_pantro5_hg38 = list(
#     rbind(chimp_common_peaks_pantro5,chimp_peaks),rbind(human_common_peaks_pantro5,human_peaks)
# )
# names(all_peaks_pantro5_hg38) =species_names[-2]

# chimp_specific_and_common_ctcf = rbindlist(peaks_with_ctcf[c(1,2)])[,peakID:=gsub("\\..*","",peakID)][,peakID:=gsub(".*C_", "",peakID)]
# human_specific_and_common_ctcf = rbindlist(peaks_with_ctcf[c(2,3)])[,peakID:=gsub(".*H_", "",peakID)]

# all_peaks_w_ctcf = list(chimp_specific_and_common_ctcf,human_specific_and_common_ctcf)
# names(all_peaks_w_ctcf) =species_names[-2]

##----------
## Some QCs
##----------
## count number peaks with ctcf and get the proportion out of all peaks
numb_peaks_w_ctcf = lapply(peaks_with_ctcf,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
numb_all_peaks = lapply(all_peaks,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

prop_peaks_w_ctcf = purrr::map2(numb_peaks_w_ctcf,numb_all_peaks,`/`)
prop_peaks_w_ctcf = lapply(prop_peaks_w_ctcf,function(x)data.table(prop=x)[,peaktype:='w_ctcf'])
prop_peaks_w_ctcf = Map(mutate,prop_peaks_w_ctcf,species=species_names)
one_minus_prop_peaks_w_ctcf = copy(prop_peaks_w_ctcf)%>%lapply(function(x)x=x[,prop:=1-prop][,peaktype:='no_ctcf'])

final_prop_peaks_ctcf = purrr::map2(prop_peaks_w_ctcf,one_minus_prop_peaks_w_ctcf,rbind)%>%rbindlist()

pdf(paste0(plot_dir,'/prop_peaks_w_ctcf.pdf',sep=''),width = 7, height = 7)
ggplot(final_prop_peaks_ctcf,aes(x=species,y=prop,fill=peaktype))+
geom_bar(stat='identity')+ 
xlab(' ')+ylab('Proportion peaks')+
theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()

## now the same but genome wide
# prop_genwide_ctcf = purrr::map2(numb_peaks_w_ctcf,genome_wide_ctcf,function(x,y)x/nrow(y)*100)  

##-----------------------------------------------------------
## count and plot the number of CTCF binding sites per peak
##-----------------------------------------------------------
## add common/specific information here
numb_ctcf_motif_per_peak = copy(peaks_with_ctcf)%>%lapply(function(z)z[,numb_ctcf_peaks:=.N,by=.(peakID)])
# numb_ctcf_motif_per_peak = purrr::map2(
#     numb_ctcf_motif_per_peak,species_specific_peaks,
#     function(x,y)x[,species:=ifelse(x$peakID %in% y$peakID,'specific','common')]
# )%>%lapply(function(z)z[,numb_ctcf_peaks:=.N,by=.(peakID)])

numb_ctcf_motif_per_peak = Map(mutate,numb_ctcf_motif_per_peak,species=species_names)%>%rbindlist()
numb_ctcf_motif_per_peak=numb_ctcf_motif_per_peak[!(species=='common' & numb_ctcf_peaks>100)] ## single outlier

pdf(paste0(plot_dir,'/numb_ctcf_motifs_per_peak.pdf',sep=''),width = 10, height = 7)
ggplot(numb_ctcf_motif_per_peak,aes(x=species,y=numb_ctcf_peaks,fill=species))+
geom_boxplot(notch=T)+ 
xlab(' ')+ylab('number CTCF motifs per peak')+
stat_compare_means(
    ref.group='common',
    method='wilcox.test'
)+
theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()

##---------------
## export files
##---------------
## write down the files with location of peaks with CTCFs
mapply(write.table,peaks_with_ctcf, file = paste(ctcf_dir,'ctcf_peaks/',names(peaks_with_ctcf),'.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## write down into another file the partial sequence matching CTCF motif for each peak for DeepBind predictions
ctcf_match_sequence = copy(peaks_with_ctcf)%>%lapply(function(x)x=x[,c('sequence','peakID')])
ctcf_match_sequence = Map(mutate,ctcf_match_sequence,file=names(ctcf_match_sequence))%>%rbindlist()

write.table(ctcf_match_sequence$sequence,paste0(ctcf_sequence_dir,'peak_sequences.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)

### do it also for all gen-wide CTCF peaks so that u can compare the affinity predictions
# get_sequence = function(x,genome){
#     peaks = copy(x)%>%as.data.table()
#     peaks = peaks[
#         ,sequence:=as.character(getSeq(genome, seqnames,start, end))
#         ]
#     return(peaks)
# }

# ## write also the sequence for all ctcf binding sites in the pantro5 and hg38 genomes 
# ## and use it to calculate/plot the background distribution of binding affinities
# hg38_ctcf_sequences = copy(genome_wide_ctcf[[1]])%>%get_sequence(Hsapiens)
# pantro5_ctcf_sequences = copy(genome_wide_ctcf[[2]])%>%get_sequence(Ptroglodytes)

# write.table(hg38_ctcf_sequences$sequence,paste0(ctcf_sequence_dir,'hg38_ctcf_sequences.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)
# write.table(pantro5_ctcf_sequences$sequence,paste0(ctcf_sequence_dir,'pantro5_ctcf_sequences.txt',sep=''),sep='\t',row.names=F,col.names=F,quote=F)

## write also a single file with all the info to make it easier to retrieve them in the next script
write.table(ctcf_match_sequence,paste0(ctcf_sequence_dir,'ctcf_match_file.txt',sep=''),sep='\t',row.names=F,col.names=T,quote=F)

