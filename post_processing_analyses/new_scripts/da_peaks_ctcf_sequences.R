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

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

ctcf_dir <- './output/homer/homer_output/ctcf/'
peakDir <- './output/DA/peaks/'
outplot_dir <- create_dir(plot_dir,'ctcf')
ctcf_outdir <- create_dir(ctcf_dir,'ctcf_peaks')

##--------------------------------------------------
## read the genomic location of all CTCF motifs
##--------------------------------------------------
genome_wide_ctcf_file <- dir(paste(ctcf_dir,'genome_wide_ctcf',sep=''),recursive=F,full.names=T,pattern=genome)
genome_wide_ctcf <- fread(genome_wide_ctcf,sep='\t',header=F,col.names=range_keys)[seqnames%in%standard_chr]


## read peaks that are predicted to have ctcf motifs 
ctcf_matches <- list.files(paste(ctcf_dir,'target_sequences/',sep=''),recursive=F,full.names=T,pattern=genome)%>%
    lapply(function(x)fread(x,sep='\t',header=T)[,c(1:3)]%>%setnames(old=c('FASTA ID','Sequence'),new=c('peakID','sequence')))

names(ctcf_matches) = species_names

## read peaks to get genomic locations of ctcf bs
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
all_peaks <- fread(da_file,sep='\t',header=T,select=c(range_keys,'logFC','DA','peakID','peak_species'))
setkeyv(all_peaks,range_keys)
all_peaks <- all_peaks%>%setorderv('peak_species',1)%>%split(by='peak_species')
 
peaks_with_ctcf = purrr::map2(ctcf_matches,all_peaks,function(x,y)x[y,on='peakID',nomatch=0])

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

pdf(paste0(outplot_dir,'prop_peaks_w_ctcf.pdf',sep=''),width = 7, height = 7)
ggplot(final_prop_peaks_ctcf,aes(x=species,y=prop,fill=peaktype))+
geom_bar(stat='identity')+ 
xlab(' ')+ylab('Proportion peaks')+
theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()

##---------------
## export files
##---------------
mapply(write.table,peaks_with_ctcf, file = paste(ctcf_outdir,names(peaks_with_ctcf),'_ctcf_peaks','.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

