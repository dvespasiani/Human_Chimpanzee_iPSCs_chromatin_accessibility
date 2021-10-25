## script used to look for overlaps between peaks and TADs (eres et al 2019)
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

go_dir <- './output/GO/'
tads_dir <-  './output/TADs/'
target_genes_dir <- './output/files/GO_enrich'
peakDir <- './output/DA/peaks/'
tmp_files_dir <- './output/temp_files/'

outplot_dir <- create_dir(plot_dir,'TADs')

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC','peak_species'))
setkeyv(da_results,range_keys)
da_results <- da_results%>%setorderv('peak_species',1)

all_peaks <- list(da_results[peak_species!='human'],da_results[peak_species!='chimp'])

##------------
## read TADs
##------------
orth_tads <- read_tads('human_chimp_orth_tads.txt')
chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL][species=='Chimp']
human_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

chimp_tads_hg38 <- copy(orth_tads)[species=='Chimp']
colnames(chimp_tads_hg38)[4:6] = range_keys
chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'species'),nomatch=0][
  ,c(range_keys):=NULL
]%>%setnames(old=c(1:3),new=range_keys)

common_tads_hg38 <- copy(human_tads_hg38)[species=='Common']

all_tads = list(rbind(chimp_tads_hg38,common_tads_hg38),human_tads_hg38)# this has both human and common peaks

all_tads <- lapply(all_tads,function(x){
  x <- setorderv(x,range_keys,1)
  x<-x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]
  setkeyv(x,range_keys)
  return(x)
  }
) 
names(all_tads) = species_names[-2]

##---------------------------------
## overlaps of peaks within TADs
##---------------------------------
peak_tads_overlap <- purrr::map2(
  all_peaks,all_tads,
  function(x,y)foverlaps(x,y,type='within')%>%
  na.omit()%>%
  setnames(
    old=c(paste('i.',range_keys[-1],sep='')),
    new=c(paste('peak_',range_keys[-1],sep=''))
    )
)
names(peak_tads_overlap) = species_names[-2]


## some numbers and %
lapply(peak_tads_overlap,function(x)x[,c('seqnames','peak_start','peak_end')]%>%unique()%>%nrow())
# $chimp
# [1] 80040
# $human
# [1] 84957

fraction_peaks_overalp = purrr::map2(
  peak_tads_overlap,all_peaks,
  function(x,y)
  x[,c(1,6,7)]%>%unique()%>%nrow()
  /
  y[,c(..range_keys)]%>%unique()%>%nrow() *100
)
# $chimp 75.0
# $human 82.0


##-------------------------------------------
## now look whether peaks and target genes 
## are in the same TAD
##-------------------------------------------
## read target genes
target_genes <- list.files(target_genes_dir,full.names= T,recursive=F)%>%lapply(
  function(x){
    x<-fread(x,sep='\t',header=T)[,c('peakID','gene','distTSS','peak_species')]
  }
)
all_target_genes <- list(rbindlist(target_genes[c(1:2)]),rbindlist(target_genes[c(2:3)]))

# target_genes = list.files(paste(go_dir,'target_genes',sep=''),full.names=T,recursive=F)%>%
# lapply(
#   function(x)fread(x,sep='\t',header=T)[,c('peakID','gene','distTSS')])
# names(target_genes) = gsub("\\_.*","",list.files(paste(go_dir,'target_genes',sep=''),full.names=F,recursive=F))
# target_genes = Map(mutate,target_genes,target_gene_species=names(target_genes))

# combined_target_genes = list(
#   rbindlist(target_genes[c(1:2)]),
#   rbindlist(target_genes[c(2:3)])
# )
# names(combined_target_genes) = species_names[-2]

peak_in_tads_w_genes <- purrr::map2(
  all_peaks,all_target_genes,
  function(p,g){
    x<-p[g,on=c('peakID','peak_species'),nomatch=0]
    return(x)
  }
)

## Not all peaks have target genes associated (~ 1%)
## you can use these to calculate the %
# purrr::map2(
#   lapply(peak_tads_overlap,function(x)x%>%dplyr::select('peakID')%>%unique()%>%nrow()),
#   lapply(peak_in_tads_w_genes,function(x)x%>%dplyr::select('peakID')%>%unique()%>%nrow()),
#   `/`
# )

## check if peaks and target genes are in the same TAD
peak_in_tads_w_genes <- 
lapply(peak_in_tads_w_genes,function(x)
  x=x[
      ,gene_start := ifelse(distTSS<0,start+distTSS,end+distTSS)
      ][
        ,distTAD := gene_start-start
        ][
          ,.SD[which.min(abs(distTAD))], by=.(peakID)
          ][
            ,peak_genes_same_tad := ifelse(distTAD>=0,'y','n')
            ]
)

## count numb peaks and target genes are in same tad + % of all peaks in tads
qc_peaks_w_genes_in_tads = copy(peak_in_tads_w_genes)%>%
lapply(
  function(x){
    x<-x%>%dplyr::select(c('peakID','peak_species','peak_genes_same_tad'))%>%unique()
    x<-x[
      ,numbpeaks_w_genes_in_out_tad:=.N,by=.(peak_genes_same_tad)
      ][
        ,numbpeaks_in_tad:=.N
        ][
          ,prop_peaks_w_genes_in_out_tad:=numbpeaks_w_genes_in_out_tad/numbpeaks_in_tad
          ]
  }
)
names(qc_peaks_w_genes_in_tads) = species_names[-2]
qc_peaks_w_genes_in_tads <- Map(mutate,qc_peaks_w_genes_in_tads,file=names(qc_peaks_w_genes_in_tads))%>%rbindlist()


pdf(paste0(outplot_dir,'prop_peaks_w_gene_in_out_tad.pdf',sep=''),width = 7, height = 7)
ggplot(
  qc_peaks_w_genes_in_tads[,c('prop_peaks_w_genes_in_out_tad','peak_genes_same_tad','file')]%>%unique(),
  aes(x=file,y=prop_peaks_w_genes_in_out_tad,fill=peak_genes_same_tad))+
geom_bar(stat='identity')+
xlab(' ')+ylab('Proportion peaks with gene')+
theme(
  legend.key = element_rect(fill = "white", colour = "black"),
  axis.line = element_blank()
)
dev.off()

# ##----------------------------------------------------------
# ## Now keep only peaks with target genes in the same TAD
# ##----------------------------------------------------------
# peak_w_genes_same_tad = copy(peak_in_tads_w_genes)%>%
# lapply(
#   function(x)
#   x=x[
#     peak_genes_same_tad=='y'
#     ][
#       ,c('seqnames','peak_start','peak_end','peakID','significant','peak_species','tadID','gene','distTSS')
#       ]%>%unique()
#       )%>%lapply(
#             function(y)y=y[
#               ,numbpeaks_per_gene_same_tad:=.N,by=.(gene,tadID)
#               ]%>%setnames(old=c(1:3),new=range_keys)%>%unique()
# )
# numb_peaks_w_genes_same_tad =lapply(peak_w_genes_same_tad,function(x)x[,c(1:3)]%>%unique()%>%nrow())
# # $chimp 18189
# # $human 21624

# numb_genes_w_peaks_same_tad=lapply(peak_w_genes_same_tad,function(x)x[,c('gene')]%>%unique()%>%nrow())
# # $chimp 8619
# # $human 10107

# peak_w_genes_same_tad = Map(mutate,peak_w_genes_same_tad,file=names(peak_w_genes_same_tad))%>%rbindlist()

# ## write down list of peaks + genes in the same TAD for the DE/DA analysis 
# ## also include the field with numb peaks per gene so that 
# ## later u can filter for those genes w multiple peaks associated for Gviz 
# list_peaks_w_genes_same_tad = copy(peak_w_genes_same_tad)[,c('numb_peaks','numb_peaks_per_regulation','prop_regulation'):=NULL]
# write.table(list_peaks_w_genes_same_tad,paste(tmp_files_dir,'list_peaks_w_genes_same_tad.txt',sep=''),sep='\t',col.names=T, quote=F,row.names=F)

