## script to test the enrichment across each chromstate cell type combination
## of common DA peaks (with different regulation) vs species specific DA peaks
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

options(width=150)
set.seed(2022)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')
source('./scripts/utils.R')

chrom_state_dir <- '../data/iPSC_chrom_states_hg38'
outplot_dir <- create_dir(plot_dir,'chromstate_annotation')
out_tabledir <- create_dir(table_dir,'funcAnnot')

## get DA peaks
da_results <- read_da_results('da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC','peaktype')]

## blacklisted regions
blacklist <- fread('../data/ENCODE_blacklisted/hg38_blacklist_v2.bed',sep='\t',header=F,col.names=range_keys)%>%makeGRangesFromDataFrame()

## get tissues 
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

## read random genomic regions
random_peaks <- fread('./output/files/random_regions/random_genomic_regions.txt',sep='\t',header=T)
random_peaks <- random_peaks[,peakID:=paste('peak_',1:nrow(random_peaks),sep='')]

## overlap ranges
annotate_peaks <- function(peaks){
  annotated <- foverlaps(copy(peaks),ipsc_chromstate,type='any')%>%na.omit()
  annotated <- annotated[
    ,overlap:=ifelse(i.start<start,i.end-start,end-i.start),by=.(cell_type)
    ][
        ,.SD[which.max(overlap)], by=.(peakID,cell_type)
        ][
            ,c(range_keys[-1],'overlap'):=NULL
            ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  return(annotated)
}

mypeaks_annotation <- annotate_peaks(da_results)
randompeaks_annotation <- annotate_peaks(random_peaks)

## OR per cell type my peaks vs random set of regions
count_peaks_chromstate <- function(peaks){
  counts <- copy(peaks)
  counts<-counts[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state,cell_type)
      ][
        ,numb_peaks:=.N,by=.(cell_type)
        ][
            ,c('cell_type','chrom_state','numb_peaks_chromstate','numb_peaks')
            ]%>%unique()
  return(counts)
}

mypeaks_counts <- count_peaks_chromstate(mypeaks_annotation)%>%split(by='cell_type')
randompeaks_counts <- count_peaks_chromstate(randompeaks_annotation)%>%split(by='cell_type')

genwide_peaks_or <- purrr::map2(mypeaks_counts,randompeaks_counts,function(x,y){
  or<-calculate_or(x,y,'chrom_state')%>%adjust_pvalues()
  or <- or[,6:=NULL]
  return(or)
})
genwide_peaks_or <- Map(mutate,genwide_peaks_or,cell_type=names(genwide_peaks_or))%>%rbindlist()

plot_or <- function(peaks){
  p <- ggplot(peaks,aes(x=factor(elements,levels=chrom_states),y=log(odds_ratio),fill=elements))+
    geom_violin(trim=T,scale = "width")+
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=0.1)+
    scale_fill_manual(values = chrom_state_colors)+
    geom_hline(yintercept=0,linetype='dashed')+
    xlab('chromatin state') + ylab('log OR')+
    theme_classic()+
    theme(
      legend.position = "bottom",
      axis.text.x =element_blank(),
      axis.ticks.x =element_blank()
    )
  return(p)
}

# pdf(paste(outplot_dir,'chromstate-enrich-mypeaks-vs-randomregion.pdf',sep=''),width=7,height = 7)
# plot_or(genwide_peaks_or)
# dev.off()

## export or results
df <- copy(genwide_peaks_or)%>%setnames(old='elements',new='chrom_state')
fwrite(df,paste(out_tabledir,'mypeaks_randomregion_or_chromstate.txt',sep=''),sep='\t',quote=F,row.names=F,col.names=T)

## or da vs non da 
da_peaks_count <- count_peaks_chromstate(mypeaks_annotation[DA=='da'])%>%split(by='cell_type')
nonda_peaks_count <- count_peaks_chromstate(mypeaks_annotation[DA!='da'])%>%split(by='cell_type')

danonda_peaks_or <- purrr::map2(da_peaks_count,nonda_peaks_count,function(x,y){
  or<-calculate_or(x,y,'chrom_state')%>%adjust_pvalues()
  or <- or[,6:=NULL]
  return(or)
})
danonda_peaks_or <- Map(mutate,danonda_peaks_or,cell_type=names(danonda_peaks_or))%>%rbindlist()

pdf(paste(outplot_dir,'chromstate-enrich-da-vs-nonda.pdf',sep=''),width=7,height = 7)
plot_or(danonda_peaks_or)
dev.off()

## export or results
df2 <- copy(danonda_peaks_or)%>%setnames(old='elements',new='chrom_state')
fwrite(df2,paste(out_tabledir,'da_nonda_peaks_or_chromstate.txt',sep=''),sep='\t',quote=F,row.names=F,col.names=T)

## get independent sets of promoters and calculate
## 1) CpG and GC content
## 2) PhastCons score

promoter_states <- 'Tss|Flnk'

# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# all_promoters <- genes(txdb,single.strand.genes.only=FALSE)%>%as.data.table()
# all_promoters <- all_promoters[!seqnames%like% "M|Un|_"][,end:=start][,start:=start-500][,width:=end-start][,c(..range_keys)][,type:='upstream_txdb']

ipsc_promoters <- copy(ipsc_chromstate)[chrom_state %like% promoter_states][,c(1:3)][,type:='ipsc_promoters']%>%unique()
# random_promoters <- copy(randompeaks_annotation)[chrom_state %like% promoter_states][,c(..range_keys)][,type:='random']%>%unique()
mypromoters <- copy(mypeaks_annotation)[chrom_state %like% promoter_states][,c(..range_keys,'DA')][,type:=DA][,DA:=NULL]%>%unique()

promoters <- rbind(mypromoters,ipsc_promoters)%>%unique()

## CpG & GC
cpg <- copy(promoters)%>%makeGRangesFromDataFrame()%>%Repitools::cpgDensityCalc(organism=Hsapiens)
promoters_cpg <- copy(promoters)[,cpg_content:=cpg][,norm_cpg_content:=cpg/(end-start)]

promoters_gc <- copy(promoters)[
    ,sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
    ][
      !sequence%like%'N' ## do this because u need to remove sequences with Ns 
      ][
        ,gc_content:=letterFrequency(DNAStringSet(sequence),"GC")/width(DNAStringSet(sequence))
        ][
          ,sequence:=NULL
]

comparisons = list(
  c('ipsc_promoters','da'),
  c('ipsc_promoters','non_da'),
  c('da','non_da')
)

## plot function
plot_promoters <- function(x,column_to_plot,ylab){

  df <- copy(x)[,column:=column_to_plot]
  
  type_order <- c('ipsc_promoters','da','non_da')
  
  p <- ggplot(df,aes(x=factor(type,level=type_order),y=column,fill=type))+
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
    scale_fill_manual(values = c('grey',da_palette),labels=c('ipsc_promoters',names(da_palette)))+
  xlab(' ')+ylab(ylab)+
  stat_compare_means(
  method = "wilcox.test",
  comparisons = comparisons,
  size=5
  )+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x =element_blank(),
    axis.ticks.x =element_blank()
  )
  return(p)
}

pdf(paste0(outplot_dir,'peaks-CpG-content.pdf',sep=''),width = 7, height = 7)
plot_promoters(promoters_cpg,promoters_cpg$norm_cpg_content,'Proportion of CpG dinucleotide')
dev.off()

pdf(paste0(outplot_dir,'peaks-GC-content.pdf',sep=''),width = 7, height = 7)
plot_promoters(promoters_gc,promoters_gc$gc_content,'Proportion of GC nucleotides')
dev.off()

# ## phastCons score
# library(phastCons7way.UCSC.hg38)

# ## expand peaks (i.e. get all bp within peaks from start to end) 
# promoters <- promoters[,peakID:=paste(type,1:nrow(promoters),sep='.')]

# get_phastCons_score <-function(df){
#     expanded_peaks <- copy(df)[,c(..range_keys,'peakID')]
#     expanded_peaks <- expanded_peaks[,list(start = seq(start, end)), by = peakID] 
#     expanded_peaks <- expanded_peaks[
#         df[,c('seqnames','peakID')],on='peakID',allow.cartesian=T
#         ][
#             ,end:=start
#     ]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

#     phastCons <- phastCons7way.UCSC.hg38

#     scores <- gscores(phastCons, expanded_peaks)%>%as.data.table()
#     rm(expanded_peaks)
#     scores <- scores[
#             ,avg_phastcons:=mean(default),by=.(peakID)
#             ][
#                 ,c('peakID','avg_phastcons')
#     ]%>%unique()%>%na.omit()
#     final_df <- copy(df)[scores,on='peakID',nomatch=0]
#     return(final_df)
# }

# promoters_scores <- promoters%>%split(by='type')%>%lapply(function(x)get_phastCons_score(x))

# promoters_scores <- rbindlist(promoters_scores)[,peakID:=NULL]

# pdf(paste0(outplot_dir,'peaks-phastcons-score.pdf',sep=''),width = 7, height = 7)
# plot_promoters(promoters_scores,promoters_scores$avg_phastcons,'mean phastCons score')
# dev.off()
