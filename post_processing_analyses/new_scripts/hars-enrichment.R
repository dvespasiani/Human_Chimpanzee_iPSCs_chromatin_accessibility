## script to test enrichment of da/all peaks on hars 
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)

options(width=150)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

source('./scripts/utils.R')

outplot_dir <- create_dir(plot_dir,'hars')
chrom_state_dir <- '../data/iPSC_chrom_states_hg38'

## read hars
hars <- fread('./data/hars-hg19-uebbing2019.txt',sep='\t',header=T)%>%setnames(old=c(2:4),new=c(range_keys))
setkeyv(hars,range_keys)

## get DA peaks
da_results <- read_da_results('da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC','peaktype')]
setkeyv(da_results,range_keys)

# ## read random genomic regions
# random_peaks <- fread('./output/files/random_regions/random_genomic_regions.txt',sep='\t',header=T)
# random_peaks <- random_peaks[,peakID:=paste('peak_',1:nrow(random_peaks),sep='')]
# setkeyv(random_peaks,range_keys)

## get tissues 
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
allpromoters <- ipsc_chromstate[chrom_state %like%  'Tss|Flnk'][,c('chrom_state','cell_type'):= NULL]%>%unique()
setkeyv(allpromoters,range_keys)

## overlap peaks hars
hars_overlap <- function(peaks){
  annotated <- foverlaps(copy(peaks),hars,type='any')%>%na.omit()
  annotated <- annotated[
    ,overlap:=ifelse(i.start<start,i.end-start,end-i.start)
    ][
        ,.SD[which.max(overlap)], by=.(Enhancer)
        ][
            ,c(range_keys[-1],'overlap'):=NULL
            ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  return(annotated)
}

mypeaks_harsoverlap <- hars_overlap(da_results)
# randompeaks_harsoverlap <- hars_overlap(random_peaks)
allpromoters_harsoverlap <- hars_overlap(allpromoters)


## OR
hars_or_matrix <- matrix(
    data=c(
        nrow(mypeaks_harsoverlap),
        nrow(allpromoters_harsoverlap),
        nrow(da_results)-nrow(mypeaks_harsoverlap),
        nrow(allpromoters)- nrow(allpromoters_harsoverlap)
        ),
    byrow=T,
    ncol=2
)
colnames(hars_or_matrix) = c('mypeaks','allpromoters')
rownames(hars_or_matrix) = c('overlap','nonoverlap')
fisher.test(hars_or_matrix)
# data:  hars_or_matrix
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.823146 4.879491
# sample estimates:
# odds ratio 
#    4.31617 

## check if DA peaks are enriched within hars compare to non-da peaks
da_hars_or_matrix <- matrix(
    data=c(
        nrow(mypeaks_harsoverlap[DA=='da']),
        nrow(mypeaks_harsoverlap[DA=='non_da']),
        nrow(da_results[DA=='da'])-nrow(mypeaks_harsoverlap[DA=='da']),
        nrow(da_results[DA=='non_da'])- nrow(mypeaks_harsoverlap[DA=='non_da'])
        ),
    byrow=T,
    ncol=2
)
colnames(da_hars_or_matrix) = c('da','nonda')
rownames(da_hars_or_matrix) = c('overlap','nonoverlap')
fisher.test(da_hars_or_matrix)
# data:  da_hars_or_matrix
# p-value = 0.1118
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7115901 1.0322931
# sample estimates:
# odds ratio 
#  0.8596088 

## then look at distribution number of hsubstitution per har based on uebbing et al data
pdf(paste(outplot_dir,'numb-hsubstitution-hars.pdf',sep=''),width=7,height = 7)
df <- rbind(
    copy(mypeaks_harsoverlap)[,c('No_hSubs','Enhancer','DA')],
    copy(hars)[,c('No_hSubs','Enhancer')][,DA:='all_hars']
    )
comparisons = list(
  c('all_hars','da'),
  c('all_hars','non_da'),
  c('da','non_da')
  )
ggplot(df,aes(x=DA,y=log10(No_hSubs+1),fill=DA))+
    geom_violin(trim=T,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
    scale_fill_manual(values = c('grey',da_palette),labels=c('all_hars',names(da_palette)))+
    xlab('Genomic element') + ylab('log10 number of human-specific substitutions')+
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
dev.off()
