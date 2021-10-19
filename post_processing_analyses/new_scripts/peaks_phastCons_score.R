## script used to calculated the phastCons score for the set of common and species-specific peaks
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(liftOver)
library(rtracklayer)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'phastCons')
# outfile_dir <- create_dir(outdir,'files/GO_enrich')
genome <- 'hg38'

## get DA peaks
da_file <- paste(da_dir,genome,'/da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC'))
setkeyv(da_results,range_keys)

all_peaks <- list(da_results[DA=='da'],da_results[DA!='da'])%>%
    lapply(function(x)x=x[,c(..range_keys,'peakID')])

## expand peaks (i.e. get all bp within peaks from start to end) 
expanded_peaks <- copy(all_peaks)%>%
    lapply(function(x)x=x[
              ,list(start = seq(start, end)), by = peakID
              ] 
)

expanded_peaks <- purrr::map2(expanded_peaks,all_peaks,full_join,by=c('peakID'))
expanded_peaks <- lapply(expanded_peaks,function(x)x=x[
    ,c('seqnames','start.x','peakID')
    ][
        ,end:=start.x
        ]%>%setnames(old='start.x',new='start')%>%makeGRangesFromDataFrame(keep.extra.columns=T)
)

## get phastcons scores
phastCons <- phastCons7way.UCSC.hg38

peak_scores <- lapply(expanded_peaks,function(x)gscores(phastCons, x)%>%as.data.table())

## divide peak ranges into deciles and for each of these compute the avg of phastcons scores 
peak_scores <- lapply(peak_scores,function(x)x=x[
    ,quantile_range := ntile(start, 10),by=.(peakID)
    ][
        ,avg_score_quantile:=round(mean(default),2),by=.(peakID,quantile_range)
        ][
            ,c('peakID','quantile_range','avg_score_quantile')
            ]%>%unique()
)

names(peak_scores) <- c('da','non_da')
peak_scores <- Map(mutate,peak_scores,peak_type = names(peak_scores))%>%rbindlist()

##------------------------------
## test for signif differences
##------------------------------
test_stat <- copy(peak_scores)%>%split(by='quantile_range')%>%
lapply(
    function(x)x <- compare_means(
    avg_score_quantile~peak_type,
    x,
    method='wilcox.test', 
    ref.group='non_da'
    )
)
test_stat <- Map(mutate,test_stat,quantile_range=names(test_stat))%>%rbindlist()
test_stat <- test_stat[,c('p.adj','p.signif'):=NULL]%>%adjust_pvalues()
test_stat <- test_stat[,c('group1','quantile_range','p.adj','p.signif')][,quantile_range:=as.numeric(quantile_range)]%>%unique()%>%setnames(old=1,new='peak_type')
nonda_mock_stat <- copy(test_stat)[,peak_type:='non_da'][,p.adj:=0][,p.signif:=' ']
final_test_stat <- rbind(nonda_mock_stat,test_stat)

##----------------
## plot results
##----------------
pdf(paste(outplot_dir,'all_peaks_hg38_phastCons7way_scores.pdf',sep=''),width=10,height=7)
ggplot(peak_scores, aes(x = as.factor(quantile_range), y=avg_score_quantile,fill=peak_type))+
    geom_boxplot(position='dodge2',notch=TRUE)+
    scale_fill_manual(values=da_palette)+
    geom_text(data = final_test_stat, aes(x = quantile_range, y = 1.01, label = p.signif),size=7)+
    xlab('peak quantile')+ylab('avg phastCons7way score per quantile range')+
    theme(
        panel.background =element_rect(fill = 'white', colour = 'black',size=1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", colour = "black"),
        axis.line = element_blank()
        )
dev.off()