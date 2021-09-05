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
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

plot_dir = './post_processing_analyses/output/plots/phastCons/'
peakDir = './post_processing_analyses/output/DA/peaks/'

## using peaks after DA analysis (both da and non DA)
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))

## convert chimp to hg38
chimp_peaks_hg38 = convert_coord(chimp_peaks,'panTro5ToHg38.over.chain')

all_peaks = list(common_regions,chimp_peaks_hg38,human_peaks)%>%
    lapply(function(x)x=x[
          ,c(..range_keys)
          ][
              ,peakID:=paste('peak_',1:nrow(x),sep='') ## need to change the peakID here
          ]
)

## expand peaks (i.e. get all bp within peaks from start to end) 
expanded_peaks = copy(all_peaks)%>%
    lapply(function(x)x=x[
              ,list(start = seq(start, end)), by = peakID
              ] 
)

expanded_peaks = purrr::map2(expanded_peaks,all_peaks,full_join,by=c('peakID'))
expanded_peaks = lapply(expanded_peaks,function(x)x=x[
    ,c('seqnames','start.x','peakID')
    ][
        ,end:=start.x
        ]%>%setnames(old='start.x',new='start')%>%makeGRangesFromDataFrame(keep.extra.columns=T)
)

## get phastcons scores
phastCons = phastCons7way.UCSC.hg38

peak_scores = lapply(expanded_peaks,function(x)gscores(phastCons, x)%>%as.data.table())

## divide peak ranges into deciles and for each of these compute the avg of phastcons scores 
peak_scores = lapply(peak_scores,function(x)x=x[
    ,quantile_range := ntile(start, 10),by=.(peakID)
    ][
        ,avg_score_quantile:=round(mean(default),2),by=.(peakID,quantile_range)
        ][
            ,c('peakID','quantile_range','avg_score_quantile')
            ]%>%unique()
)

names(peak_scores)= c('common','chimp','human')
peak_scores = Map(mutate,peak_scores,'peak'=names(peak_scores))%>%rbindlist()

##------------------------------
## test for signif differences
##------------------------------
test_stat=copy(peak_scores)%>%split(by='quantile_range')%>%
lapply(
    function(x)x=compare_means(
    avg_score_quantile~peak,
    x,
    method='wilcox.test', 
    ref.group='common'
    )
)
test_stat = Map(mutate,test_stat,quantile_range=names(test_stat))%>%rbindlist()
test_stat = test_stat[,c('p.adj','p.signif'):=NULL]%>%adjust_pvalues()

test_stat=test_stat[,c('group1','quantile_range','p.adj','p.signif')][,quantile_range:=as.numeric(quantile_range)]%>%unique()%>%setnames(old=1,new='peak')

human_mock_stat = copy(test_stat)[,peak:='human'][,p.adj:=0][,p.signif:=' ']
chimp_mock_stat = copy(test_stat)[,peak:='chimp'][,p.adj:=0][,p.signif:=' ']

final_test_stat = rbind(human_mock_stat,chimp_mock_stat,test_stat)

##----------------
## plot results
##----------------

pdf(paste(plot_dir,'all_peaks_hg38_phastCons7way_scores.pdf',sep=''),width=10,height=7)
ggplot(peak_scores, aes(x = as.factor(quantile_range), y=avg_score_quantile,fill=peak))+
    geom_boxplot(position='dodge2',notch=TRUE)+
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