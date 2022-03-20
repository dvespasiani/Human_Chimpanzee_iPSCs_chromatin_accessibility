# ## script used to calculated the phastCons score for the set of common and species-specific peaks
# library(dplyr)
# library(data.table)
# library(magrittr)
# library(GenomicRanges)
# library(phastCons7way.UCSC.hg38)
# library(GenomicScores)
# library(ggthemes)
# library(ggplot2)
# library(ggpubr)

# options(width=150)
# setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

# scripts_dir <- './scripts/'
# source(paste(scripts_dir,'utils.R',sep=''))

# outplot_dir <- create_dir(plot_dir,'phastCons')

# ## get DA peaks
# da_file <- paste(da_dir,genome,'/da_results.txt',sep='')
# da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))%>%setorderv('peak_species',1)

# all_peaks <- split(da_results,by='DA')%>%lapply(function(x)x=x[,c(..range_keys,'peakID')]%>%unique())

# ## expand peaks (i.e. get all bp within peaks from start to end) 
# expanded_peaks <- copy(all_peaks)%>%
#     lapply(function(x)x=x[
#               ,list(start = seq(start, end)), by = peakID
#               ] 
# )

# expanded_peaks <- purrr::map2(expanded_peaks,all_peaks,full_join,by=c('peakID'))
# expanded_peaks <- lapply(expanded_peaks,function(x)x=x[
#     ,c('seqnames','start.x','peakID')
#     ][
#         ,end:=start.x
#         ]%>%setnames(old='start.x',new='start')%>%makeGRangesFromDataFrame(keep.extra.columns=T)
# )

# ## get phastcons scores
# phastCons <- phastCons7way.UCSC.hg38

# peak_scores <- lapply(expanded_peaks,function(x)gscores(phastCons, x)%>%as.data.table())
# peak_scores <- Map(mutate,peak_scores,DA = names(peak_scores))%>%rbindlist()

# peak_scores <- peak_scores[
#         ,avg_phastcons:=round(mean(default),2),by=.(peakID,DA)
#         ][
#             ,c('peakID','DA','avg_phastcons')
# ]%>%unique()%>%na.omit()

# all_peaks <- rbindlist(all_peaks)

# peak_scores <- peak_scores[
#   all_peaks,on='peakID',nomatch=0
#   ][
#     ,width:=end-start][
#       ,avg_phastcons:=avg_phastcons/width
# ]

# phastcons_permutation <- permute_data(peak_scores$avg_phastcons,peak_scores$DA)

# ## permutation plot
# pdf(paste0(outplot_dir,'permuted_phastcons_differences.pdf',sep=''),width = 7, height = 7)
# ggplot(phastcons_permutation$permuted_scores,aes(x=permuted_scores))+
# geom_histogram(binwidth=0.0005)+
# geom_vline(xintercept=phastcons_permutation$observed_value,linetype='dashed',color='red')+
# annotate(
#   "text", x=phastcons_permutation$observed_value-0.002, 
#   y=600, 
#   label = paste('zscore = ',round(phastcons_permutation$stat_results[[1]],2),' ', '-log10 pval = ',round(-log10(phastcons_permutation$stat_results[[2]]),2),sep='')) + 
#   xlab('mean difference phastCons7way score non-da vs da peaks')+ylab('counts')+theme_bw()
# dev.off()

# # test$avg_phastcons[is.na(test$avg_phastcons)]<-0
# ## divide peak ranges into deciles and for each of these compute the avg of phastcons scores 
# # peak_scores <- lapply(peak_scores,function(x)x=x[
# #     ,quantile_range := ntile(start, 10),by=.(peakID)
# #     ][
# #         ,avg_score_quantile:=round(mean(default),2),by=.(peakID,quantile_range)
# #         ][
# #             ,c('peakID','quantile_range','avg_score_quantile')
# #             ]%>%unique()
# # )

# # peak_scores <- Map(mutate,peak_scores,DA = names(peak_scores))%>%rbindlist()

# # ##------------------------------
# # ## test for signif differences
# # ##------------------------------
# # test_stat <- copy(peak_scores)%>%split(by='quantile_range')%>%
# # lapply(
# #     function(x)x <- compare_means(
# #     avg_score_quantile~DA,
# #     x,
# #     method='wilcox.test',
# #     ref.group='non_da'
# #     )
# # )
# # test_stat <- Map(mutate,test_stat,quantile_range=names(test_stat))%>%rbindlist()
# # test_stat <- test_stat[,c('p.adj','p.signif'):=NULL]%>%adjust_pvalues()
# # test_stat <- test_stat[,c('group2','quantile_range','p.adj','p.signif')][,quantile_range:=as.numeric(quantile_range)]%>%unique()%>%setnames(old=1,new='DA')
# # mock_stat <- copy(test_stat)[,DA:='non_da'][,p.adj:=0][,p.signif:=' ']
# # final_test_stat <- rbind(mock_stat,test_stat)

# ##----------------
# ## plot results
# ##----------------
# # pdf(paste(outplot_dir,'all_peaks_hg38_phastCons7way_scores.pdf',sep=''),width=10,height=7)
# # ggplot(peak_scores, aes(x = as.factor(quantile_range), y=avg_score_quantile,fill=DA))+
# #     geom_boxplot(position='dodge2',notch=TRUE)+
# #     scale_fill_manual(values=da_palette)+
# #     geom_text(data = final_test_stat, aes(x = quantile_range, y = 1.01, label = p.signif),size=7)+
# #     xlab('peak quantile')+ylab('avg phastCons7way score per quantile range')+
# #     theme(
# #         panel.background =element_rect(fill = 'white', colour = 'black',size=1),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major = element_blank(),
# #         legend.position = "bottom",
# #         legend.key = element_rect(fill = "white", colour = "black"),
# #         axis.line = element_blank()
# #         )
# # dev.off()

# # ## 
# # test <- copy(peak_scores)[
# #     ,quantile_range:=NULL
# #     ][
# #         ,avg_phastcons:=round(mean(avg_score_quantile),2),by=.(peakID,DA)
# #         ][
# #             ,avg_score_quantile:=NULL
# # ]%>%unique()      

# pdf(paste(outplot_dir,'single_score_peak_hg38_phastCons7way.pdf',sep=''),width=7,height=7)
# ggplot(peak_scores,aes(x=DA,y=avg_phastcons,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# scale_fill_manual(values=da_palette)+
# xlab('')+ylab('avg phastCons7way score per peak')+
# geom_hline(yintercept=median(peak_scores[DA=='da']$avg_phastcons),linetype='dashed')+
# theme_bw()+
# theme(
#   legend.position='none',
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()