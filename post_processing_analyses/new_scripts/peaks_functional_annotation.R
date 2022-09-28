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

options(width=150)
set.seed(2022)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

chrom_state_dir <- '../data/iPSC_chrom_states_hg38'
outplot_dir <- create_dir(plot_dir,'chromstate_annotation')
out_tabledir <- create_dir(table_dir,'funcAnnot')

## get DA peaks
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC','peaktype')]

## blacklisted regions
blacklist <- fread('../data/ENCODE_blacklisted/hg38_blacklist_v2.bed',sep='\t',header=F,col.names=range_keys)%>%makeGRangesFromDataFrame()

## get tissues 
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

# ## generate random genomic region
# numb_random_sets = 2
# random_peaks <- list()

# for (i in 1:numb_random_sets){
#   random_peaks[[i]] <- randomizeRegions(makeGRangesFromDataFrame(da_results),mask=blacklist,genome='hg38',allow.overlaps=F)%>%as.data.table()
# }
# random_peaks <- rbindlist(random_peaks)
# random_peaks <- random_peaks[,peakID:=paste('peak_',1:nrow(random_peaks),sep='')][seqnames%in% standard_chr]

## read random genomic regions

random_peaks <- fread('./output/random_regions/random_genomic_regions.txt',sep='\t',header=T)
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

## counts peaks per state per cell type 
## if start peak < start state then overlap = end peak - start state; 
# ## otherwise if start peak > start state then overlap = end state - start peak
# peak_annotation <- peak_annotation[
#     ,overlap:=ifelse(i.start<start,i.end-start,end-i.start),by=.(cell_type)
#     ][
#         ,.SD[which.max(overlap)], by=.(peakID,cell_type)
#         ][
#             ,c(range_keys[-1],'overlap'):=NULL
# ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
# # peak_annotation <- foverlaps(da_results,ipsc_chromstate,type='any')[,c('start','end'):=NULL]%>%na.omit()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

# peakid_chromstate <- copy(peak_annotation)[,c('peakID','DA','chrom_state')]%>%unique()

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

pdf(paste(outplot_dir,'mypeaks_vs_randomregion_enrichment_chromstates.pdf',sep=''),width=7,height = 7)
plot_or(genwide_peaks_or)
dev.off()

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

pdf(paste(outplot_dir,'da_vs_nonda_enrichment_chromstates.pdf',sep=''),width=7,height = 7)
plot_or(danonda_peaks_or)
dev.off()

## export or results
df2 <- copy(danonda_peaks_or)%>%setnames(old='elements',new='chrom_state')
fwrite(df2,paste(out_tabledir,'da_nonda_peaks_or_chromstate.txt',sep=''),sep='\t',quote=F,row.names=F,col.names=T)


## get independent sets of promoters and calculate
## 1) CpG and GC content
## 2) PhastCons score

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

promoter_states <- 'Tss|Flnk'

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

all_promoters <- genes(txdb,single.strand.genes.only=FALSE)%>%as.data.table()
all_promoters <- all_promoters[!seqnames%like% "M|Un|_"][,end:=start][,start:=start-500][,width:=end-start][,c(..range_keys)][,type:='upstream_txdb']

ipsc_promoters <- copy(ipsc_chromstate)[chrom_state %like% promoter_states][,c(1:3)][,type:='ipsc_promoters']%>%unique()
random_promoters <- copy(randompeaks_annotation)[chrom_state %like% promoter_states][,c(..range_keys)][,type:='random']%>%unique()
mypromoters <- copy(mypeaks_annotation)[chrom_state %like% promoter_states][,c(..range_keys,'DA')][,type:=DA][,DA:=NULL]%>%unique()

promoters <- rbind(mypromoters,ipsc_promoters,random_promoters,all_promoters)%>%unique()

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
  c('da','ipsc_promoters'),
  c('non_da','ipsc_promoters'),
  c('da','upstream_txdb'),
  c('non_da','upstream_txdb'),
  c('da','random'),
  c('non_da','random'),
  c('non_da','da')
)
## plot function
plot_promoters <- function(x,column_to_plot,ylab){

  df <- copy(x)[,column:=column_to_plot]
  
  type_order <- c('ipsc_promoters', 'upstream_txdb','random', 'da','non_da')
  
  p <- ggplot(df,aes(x=factor(type,level=type_order),y=column,fill=type))+
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
  scale_fill_manual(values=promoters_palette)+
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

pdf(paste0(outplot_dir,'peaks_CpG_content.pdf',sep=''),width = 7, height = 7)
plot_promoters(promoters_cpg,promoters_cpg$norm_cpg_content,'Normd CpG content')
dev.off()

pdf(paste0(outplot_dir,'peaks_gc_content.pdf',sep=''),width = 7, height = 7)
plot_promoters(promoters_gc,promoters_gc$gc_content,'Normd GC content')
dev.off()

## phastCons score
library(phastCons7way.UCSC.hg38)

## expand peaks (i.e. get all bp within peaks from start to end) 
promoters <- promoters[,peakID:=paste(type,1:nrow(promoters),sep='.')]

get_phastCons_score <-function(df){
    expanded_peaks <- copy(df)[,c(..range_keys,'peakID')]
    expanded_peaks <- expanded_peaks[,list(start = seq(start, end)), by = peakID] 
    expanded_peaks <- expanded_peaks[
        df[,c('seqnames','peakID')],on='peakID',allow.cartesian=T
        ][
            ,end:=start
    ]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

    phastCons <- phastCons7way.UCSC.hg38

    scores <- gscores(phastCons, expanded_peaks)%>%as.data.table()
    rm(expanded_peaks)
    scores <- scores[
            ,avg_phastcons:=mean(default),by=.(peakID)
            ][
                ,c('peakID','avg_phastcons')
    ]%>%unique()%>%na.omit()
    final_df <- copy(df)[scores,on='peakID',nomatch=0]
    return(final_df)
}

promoters_scores <- promoters%>%split(by='type')%>%lapply(function(x)get_phastCons_score(x))

promoters_scores <- rbindlist(promoters_scores)[,peakID:=NULL]

pdf(paste0(outplot_dir,'peaks_phastcons_score.pdf',sep=''),width = 7, height = 7)
plot_promoters(promoters_scores,promoters_scores$avg_phastcons,'avg phastCons score')
dev.off()

##--------------------------------------------------------
## Empirical permutation to assess genome-wide enrichment 
## of our set of peaks within any particular chrom state 
##------------------------------------------------------
# numb_permutations =10000

# ## number of peaks annotated in each chrom state
# genwide_observed_value = as.data.table(table(copy(peakid_chromstate)$chrom_state))%>%setnames(old=c(1:2),new=c('chrom_state','score'))%>%split(by='chrom_state')

# ## create fake list of genomic regions of same size distribution of our set of peaks
# ## use regioneR
# ## get blacklist regions for masking
# blacklist <- fread('../data/ENCODE_blacklisted/hg38_blacklist_v2.bed',sep='\t',header=F,col.names=range_keys)%>%makeGRangesFromDataFrame()

# da_results_gr <- makeGRangesFromDataFrame(da_results,keep.extra.columns=T)

# genwide_list_permutations = list()
# for(i in 1:(numb_permutations/100)){
#       random_gr <- randomizeRegions(da_results_gr,mask=blacklist,genome='hg38',allow.overlaps=F)%>%as.data.table()
#       random_gr <- random_gr[seqnames %in% standard_chr]
#       random_gr_annotation <- foverlaps(random_gr,ipsc_chromstate,type='any')
#       random_gr_annotation <- random_gr_annotation[,c(range_keys[-1],'cell_type'):=NULL]%>%na.omit()%>%unique()
#       permuted_scores=table(random_gr_annotation$chrom_state)
#       genwide_list_permutations[[i]] = as.data.table(permuted_scores)%>%setnames(old=c(1:2),new=c('chrom_state','permuted_score'))
# }
# genwide_permutation = rbindlist(genwide_list_permutations)%>%split(by='chrom_state')

# genwide_enrichment = purrr::map2(genwide_permutation,genwide_observed_value,function(p,o){
#   zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
#   pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#   stat =  data.table(zscore = zscore,pval=pvalue)
#   return(stat)
#   }
# )
# genwide_enrichment <- Map(mutate,genwide_enrichment,chrom_state=names(genwide_enrichment))%>%rbindlist()%>%setorderv('zscore',1)
# genwide_enrichment <- genwide_enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]


# ipsc_largest_state <-copy(ipsc_chromstate)[,cell_type:=NULL]%>%unique()
# ipsc_largest_state <- ipsc_largest_state[
#   ,.SD[which.max(end)], by=.(seqnames,start)
#   ][
#     ,.SD[which.min(start)], by=.(seqnames,end)
# ]

# observed_value = as.data.table(table(copy(peak_annotation)$chrom_state))%>%setnames(old=c(1:2),new=c('chrom_state','score'))%>%split(by='chrom_state')

# list_df_permuted_scores = list()
# for(i in 1:10000){
#       permuted_scores=table(sample(ipsc_largest_state$chrom_state, nrow(peak_annotation), replace = F))
#       list_df_permuted_scores[[i]] = as.data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('chrom_state','permuted_score'))
# }
# list_df_permuted_scores = rbindlist(list_df_permuted_scores)%>%split(by='chrom_state')

# enrichment = purrr::map2(list_df_permuted_scores,observed_value,function(p,o){
#   zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
#   pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#   stat =  data.table(zscore = zscore,pval=pvalue)
#   return(stat)
#   }
# )
# enrichment <- Map(mutate,enrichment,chrom_state=names(enrichment))%>%rbindlist()%>%setorderv('zscore',1)
# enrichment <- enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]

## get chromHMM colors (from utils.R)
# plot_enrichments <- function(enrichment){
#  p <- ggplot(enrichment,aes(x=reorder(chrom_state,-log_zscore),y=log_zscore,fill=chrom_state))+
#     geom_bar(stat='identity',position =  position_dodge())+
#     scale_fill_manual(values = chrom_state_colors)+
#     geom_hline(yintercept=0,linetype='dashed')+
#     xlab('chromatin state') + ylab('log zscore')+
#     theme_classic()+
#     theme(
#       legend.position = "bottom",
#       axis.text.x =element_blank(),
#       axis.ticks.x =element_blank()
#     )
#   return(p)
# }
# pdf(paste(outplot_dir,'genwide_enrichment_chromstates.pdf',sep=''),width=8,height = 5)
# plot_enrichments(genwide_enrichment)
# dev.off()


# ## permutation again to check enrichment vs DA vs nonDA 
# da_observed_value = as.data.table(table(copy(peakid_chromstate[DA=='da'])$chrom_state))%>%setnames(old=c(1:2),new=c('chrom_state','score'))%>%split(by='chrom_state')

# da_list_permutations = list()
# for(i in 1:numb_permutations){
#       new_da_col <- sample(peakid_chromstate$DA,nrow(peakid_chromstate),replace=F)
#       new_df <- copy(peakid_chromstate)[,DA:=new_da_col]
#       permuted_scores=table(new_df[DA=='da']$chrom_state)
#       da_list_permutations[[i]] = data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('chrom_state','permuted_score'))
# }
# da_permutation <- rbindlist(da_list_permutations)%>%split(by='chrom_state')

# da_enrichment = purrr::map2(da_permutation,da_observed_value,function(p,o){
#   zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
#   pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#   stat =  data.table(zscore = zscore,pval=pvalue)
#   return(stat)
#   }
# )
# da_enrichment <- Map(mutate,da_enrichment,chrom_state=names(da_enrichment))%>%rbindlist()%>%setorderv('zscore',1)
# da_enrichment <- da_enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]


# pdf(paste(outplot_dir,'da_enrichment_chromstates.pdf',sep=''),width=8,height = 5)
# plot_enrichments(da_enrichment)
# dev.off()

# peak_counts_chromstate <- copy(peak_annotation)[
#     ,cell_type:=NULL
# ]%>%unique()

# peak_counts_chromstate <- peak_counts_chromstate[
#       ,numb_peaks_chromstate:=.N,by=.(chrom_state,DA)
#       ][
#         ,numb_peaks:=.N,by=.(DA)
#         ][
#             ,c('DA','chrom_state','numb_peaks_chromstate','numb_peaks')
# ]%>%unique()

# nonda_peaks <- copy(peak_counts_chromstate)[DA =='non_da'] 
# da_peaks <- copy(peak_counts_chromstate)[DA =='da'] 

# peaks_or <- calculate_or(da_peaks,nonda_peaks,'chrom_state')%>%adjust_pvalues()
# peaks_or <-peaks_or[,6:=NULL]

# ## plot the results
# plot_or <- function(or,ylab){
#   p <- ggplot(or, aes(x=factor(elements,levels=chrom_states), y=odds_ratio,label = p.signif)) + 
#   geom_point(aes(colour = elements))+
#   geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
#   scale_colour_manual(values = chrom_state_colors)+
#   geom_hline(yintercept=1,linetype='dashed',size=.5)+
#   geom_text(aes(y = or$upper_ci+0.1),size=5)+
#   xlab(' ')+ylab(ylab)+
#   theme_classic()+
#   theme(
#     legend.position = "none",
#     axis.text.x =element_blank(),
#     axis.ticks.x =element_blank()
#   )
#   return(p)
# }

# pdf(paste(outplot_dir,'OR_da_vs_nonda_peaks_chromstate.pdf',sep=''),width=8,height = 5)
# plot_or(or = peaks_or,ylab = 'OR DA vs non-DA peaks')
# dev.off()

# ## permutation
# calculate_or_permutations <- function(peaks){
#   df <- copy(peaks)
#   df <- df[
#     ,numb_peaks_chromstate:=.N,by=.(chrom_state,DA)
#     ][
#         ,numb_peaks:=.N,by=.(DA)
#         ][
#             ,c('DA','chrom_state','numb_peaks_chromstate','numb_peaks')
#             ]%>%unique()
#   non_da_peaks <- copy(df)[DA =='non_da'] 
#   da_peaks <- copy(df)[DA =='da'] 
#   odds_ratio_df <- calculate_or(non_da_peaks,da_peaks,'chrom_state')
#   odds_ratio_vector <- odds_ratio_df$odds_ratio
#   names(odds_ratio_vector) = odds_ratio_df$elements
#   return(odds_ratio_vector)
# }

# score_distribution = list()

# for(i in 1:10000){
#     new_da_col = sample(peak_annotation$DA, length(peak_annotation$DA), replace = F)
#     new_peak_annot_df <- copy(peak_annotation)[,DA:=new_da_col]
#     permuted_odds_ratio <- calculate_or_permutations(new_peak_annot_df)
#     dat <- data.table(chrom_state = names(permuted_odds_ratio), odds_ratio = permuted_odds_ratio)
#     score_distribution[[i]] <- dat # add it to your list
    
# }
# score_distribution =  rbindlist(score_distribution)%>%split(by='chrom_state')

# observed_values = copy(peaks_or)[,c('odds_ratio','elements')]%>%split(by='elements')

# permutation_results <- purrr::map2(observed_values,score_distribution,function(x,y){
#    zscore <- (x$odds_ratio-mean(y$odds_ratio))/sd(y$odds_ratio)
#    pvalue <- 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#    return <- copy(y)[,zscore:=zscore][,pval:=pvalue][,observed_value:=x$odds_ratio]
#    return(return)
# })

# permutation_results <- Map(mutate,permutation_results,chrom_state = names(permutation_results))%>%rbindlist()

# names(chrom_state_colors) = chrom_states

# pdf(paste(outplot_dir,'permuted_odds_ratios.pdf',sep=''),width=10,height = 7)
# ggplot(permutation_results,aes(x=odds_ratio,col=chrom_state))+
# scale_colour_manual(values=chrom_state_colors)+
# geom_histogram(binwidth=0.005)+
# xlab('permutation odds ratio')+ 
# geom_vline(data = permutation_results, aes(xintercept=observed_value),linetype='dashed',color='red')+
# facet_wrap(chrom_state~.,scale='free_y',ncol=3)+
# theme_bw()+theme(legend.position='none')
# dev.off()


# ## get independent sets of promoters and calculate
# ## 1) CpG and GC content
# ## 2) PhastCons score

# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(BSgenome.Hsapiens.UCSC.hg38)

# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# all_promoters <- transcripts(txdb)%>%as.data.table()
# all_promoters <- all_promoters[!seqnames%like% "M|Un|_"][,end:=start][,start:=start-500][,width:=end-start][,c(1:3)][,type:='upstream_txdb']

# ipsc_promoters <- copy(ipsc_chromstate)[chrom_state %like% 'Tss|Flnk'][,c(1:3)][,type:='ipsc_promoters']%>%unique()

# mypromoters <- copy(peak_annotation)[chrom_state %like% 'Tss|Flnk'][,c(..range_keys,'DA')][,type:=DA][,DA:=NULL]

# promoters <- rbind(mypromoters,ipsc_promoters,all_promoters)

# ## CpG & GC
# cpg <- copy(promoters)%>%makeGRangesFromDataFrame()%>%Repitools::cpgDensityCalc(organism=Hsapiens)
# promoters_cpg <- copy(promoters)[,cpg_content:=cpg][,norm_cpg_content:=cpg/(end-start)]

# promoters_gc <- copy(promoters)[
#     ,sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
#     ][
#       !sequence%like%'N' ## do this because u need to remove sequences with Ns 
#       ][
#         ,gc_content:=letterFrequency(DNAStringSet(sequence),"GC")/width(DNAStringSet(sequence))
#         ][
#           ,sequence:=NULL
# ]

# comparisons = list(
#   c('da','ipsc_promoters'),
#   c('non_da','ipsc_promoters'),
#   c('da','upstream_txdb'),
#   c('non_da','upstream_txdb'),
#   c('non_da','da')
# )
# ## plot function
# plot_promoters <- function(x,column_to_plot,ylab){
#   promoters_palette = c('#94D2BD','#F4A261',da_palette)
#   names(promoters_palette) = c('ipsc_promoters','upstream_txdb',names(da_palette))
#   df <- copy(x)[,column:=column_to_plot]
  
#   type_order <- c('ipsc_promoters', 'upstream_txdb', 'da','non_da')
  
#   p <- ggplot(df,aes(x=factor(type,level=type_order),y=column,fill=type))+
#   geom_violin(trim=T,scale = "width")+
#   geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
#   scale_fill_manual(values=promoters_palette)+
#   xlab(' ')+ylab(ylab)+
#   stat_compare_means(
#   method = "wilcox.test",
#   comparisons = comparisons,
#   size=5
#   )+
#   theme_classic()+
#   theme(
#     legend.position = "bottom",
#     axis.text.x =element_blank(),
#     axis.ticks.x =element_blank()
#   )
#   return(p)
# }

# pdf(paste0(outplot_dir,'peaks_CpG_content.pdf',sep=''),width = 7, height = 7)
# plot_promoters(promoters_cpg,promoters_cpg$norm_cpg_content,'Normd CpG content')
# dev.off()

# pdf(paste0(outplot_dir,'peaks_gc_content.pdf',sep=''),width = 7, height = 7)
# plot_promoters(promoters_gc,promoters_gc$gc_content,'Normd GC content')
# dev.off()

# ## phastCons score
# library(phastCons7way.UCSC.hg38)

# ## expand peaks (i.e. get all bp within peaks from start to end) 
# promoters <- promoters[,peakID:=paste(type,1:nrow(promoters),sep='.')]

# get_phastCons_score <-function(df){
#     expanded_peaks <- copy(df)[,c(..range_keys,'peakID')]%>%unique()
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

# promoters_scores = get_phastCons_score(promoters)

# # expanded_peaks <- copy(promoters)[,list(start = seq(start, end)), by = peakID] 
# # expanded_peaks <- expanded_peaks[promoters[,c('seqnames','peakID')],on='peakID',allow.cartesian=T]
# # expanded_peaks <- expanded_peaks[,end:=start]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

# # phastCons <- phastCons7way.UCSC.hg38

# # promoter_scores <- gscores(phastCons, expanded_peaks)%>%as.data.table()

# # promoter_scores <- promoter_scores[
# #         ,avg_phastcons:=mean(default),by=.(peakID)
# #         ][
# #             ,c('peakID','avg_phastcons')
# # ]%>%unique()%>%na.omit()

# promoters_scores <- promoters_scores[,type:=gsub('\\..*','',peakID)][,peakID:=NULL]

# pdf(paste0(outplot_dir,'peaks_phastcons_score.pdf',sep=''),width = 7, height = 7)
# plot_promoters(promoters_scores,promoters_scores$avg_phastcons,'avg phastCons score')
# dev.off()