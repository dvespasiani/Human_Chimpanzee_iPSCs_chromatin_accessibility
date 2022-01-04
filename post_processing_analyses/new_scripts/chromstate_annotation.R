## script to test the enrichment across each chromstate cell type combination
## of common DA peaks (with different regulation) vs species specific DA peaks
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

chrom_state_dir <- '../data/iPSC_chrom_states_hg38'
outplot_dir <- create_dir(plot_dir,'chromstate_annotation')

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))
setkeyv(da_results,range_keys)

## get tissues 
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

## overlap ranges
peak_annotation <- foverlaps(da_results,ipsc_chromstate,type='any')[
  ,c(range_keys[-1]):=NULL
]%>%na.omit()%>%
  setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  
peak_annotation <- peak_annotation[
    ,pleiotropy:=.N,by=.(peakID,chrom_state)
    ][
      ,.SD[which.max(abs(pleiotropy))], by=.(peakID)
]

##-------------------------------------------------------
## Calculate odds ratio of 
## DA(species-specific) vs non-DA (common) windows(peaks)
##-------------------------------------------------------
peak_counts_chromstate <- copy(peak_annotation)[
    ,cell_type:=NULL
]%>%unique()

peak_counts_chromstate <- peak_counts_chromstate[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state,DA)
      ][
        ,numb_peaks:=.N,by=.(DA)
        ][
            ,c('DA','chrom_state','numb_peaks_chromstate','numb_peaks')
]%>%unique()

common_peaks <- copy(peak_counts_chromstate)[DA =='non_da'] 
species_specific_peaks <- copy(peak_counts_chromstate)[DA =='da'] 

peaks_or <- calculate_or(common_peaks,species_specific_peaks,'chrom_state')


## plot the results
## get chromHMM colors
chromHMM_colors <- nihroadmap_colors(peaks_or,'elements')

plot_or <- function(or,ylab){
  p <- ggplot(or, aes(x=factor(elements,levels=chrom_states), y=odds_ratio,label = p.signif)) + 
  geom_point(aes(colour = elements))+
  geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = elements),width=0,position=position_dodge(0.05))+
  scale_colour_manual(values = chromHMM_colors)+
  geom_hline(yintercept=1,linetype='dashed',size=.5)+
  geom_text(aes(y = or$upper_ci+0.1),size=5)+
  xlab(' ')+ylab(ylab)+
  theme(
    panel.background =element_rect(fill = 'white', colour = 'black',size=1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "none",
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
  )
  return(p)
}

pdf(paste(outplot_dir,'OR_common_vs_species_specific_peaks_chromstate.pdf',sep=''),width=8,height = 5)
plot_or(or = peaks_or,ylab = 'odds ratio \n non da vs da peaks')
dev.off()

## permutation
calculate_or_permutations <- function(peaks){
  df <- copy(peaks)
  df <- df[
    ,numb_peaks_chromstate:=.N,by=.(chrom_state,DA)
    ][
        ,numb_peaks:=.N,by=.(DA)
        ][
            ,c('DA','chrom_state','numb_peaks_chromstate','numb_peaks')
            ]%>%unique()
  non_da_peaks <- copy(df)[DA =='non_da'] 
  da_peaks <- copy(df)[DA =='da'] 
  odds_ratio_df <- calculate_or(non_da_peaks,da_peaks,'chrom_state')
  odds_ratio_vector <- odds_ratio_df$odds_ratio
  names(odds_ratio_vector) = odds_ratio_df$elements
  return(odds_ratio_vector)
}

score_distribution = list()

for(i in 1:10000){
    new_da_col = sample(peak_annotation$DA, length(peak_annotation$DA), replace = F)
    new_peak_annot_df <- copy(peak_annotation)[,DA:=new_da_col]
    permuted_odds_ratio <- calculate_or_permutations(new_peak_annot_df)
    dat <- data.table(chrom_state = names(permuted_odds_ratio), odds_ratio = permuted_odds_ratio)
    score_distribution[[i]] <- dat # add it to your list
    
}
score_distribution =  rbindlist(score_distribution)%>%split(by='chrom_state')

observed_values = copy(peaks_or)[,c('odds_ratio','elements')]%>%split(by='elements')

permutation_results <- purrr::map2(observed_values,score_distribution,function(x,y){
   zscore <- (x$odds_ratio-mean(y$odds_ratio))/sd(y$odds_ratio)
   pvalue <- 2*pnorm(q=abs(zscore), lower.tail=FALSE)
   return <- copy(y)[,zscore:=zscore][,pval:=pvalue][,observed_value:=x$odds_ratio]
   return(return)
})

permutation_results <- Map(mutate,permutation_results,chrom_state = names(permutation_results))%>%rbindlist()

names(chrom_state_colors) = chrom_states

pdf(paste(outplot_dir,'permuted_odds_ratios.pdf',sep=''),width=10,height = 7)
ggplot(permutation_results,aes(x=odds_ratio,col=chrom_state))+
scale_colour_manual(values=chrom_state_colors)+
geom_histogram(binwidth=0.005)+
xlab('permutation odds ratio')+ 
geom_vline(data = permutation_results, aes(xintercept=observed_value),linetype='dashed',color='red')+
facet_wrap(chrom_state~.,scale='free_y',ncol=3)+
theme_bw()+theme(legend.position='none')
dev.off()
