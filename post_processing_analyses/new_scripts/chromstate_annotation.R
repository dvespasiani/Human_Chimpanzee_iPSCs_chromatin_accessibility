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
# outfile_dir <- create_dir(da_dir,genome)
outplot_dir <- create_dir(plot_dir,'chromstate_annotation')

genome <- 'hg38'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC'))
setkeyv(da_results,range_keys)

## get tissues 
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

## overlap ranges
peak_annotation <- foverlaps(da_results,ipsc_chromstate,type='any')[
  ,c(range_keys[-1]):=NULL
]%>%na.omit()%>%
  setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

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

da_peaks_counts <- copy(peak_counts_chromstate)[DA =='da'] 
nonda_peaks_counts <- copy(peak_counts_chromstate)[DA =='non_da'] 

calculate_or <- function(peaks_oi,peaks_noi,merging_keys){
  peaksoi_vs_peaksnoi_or <- merge(peaks_oi,peaks_noi,by=c(merging_keys))%>%
  dplyr::select(c(contains('numb'),contains(all_of(merging_keys))))

  fisher_test <- copy(peaksoi_vs_peaksnoi_or)%>%split(by=c(merging_keys))%>%
  lapply(
    function(x){
    x <- x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%fisher.test()
    x <- data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
        )
    }
  )
  or_results <- Map(mutate,fisher_test,elements=names(fisher_test))%>%rbindlist()%>%adjust_pvalues()
  return(or_results)
}

da_vs_nonda_or <- calculate_or(da_peaks_counts,nonda_peaks_counts,'chrom_state')

## plot the results
## get chromHMM colors
chromHMM_colors <- nihroadmap_colors(da_vs_nonda_or,'elements')

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

pdf(paste(outplot_dir,'OR_da_vs_nonda_peaks_chromstate.pdf',sep=''),width=8,height = 5)
plot_or(or = da_vs_nonda_or,ylab = 'odds ratio \n DA vs non DA peaks')
dev.off()

## hypothesis: enrichment in het because chimp ipscs are > primed than humans 
## check whether this enrichment is caused by peaks that are not accessible in humans 
species_accessibility <- copy(peak_annotation)[
  DA=='da'
  ][
    ,species:=ifelse(logFC<0,'chimp','human')
    ][
      ,c('cell_type','logFC'):=NULL
]%>%unique()

species_accessibility <- species_accessibility[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state,species)
      ][
        ,numb_peaks:=.N,by=.(species)
        ][
            ,c('species','chrom_state','numb_peaks_chromstate','numb_peaks')
]%>%unique()

chimp_peaks_counts <- copy(species_accessibility)[species =='chimp'] 
human_peaks_counts <- copy(species_accessibility)[species =='human'] 

human_vs_chimp_or <- calculate_or(human_peaks_counts,chimp_peaks_counts,'chrom_state')

pdf(paste(outplot_dir,'OR_human_vs_chimp_da_peaks_chromstate.pdf',sep=''),width=8,height = 5)
plot_or(or = human_vs_chimp_or,ylab = 'odds ratio \n humans vs chimp DA peaks')
dev.off()

