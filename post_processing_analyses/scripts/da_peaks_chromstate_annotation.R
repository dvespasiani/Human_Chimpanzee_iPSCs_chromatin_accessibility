## script to test the enrichment across each chromstate cell type combination
## of common DA peaks (with different regulation) vs species specific DA peaks
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(liftOver)
library(rtracklayer)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

chrom_state_dir = './data/iPSC_chrom_states_hg38'
plot_dir = './post_processing_analyses/output/plots/chromstate_annotation/'
peakDir = './post_processing_analyses/output/DA/peaks/'

## using peaks after DA analysis (both da and non DA)
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))

## convert pantro5 to hg38
chimp_peaks_hg38 = convert_coord(chimp_peaks,'panTro5ToHg38.over.chain')

species_specific_peaks = rbind(chimp_peaks_hg38,human_peaks)

all_peaks = list(common_regions,species_specific_peaks)
names(all_peaks) = c('common','specific')

## get tissues 
ipsc_chromstate = read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

## overlap ranges
peak_annotation = lapply(all_peaks,function(x)foverlaps(x,ipsc_chromstate,type='any')%>%na.omit()) ## removes peaks in sex chromosomes
peak_annotation = lapply(peak_annotation,function(x)x=x[,c(range_keys[-1]):=NULL]%>%
    setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
)
##----------------------------------------------------
## Calculate odds ratio of 
## common vs species specific peaks across each state
##----------------------------------------------------
## to make the odds ratio calculations simpler dont split by cell type but remove that info
## and simply count the annotations
count_peaks_chromstate_celltype = function(x){
  df=copy(x)
  df =df[
    ,cell_type:=NULL
    ]%>%unique()
    df=df[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state)
      ][
        ,numb_peaks:=.N
        ][
          ,prop_peaks_chromstate:=numb_peaks_chromstate/numb_peaks
          ][
            ,c('numb_peaks_chromstate','numb_peaks','prop_peaks_chromstate','chrom_state')
            ]%>%unique()
    return(df)
}

all_peak_annotation = peak_annotation%>%lapply(function(y)count_peaks_chromstate_celltype(y))

## Calculate hypergeometric enrichment and odds ratio of common vs species-specific peaks across chromstate/cell types
make_enrichment_df = function(common,specific){
  df = merge(common,specific,by=c('chrom_state'))

  df = df[
          ,enrichment := prop_peaks_chromstate.x/prop_peaks_chromstate.y ## this is the odds ratio common vs specific
          ][
              order(readr::parse_number(gsub("^.*\\.", "",df$chrom_state))),
              ]

    df_enrichment = copy(df)[,c('enrichment','chrom_state')]%>%unique()

    df_numb_peaks = copy(df)[
      ,c('numb_peaks_chromstate.x','numb_peaks_chromstate.y','numb_peaks.x','numb_peaks.y','chrom_state')
      ][
        ,numb_common_peaks_non_chromstate:=numb_peaks.x-numb_peaks_chromstate.x
        ][
          ,numb_species_peaks_non_chromstate:=numb_peaks.y-numb_peaks_chromstate.y
          ]%>%
          unique()%>%
          setnames(
            old=c('numb_peaks_chromstate.x','numb_peaks_chromstate.y'),
            new=c('numb_common_peaks_chromstate','numb_species_chromstate')
            )%>%
            dplyr::select(
              -c(contains('.x'),contains('.y'))
            )

    # colnames(df_numb_peaks)[1:4]=c('numb_common_peaks_celltype_chromstate','numb_species_peaks_celltype_chromstate','numb_common_peaks_celltype','numb_species_peaks_celltype')

  return(list(df_enrichment,df_numb_peaks))
}

peak_enrichment_df =  make_enrichment_df(all_peak_annotation$common,all_peak_annotation$specific)

fisher_test = copy(peak_enrichment_df[[2]])%>%
  split(by=c('chrom_state'))%>%
    lapply(function(x)x=x[
      ,c('chrom_state'):=NULL]%>%
      summarize_if(is.numeric, sum, na.rm=TRUE)%>%
      as.numeric()%>% ## keep this otherwise fisher.test doesnt accept matrix
      matrix(nrow=2,byrow=T)%>%
      fisher.test()
)

fisher_test_results = copy(fisher_test)%>%
  lapply(function(x)
    x=data.table(
      'p'=x$p.value,
      'odds_ratio'=x$estimate,
      'lower_ci'=x$conf.int[[1]],
      'upper_ci'=x$conf.int[[2]]
    )
)
fisher_test_results = Map(mutate,fisher_test_results,chrom_state=names(fisher_test_results))%>%rbindlist()
fisher_test_results = fisher_test_results[
  ,significance:=ifelse(p<0.05,'*','')
]

## plot the results
## get chromHMM colors
chromHMM_colors=nihroadmap_colors(fisher_test_results,'chrom_state')

state_levels=fisher_test_results$chrom_state

pdf(paste(plot_dir,'odds_ratio_common_vs_all_specific_peaks.pdf',sep=''),width=8,height = 5)
ggplot(fisher_test_results, aes(x=factor(chrom_state,levels=state_levels), y=odds_ratio,label = significance)) + 
geom_point(aes(colour = chrom_state))+
geom_errorbar(aes(ymin=lower_ci, ymax=upper_ci,colour = chrom_state),width=0,position=position_dodge(0.05))+
scale_colour_manual(values = chromHMM_colors)+
geom_hline(yintercept=1,linetype='dashed',size=.5)+
geom_text(aes(y= fisher_test_results$upper_ci+0.1),size=5)+
xlab(' ')+ylab('odds ratio \n common vs species specific peaks')+
theme(
        panel.background =element_rect(fill = 'white', colour = 'black',size=1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white", colour = "black"),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()
