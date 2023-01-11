library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)
library(purrr)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

# ctcf_dir <- './output/homer/homer_output/ctcf/'
outplot_dir <- create_dir(plot_dir,'ctcf')
peakDir = './output/DA/peaks/'
tads_dir =  './output/TADs/'
files_dir <- './output/files/'
homer_outdir <-'output/homer/homer_output/'

##-----------------------
## read CTCF motifs
##-----------------------
ctcf_motif_files <- list.files(homer_outdir,full.names=T,recursive=T,pattern=c('discovered_ctcf.txt'))
names(ctcf_motif_files) = c('all_peaks','random')

ctcf_motifs <- lapply(ctcf_motif_files,function(x){
  results <- fread(x,sep='\t',header=T,drop=c("Motif Name"))%>%setnames(old=1,new='peakID')
  results <- results[,motif_length:=nchar(Sequence)][,rank_score:=round(MotifScore/max(MotifScore),2)]
  return(results)
  }
)
ctcf_motifs <- Map(mutate,ctcf_motifs,tf='ctcf')


## add genomic region
da_results <- read_da_results('da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]
 
rand_peaks <- fread('./output/files/random_regions/random_genomic_regions.txt',sep='\t',header=T)
rand_peaks <- rand_peaks[,peakID:=paste('peak_',1:nrow(rand_peaks),sep='')]

peaks <- list(da_results,rand_peaks)

ctcf_peaks <- purrr::map2(ctcf_motifs,peaks,function(c,p){
    c <- c[p,on='peakID',nomatch=0
    ][
        ,peak_center:=start+(end-start)
        ][
            ,motif_start:=peak_center+Offset
            ][
                ,motif_end:=motif_start+motif_length
            ][
                ,.SD[which.max(MotifScore)], by=.(peakID)
                ]
    return(c)
})


nonctcf_peaks <- purrr::map2(peaks,ctcf_peaks,function(p,c){
    nonctcf <- copy(p)[!peakID %in% c$peakID][,tf := 'other']
    return(nonctcf)
})
names(nonctcf_peaks) = c('all_peaks','random')
##------------------------------------------------------------------
## compare the overlap of da/non-da ctcf at common TAD boundaries
## relative to all genome-wide ctcf overlap
##------------------------------------------------------------------
orth_tads <- read_tads('human_chimp_orth_tads.txt')
chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL][species=='Chimp']
human_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

chimp_tads_hg38 <- copy(orth_tads)[species=='Chimp']
colnames(chimp_tads_hg38)[4:6] = range_keys
chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'species'),nomatch=0][
  ,c(range_keys):=NULL
]%>%setnames(old=c(1:3),new=range_keys)

common_tads_hg38 <- copy(human_tads_hg38)[species=='Common']

ssp_tads <- rbind(chimp_tads_hg38,human_tads_hg38[species!='Common'])[,species:='ssp']
all_tads <- list(common_tads_hg38,ssp_tads)

all_tads <- lapply(all_tads,function(x){
  x <- setorderv(x,range_keys,1)
  x<-x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]
  setkeyv(x,range_keys)
  return(x)
  }
) 
names(all_tads) = c('common','ssp')

## define TAD boundaries
boundary =  50000
boundaries = copy(all_tads)%>%
lapply(
    function(x)x=x[
    ,tadID:=paste('tad_',seq(1:nrow(x)),sep='')
    ][
        ,upstream_start:= start-boundary
        ][
            ,upstream_end:= start
            ][
                ,downstream_start:= end 
                ][
                    ,downstream_end:= end+ boundary
            ]
)

# boundaries = copy(all_tads)[
#         ,upstream_start:= start-boundary
#         ][
#             ,upstream_end:= start
#             ][
#                 ,downstream_start:= end 
#                 ][
#                     ,downstream_end:= end+ boundary
# ]

upstream_boundaries = copy(boundaries)%>%
lapply(
    function(x)x=x%>%dplyr::select(c('seqnames',contains('upstream'),'species','tadID'))%>%
    setnames(old=c('upstream_start','upstream_end'),new=range_keys[-1])%>%
    mutate('boundary'='upstream')
)
downstream_boundaries = copy(boundaries)%>%
lapply(
    function(x)x=x%>%dplyr::select(c('seqnames',contains('downstream'),'species','tadID'))%>%
    setnames(old=c('downstream_start','downstream_end'),new=range_keys[c(2,3)])%>%
    mutate('boundary'='downstream')
)

tad_boundaries = map2(upstream_boundaries,downstream_boundaries,rbind)%>%rbindlist()
# lapply(tad_boundaries,function(x)setkeyv(x,range_keys))
setkeyv(tad_boundaries,range_keys)

## look at overlap @ TAD boundaries
tadbound_overlap <-function(peaks,tadbound){
    overlap <- foverlaps(peaks,tadbound,type='any')%>%na.omit()
    overlap <- overlap[
        ,distance_peak_tad:= ifelse(boundary=='upstream',i.end-end,i.start-start)
        ][
            ,.SD[which.min(abs(distance_peak_tad))], by=.(peakID,species)
            ][
                ,c('start','end','boundary'):=NULL
                ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
    return(overlap)
}

# ctcf_tadbound_overlap <- lapply(tad_boundaries,function(x)tadbound_overlap(ctcf_motifs,x))
# nonctct_tadbound_overlap <- lapply(tad_boundaries,function(x)tadbound_overlap(nonctcf_peaks,x))

ctcf_tadbound_overlap <- lapply(ctcf_peaks,function(c)tadbound_overlap(c,tad_boundaries))
nonctct_tadbound_overlap <- lapply(nonctcf_peaks,function(n)tadbound_overlap(n,tad_boundaries)) 
rand_tadbound_overlap <- tadbound_overlap(rand_peaks,tad_boundaries)

## count number of ctcf (da & non da) associated to tad boundaries
associatedCTCF <- map2(copy(ctcf_tadbound_overlap),copy(nonctct_tadbound_overlap),function(a,n){
    numbAssociated <- length(unique(a$peakID))
    numbNonAssociated <- length(unique(n$peakID))
    total = numbAssociated + numbNonAssociated
    prop = numbAssociated/total
    return <- data.table(numbAssociated = numbAssociated,numbNonAssociated = numbNonAssociated, proportion = prop)
    return(return)
})
# $all_peaks
#    numbAssociated numbNonAssociated proportion
# 1:          12926             31428  0.2914281

# $random
#    numbAssociated numbNonAssociated proportion
# 1:           2089             35372 0.05576466

## CALCULATE ENRICHMNET OF CTCF AT TAD BOUNDARIES VIA FISHER'S OR
fisher_matrix <- data.frame(
    tad_associated = c(associatedCTCF$all_peaks$numbAssociated,associatedCTCF$random$numbAssociated),
    tad_nonassociated = c(associatedCTCF$all_peaks$numbNonAssociated,associatedCTCF$random$numbNonAssociated)
)
rownames(fisher_matrix) = c('myctcf','random')
fisher_matrix = as.matrix(fisher_matrix)
fisher.test(fisher_matrix)
# data:  fisher_matrix
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  6.632578 7.314685
# sample estimates:
# odds ratio 
#   6.964355 

## test whether da ctcf are enriched @ species specific TADs 
da_ctcf_ssp <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA=='da' & species == 'ssp']$peakID))
nonda_ctcf_ssp <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA!='da' & species == 'ssp']$peakID))
da_ctcf_common <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA=='da' & species != 'ssp' ]$peakID))
nonda_ctcf_common <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA!='da' & species != 'ssp' ]$peakID))

## enrichment via OR
fisher_matrix_dactcf <- data.frame(
    sp_associated=c(da_ctcf_ssp,nonda_ctcf_ssp),
    common_associated=c(da_ctcf_common,nonda_ctcf_common)
)%>%as.matrix()

fisher.test(fisher_matrix_dactcf)
# data:  fisher_matrix_dactcf
# p-value = 0.912
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9253411 1.0705861
# sample estimates:
# odds ratio 
#  0.9953094 

## Now consider all CTCFs that overlap tad boundaries 
## and w/o dividing them into common/sp-specific
## plot cumulative distribution ctcf peaks and other peaks @ tad boundaries
calculate_cumdist <- function(peaks,tf){
    cumdist <- copy(peaks)[
        ,.SD[which.min(abs(distance_peak_tad))], by=.(peakID)
        ][
            ,binned_dist:=plyr::round_any(abs(distance_peak_tad), 10)
            ][
                abs(distance_peak_tad)<=100000
                ][
                    ,tf:=tf
                ][
                    ,c('tf','binned_dist','peakID')
                    ]
    return(cumdist)
}
cumdist_ctcfpeaks_tadbound = calculate_cumdist(ctcf_tadbound_overlap[[1]],'ctcf')
cumdist_nonctcfpeaks_tadbound = calculate_cumdist(nonctct_tadbound_overlap[[1]],'other')
cumditst_rand_tadbound_overlap <- calculate_cumdist(rbindlist(rand_tadbound_overlap),'random')

cumdist <- rbind(cumdist_ctcfpeaks_tadbound,cumdist_nonctcfpeaks_tadbound,cumditst_rand_tadbound_overlap)

pdf(paste(outplot_dir,'cumulative_distribution_peaks_tad_distances.pdf',sep=''),width = 7,height = 7)
palette <- c(promoters_palette[3],plurip_tf_palette[c(4,5)])
ggplot(cumdist, aes(binned_dist,col=tf)) + 
stat_ecdf(pad=F)+
scale_color_manual(values=palette)+
xlab('Distance from TAD (bp)')+ylab('Cumulative distribution')+
theme_classic()
dev.off()

## calculate correlation number of peaks distance from TAD boundary
cor_ctcf <- copy(cumdist_ctcfpeaks_tadbound)[
    ,binned_dist:=plyr::round_any(abs(binned_dist), 10000)
    ][
        ,numbpeaks:=.N,by=.(binned_dist)
]
cor.test(cor_ctcf$numbpeaks,cor_ctcf$binned_dist)
# data:  cor_ctcf$numbpeaks and cor_ctcf$binned_dist
# t = -481.82, df = 12924, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9741696 -0.9723509
# sample estimates:
#        cor 
# -0.9732755 


cor_nonctcf <- copy(cumdist_nonctcfpeaks_tadbound)[
    ,binned_dist:=plyr::round_any(abs(binned_dist), 10000)
    ][
        ,numbpeaks:=.N,by=.(binned_dist)
]
cor.test(cor_nonctcf$numbpeaks,cor_nonctcf$binned_dist)

# 	Pearson's product-moment correlation

# data:  cor_nonctcf$numbpeaks and cor_nonctcf$binned_dist
# t = -242.23, df = 42692, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7647803 -0.7567899
# sample estimates:
#       cor 
# -0.7433908 

cor_rand <- copy(cumditst_rand_tadbound_overlap)[
    ,binned_dist:=plyr::round_any(abs(binned_dist), 10000)
    ][
        ,numbpeaks:=.N,by=.(binned_dist)
]
cor.test(cor_rand$numbpeaks,cor_rand$binned_dist)
# 	Pearson's product-moment correlation

# data:  cor_rand$numbpeaks and cor_rand$binned_dist
# t = -267.93, df = 81601, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6877454 -0.6804451
# sample estimates:
#        cor 
# -0.7089402 

# simplify_motifs <- function(motifs){
#     simplified <- copy(motifs)[
#         ,c('peakID','seqnames','motif_start','motif_end','Strand')
#         ]%>%setnames(old=c(3:4),new=range_keys[-1])
#     simplified <- simplified[,.SD[which.max(end)], by=.(peakID)][,.SD[which.min(start)], by=.(peakID)]
#     return(simplified)
# }


simplifiedCTCFMotifs <- lapply(copy(ctcf_peaks),function(x){
    simplified <- copy(x)[
        ,c('peakID','seqnames','motif_start','motif_end','Strand')
        ]%>%setnames(old=c(3:4),new=range_keys[-1])
    simplified <- simplified[,.SD[which.max(end)], by=.(peakID)][,.SD[which.min(start)], by=.(peakID)]
    return(simplified)
})

## calculate conservation score for these ctcfs and compare across peaks
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

ctcfpeaks_conservation <- lapply(simplifiedCTCFMotifs,function(x)get_phastCons_score(x))

myctcfpeaks_conservation <- copy(ctcfpeaks_conservation[[1]])[,tadbound:=ifelse(peakID %in% cumdist_ctcfpeaks_tadbound$peakID,'tadbound','no_tadbound')][,tf:='ctcf']

## plot conservation
# conservation <- rbind(ctcfpeaks_conservation)

## plot function
plot_ctcf_tad = function(x,column,ylab,xlab){
 tadass_palette <- c('#577399','#BDD5EA')
    df=copy(x)[
        ,column_to_plot:=column
    ]
    p = ggplot(df, aes(x=tadbound, y=column_to_plot,fill=tadbound)) +
    geom_violin(trim=T,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
    ylab(ylab)+xlab(xlab)+
    scale_fill_manual(values=tadass_palette)+
    stat_compare_means(
        method = "wilcox.test",
        label.y = max(df$column_to_plot+0.1),
        size=5
        )+
    theme_classic()+
    theme(
        legend.position = 'none',
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
    return(p)
}

pdf(paste(outplot_dir,'conservation_score_tf_tadbound.pdf',sep=''),width = 7,height = 7)
plot_ctcf_tad(myctcfpeaks_conservation,myctcfpeaks_conservation$avg_phastcons,ylab='avg phastCons score',xlab=' ') 
dev.off()

## read and plot Deepbind affinity predictions
deepbind_results <- fread('./output/files/deepbind_output/ctcf_affinity_predictions.txt',sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(1,2,3),new=c('peakID','Sequence','prediction'))

ctcf_motifs_affinitites <- copy(ctcf_peaks[[1]])[
    deepbind_results,on=c('peakID'),nomatch=0
    ][
        ,tadbound:=ifelse(peakID %in% cumdist_ctcfpeaks_tadbound$peakID,'tadbound','no_tadbound')
]

pdf(paste(outplot_dir,'ctcf_peaks_deepbind.pdf',sep=''),width = 7,height = 7)
plot_ctcf_tad(ctcf_motifs_affinitites,ctcf_motifs_affinitites$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
dev.off()
