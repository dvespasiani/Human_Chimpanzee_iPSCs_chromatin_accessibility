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
ctcf_motif_files <- list.files(homer_outdir,full.names=T,recursive=T,pattern=c('ctcf_discovered.txt'))

ctcf_motifs<- lapply(ctcf_motif_files,function(x){
  results <- fread(x,sep='\t',header=T,drop=c("Motif Name"))%>%setnames(old=1,new='peakID')
  results <- results[,motif_length:=nchar(Sequence)][,rank_score:=round(MotifScore/max(MotifScore),2)]
  return(results)
  }
)
ctcf_motifs <- Map(mutate,ctcf_motifs,tf='ctcf')%>%rbindlist()

## add genomic region
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

ctcf_motifs <- ctcf_motifs[
  da_results,on='peakID',nomatch=0
  ][
    ,peak_center:=start+(end-start)
    ][
      ,motif_start:=peak_center+Offset
      ][
        ,motif_end:=motif_start+motif_length
]

nonctcf_peaks <- copy(da_results)[!peakID %in% ctcf_motifs$peakID]

## read random peaks 
rand_peaks <- fread('./output/random_regions/random_genomic_regions.txt',sep='\t',header=T)
rand_peaks <- rand_peaks[,peakID:=paste('peak_',1:nrow(rand_peaks),sep='')]

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

tad_boundaries = purrr::map2(upstream_boundaries,downstream_boundaries,rbind)
lapply(tad_boundaries,function(x)setkeyv(x,range_keys))

## look at overlap @ TAD boundaries
tadbound_overlap <-function(peaks,tadbound){
    overlap <- foverlaps(peaks,tadbound,type='any')%>%na.omit()
    overlap <- overlap[
        ,distance_peak_tad:= ifelse(boundary=='upstream',i.end-end,i.start-start)
        ][
            ,.SD[which.min(abs(distance_peak_tad))], by=.(peakID,species)
            ][
                ,c('start','end','species','boundary'):=NULL
                ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
    return(overlap)
}

ctcf_tadbound_overlap <- lapply(tad_boundaries,function(x)tadbound_overlap(ctcf_motifs,x))
nonctct_tadbound_overlap <- lapply(tad_boundaries,function(x)tadbound_overlap(nonctcf_peaks,x))
rand_tadbound_overlap <- lapply(tad_boundaries,function(x)tadbound_overlap(rand_peaks,x))

## count number of ctcf (da & non da) associated to tad boundaries
ctcf_overlap_peakIDs <- unique(rbindlist(copy(ctcf_tadbound_overlap))$peakID)
numb_ctcf_associated<- length(ctcf_overlap_peakIDs)
numb_ctcf_nonassociated<- length(unique(ctcf_motifs[!peakID %in% ctcf_overlap_peakIDs]$peakID))

round(numb_ctcf_associated/(numb_ctcf_nonassociated+numb_ctcf_associated)*100,2)
# [1] 49.22

## nonctcf
nonctcf_overlap_peakIDs <- unique(rbindlist(copy(nonctct_tadbound_overlap))$peakID)
numb_nonctcf_associated<- length(nonctcf_overlap_peakIDs)
numb_nonctcf_nonassociated<- length(unique(da_results[!peakID %in% nonctcf_overlap_peakIDs])$peakID)

round(numb_nonctcf_associated/(numb_nonctcf_nonassociated+numb_nonctcf_associated)*100,2)
# [1] 28.63

## random
random_overlap_peakIDs <- unique(rbindlist(copy(rand_tadbound_overlap))$peakID)
numb_random_associated<- length(random_overlap_peakIDs)
numb_random_nonassociated<- length(unique(rand_peaks[!peakID %in% random_overlap_peakIDs])$peakID)

round(numb_random_associated/(numb_random_nonassociated+numb_random_associated)*100,2)
# [1] 28.63

## enrichment via OR
fisher_matrix <- data.frame(
    tad_associated=c(numb_ctcf_associated,numb_nonctcf_associated),
    tad_nonassociated=c(numb_ctcf_nonassociated,numb_nonctcf_nonassociated)
)
rownames(fisher_matrix) = c('ctcf','other')
fisher_matrix = as.matrix(fisher_matrix)
fisher.test(fisher_matrix)
# 	Fisher's Exact Test for Count Data

# data:  fisher_matrix
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.313282 2.522777
# sample estimates:
# odds ratio 
#    2.41574 

## test whether da ctcf are enriched @ species specific TADs 
da_ctcf_spsp <- length(unique(copy(ctcf_tadbound_overlap[[2]])[DA=='da']$peakID))
nonda_ctcf_spsp <- length(unique(copy(ctcf_tadbound_overlap[[2]])[DA!='da']$peakID))
da_ctcf_common <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA=='da']$peakID))
nonda_ctcf_common <- length(unique(copy(ctcf_tadbound_overlap[[1]])[DA!='da']$peakID))

## enrichment via OR
fisher_matrix_dactcf <- data.frame(
    sp_associated=c(da_ctcf_spsp,nonda_ctcf_spsp),
    common_associated=c(da_ctcf_common,nonda_ctcf_common)
)%>%as.matrix()

fisher.test(fisher_matrix_dactcf)
# 	Fisher's Exact Test for Count Data

# data:  fisher_matrix_dactcf
# p-value = 0.1961
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.9657734 1.1860531
# sample estimates:
# odds ratio 
#   1.070214

## Now consider all CTCFs that overlap tad boundaries 
## and w/o dividing them into common/sp-specific
## plot cumulative distribution ctcf peaks and other peaks @ tad boundaries
calculate_cumdist <- function(peaks,tf){
    cumdist<-rbindlist(copy(peaks))[
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
cumdist_ctcfpeaks_tadbound = calculate_cumdist(ctcf_tadbound_overlap,'ctcf')
cumdist_nonctcfpeaks_tadbound = calculate_cumdist(nonctct_tadbound_overlap,'other')
cumditst_rand_tadbound_overlap <- calculate_cumdist(rand_tadbound_overlap,'random')

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

# 	Pearson's product-moment correlation

# data:  cor_ctcf$numbpeaks and cor_ctcf$binned_dist
# t = -211.06, df = 4970, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9512044 -0.9456203
# sample estimates:
#       cor 
# -0.948486 


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
# -0.760814 

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
# -0.6841124 

simplify_motifs <- function(motifs){
    simplified <- copy(motifs)[
        ,c('peakID','seqnames','motif_start','motif_end','Strand')
        ]%>%setnames(old=c(3:4),new=range_keys[-1])
    simplified <- simplified[,.SD[which.max(end)], by=.(peakID)][,.SD[which.min(start)], by=.(peakID)]
    return(simplified)
}

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

ctcfpeaks_conservation <- get_phastCons_score(simplify_motifs(ctcf_motifs))[,tadbound:=ifelse(peakID %in% cumdist_ctcfpeaks_tadbound$peakID,'tadbound','no_tadbound')][,tf:='ctcf']

## plot conservation
conservation <- rbind(ctcfpeaks_conservation)

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
plot_ctcf_tad(conservation,conservation$avg_phastcons,ylab='avg phastCons score',xlab=' ') 
dev.off()

## read and plot Deepbind affinity predictions
deepbind_results <- fread('./output/files/deepbind_output/ctcf_affinity_predictions.txt',sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(1,2,3),new=c('peakID','Sequence','prediction'))

ctcf_motifs_affinitites <- copy(ctcf_motifs)[
    deepbind_results,on=c('peakID','Sequence'),nomatch=0
    ][
        ,tadbound:=ifelse(peakID %in% cumdist_ctcfpeaks_tadbound$peakID,'tadbound','no_tadbound')
]

pdf(paste(outplot_dir,'ctcf_peaks_deepbind.pdf',sep=''),width = 7,height = 7)
plot_ctcf_tad(ctcf_motifs_affinitites,ctcf_motifs_affinitites$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
dev.off()



# ## keep only da and non-da ctcf respectively associated to ssp and common TAD boundaries
# common_ctcf_common_tads <- copy(ctcf_tadbound_overlap[[1]])[DA=='non_da']
# spsp_ctcf_spsp_tads <- copy(ctcf_tadbound_overlap[[2]])[DA=='da']

# refined_ctcf_list <- list(common_ctcf_common_tads,spsp_ctcf_spsp_tads)
# names(refined_ctcf_list) = names(ctcf_tadbound_overlap)

# density_myctcf_boundary_overlap <- copy(refined_ctcf_list)%>%lapply(function(x){
#     x<-x[
#         ,c(..range_keys,'distance_ctcf_tad','peakID','DA')
#             ][
#                 ,binned_abs_dist:=plyr::round_any(abs(distance_ctcf_tad), 200)
#                 ][
#                     ,numb_peaks_binned_dist:=.N,by=.(binned_abs_dist)
#                     ]
# })
# lapply(density_myctcf_boundary_overlap,function(x)cor.test(x$numb_peaks_binned_dist,x$binned_abs_dist,method='spearman'))
# # $common
# # 	Spearman's rank correlation rho
# # data:  x$numb_peaks_binned_dist and x$binned_abs_dist
# # S = 1.5155e+10, p-value < 2.2e-16
# # alternative hypothesis: true rho is not equal to 0
# # sample estimates:
# #        rho 
# # -0.8776644 

# # $ssp
# # 	Spearman's rank correlation rho
# # data:  x$numb_peaks_binned_dist and x$binned_abs_dist
# # S = 688545, p-value = 2.946e-09
# # alternative hypothesis: true rho is not equal to 0
# # sample estimates:
# #       rho 
# # -0.473833 

# ## plot correlation CTCF densities and distance from tad
# pdf(paste(outplot_dir,'ctcf_corr_peaks_tad_distance.pdf',sep=''),width= 7,height = 7)
# df <- copy(density_myctcf_boundary_overlap)%>%rbindlist()
# ggplot(df,aes(x=binned_abs_dist,y=numb_peaks_binned_dist,colour=DA))+
# geom_point(size=2)+
# geom_smooth(method=loess, se=FALSE, linetype="dashed",size=1)+
# scale_color_manual(values=da_palette)+
# # facet_wrap(DA~.,ncol=1)+
# xlab('Distance from TAD @ 200bp resolution')+ylab('Number of CTCF binding sites')+
# theme_classic()+
# theme(
#     legend.position = 'bottom',
#      strip.background = element_blank()
# )
# dev.off()

# ##-------------------------------------------------------------------------
# ## perform comparison between CTCF associated and non to TAD boundaries
# ##-------------------------------------------------------------------------
# ctcf_peaks_tad_info <- copy(ctcf_motifs)[
#     ,overlap:=ifelse(peakID %in% rbindlist(copy(ctcf_tadbound_overlap))$peakID,'tad_bound','not_associated')
# ]%>%setnames(old='Sequence',new='sequence')

# ## get position ctcf motifs
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(stringr)

# ctcf_motif_positions <- copy(ctcf_peaks_tad_info)[
#     ,c(..range_keys,'peakID','sequence')
#     ][
#         ,full_sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
# ]
# ctcf_motif_positions <- ctcf_motif_positions%>%split(by='sequence')%>%lapply(function(x){
#     ctcf_seq <-x$sequence
#     full_seq <- x$full_sequence
#     locations <- str_locate(full_seq,ctcf_seq)%>%as.data.table()
#     return(locations)
# })%>%rbindlist()
# colnames(ctcf_motif_positions) =c('offset_start','offset_end')

# myctcf_peaks_tad_info <- cbind(copy(ctcf_peaks_tad_info),ctcf_motif_positions)
# myctcf_peaks_tad_info <- myctcf_peaks_tad_info[
#     ,start:=start+offset_start
#     ][
#         ,end:=start+20
#         ][
#             ,c('offset_end','offset_start'):=NULL
# ]

# ## Get conservation score 
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

# ctcf_peaks_conScore <- copy(myctcf_peaks_tad_info)[,.SD[which.max(abs(end))], by=.(peakID)]%>%get_phastCons_score()

# ## plot function
# plot_ctcf_tad = function(x,column,ylab,xlab){
#  tadass_palette <- c('#577399','#BDD5EA')
#     df=copy(x)[
#         ,column_to_plot:=column
#     ]
#     p = ggplot(df, aes(x=overlap, y=column_to_plot,fill=overlap)) +
#     geom_violin(trim=T,scale = "width")+
#     geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
#     ylab(ylab)+xlab(xlab)+
#     scale_fill_manual(values=tadass_palette)+
#     facet_wrap(DA~.,ncol=3)+
#     stat_compare_means(
#         method = "wilcox.test",
#         label.y = max(df$column_to_plot+0.1),
#         size=5
#         )+
#     theme_classic()+
#     theme(
#         legend.position = 'none',
#         strip.background = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
#     return(p)
# }

# ## plot cons scores
# pdf(paste(outplot_dir,'ctcf_tad_phastcons.pdf',sep=''),width= 7,height = 7)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$avg_phastcons,ylab='PhastCons score',xlab=' ') 
# dev.off()

# ## plot Deepbind CTCF affinities
# pdf(paste(outplot_dir,'ctcf_peaks_deepbind.pdf',sep=''),width= 7,height = 7)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
# dev.off()






















# ## read deepbind input with all the info 
# deepbind_output_file <- dir(paste(files_dir,'deepbind_output',sep=''),recursive=F,full.names=T)
# ctcf_peak_files <- dir(homer_ctcf_dir,recursive=T,full.names=T,pattern='all')

# ctcf_peaks <- lapply(ctcf_peak_files,function(x)fread(x,header=T,sep='\t'))

# ## add DA info and genomic coordinates
# da_file <- paste(da_dir,'da_results.txt',sep='')
# da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'logFC','DA','peakID','peak_species'))

# myctcf_peaks <- copy(ctcf_peaks[[2]])[,'Motif Name':=NULL]%>%setnames(old=c('FASTA ID','Sequence'),new=c('peakID','sequence'))
# myctcf_peaks <- myctcf_peaks[da_results,on='peakID',nomatch=0]

# ctcf_peaks <- list(ctcf_peaks[[1]],myctcf_peaks)

# ctcf_affinities <- lapply(
#     deepbind_output_file,function(x)
#     fread(x,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(1,2,3),new=c('peakID','sequence','prediction'))
# )

# ctcf_peaks <- purrr::map2(ctcf_peaks,ctcf_affinities,function(x,y){
#     z<-copy(x)[copy(y),on=c('peakID','sequence'),nomatch=0]
#     setkeyv(z,range_keys)
#     return(z)
# })

# names(ctcf_peaks) = c('genome_wide','mypeaks')

# ## get number and % peaks with CTCF compared to genome wide
# totnumb_ctcf_genwide = copy(ctcf_peaks[[1]])[,c('peakID')]%>%unique()%>%nrow()
# totnumb_ctcf_mypeaks = copy(ctcf_peaks[[2]])[,c('peakID','DA')]%>%unique()%>%setorderv('DA')%>%split(by='DA')%>%lapply(function(x)nrow(x))
# prop_ctcf_mypeaks = lapply(totnumb_ctcf_mypeaks,function(x)round((x/totnumb_ctcf_genwide)*100,2))

# ## plot number of CTCF binding sites per peak (normalised by peak size)
# pdf(paste(outplot_dir,'homer_numb_ctcfs_per_peaktype.pdf',sep=''),width= 7,height = 7)
# df <- copy(ctcf_peaks[[2]])[,c('peakID','DA','sequence','prediction')]%>%unique()
# df <- df[,numbctcf:=.N,by=.(peakID,DA)]
# ggplot(df,aes(x=numbctcf,fill=DA))+
# geom_bar()+
# facet_wrap(DA~.,ncol=2)+
# scale_fill_manual(values=da_palette)+
# xlab('')+ylab('Number CTCF motifs per peak')+
# theme_classic()+
# theme(legend.position='bottom')
# dev.off()

# ## check if there is significant enrichment of DA CTCFs compared to non-DA CTCFs 
# totnumb_nonctcf_mypeaks = copy(da_results)[,numb:=.N,by='DA'][,c('DA','numb')]%>%unique()%>%setorderv('DA')%>%split(by='DA')
# totnumb_nonctcf_mypeaks = purrr::map2(totnumb_nonctcf_mypeaks,totnumb_ctcf_mypeaks,function(x,y){
#     return <- copy(x)$numb-y
#     return(return)
# })

# fisher_matrix <- matrix(
#     c(
#         totnumb_ctcf_mypeaks[[1]],totnumb_nonctcf_mypeaks[[1]],
#         totnumb_ctcf_mypeaks[[2]],totnumb_nonctcf_mypeaks[[2]]
#     ),nrow=2,ncol=2,byrow=TRUE
# )
# fisher.test(fisher_matrix)
# ## NB: tested DA vs non DA
# # Fisher's Exact Test for Count Data
# # data:  fisher_matrix
# # p-value = 1.901e-08
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  0.7005606 0.8462347
# # sample estimates:
# # odds ratio 
# #  0.7706156 


# ## plot proportion CTCF in my peaks relative to genome-wide
# pdf(paste(outplot_dir,'homer_prop_of_ctcfs_genwide.pdf',sep=''),width= 7,height = 7)
# df <- copy(prop_ctcf_mypeaks)%>%lapply(function(x)data.table(prop=x))
# df <- Map(mutate,df,DA=names(df))%>%rbindlist()
# ggplot(df,aes(x=DA,y=prop,fill=DA))+
# geom_bar(stat='identity',position=position_dodge())+
# # geom_density(alpha=0.5)+
# scale_fill_manual(values=da_palette)+
# xlab('')+ylab('Proportion of all genome-wide predicted CTCF binding sites')+
# theme_classic()+
# theme()
# dev.off()

# ##------------------------------------------------------------------
# ## compare the overlap of da/non-da ctcf at common TAD boundaries
# ## relative to all genome-wide ctcf overlap
# ##------------------------------------------------------------------
# orth_tads <- read_tads('human_chimp_orth_tads.txt')
# chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL][species=='Chimp']
# human_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

# chimp_tads_hg38 <- copy(orth_tads)[species=='Chimp']
# colnames(chimp_tads_hg38)[4:6] = range_keys
# chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'species'),nomatch=0][
#   ,c(range_keys):=NULL
# ]%>%setnames(old=c(1:3),new=range_keys)

# common_tads_hg38 <- copy(human_tads_hg38)[species=='Common']

# ssp_tads <- rbind(chimp_tads_hg38,human_tads_hg38[species!='Common'])[,species:='ssp']
# # all_tads <- list(chimp_tads_hg38,common_tads_hg38,human_tads_hg38[species!='Common'])
# all_tads <- list(common_tads_hg38,ssp_tads)

# all_tads <- lapply(all_tads,function(x){
#   x <- setorderv(x,range_keys,1)
#   x<-x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]
#   setkeyv(x,range_keys)
#   return(x)
#   }
# ) 
# names(all_tads) = c('common','ssp')

# ## define TAD boundaries
# boundary =  50000
# boundaries = copy(all_tads)%>%
# lapply(
#     function(x)x=x[
#     ,tadID:=paste('tad_',seq(1:nrow(x)),sep='')
#     ][
#         ,upstream_start:= start-boundary
#         ][
#             ,upstream_end:= start
#             ][
#                 ,downstream_start:= end 
#                 ][
#                     ,downstream_end:= end+ boundary
#             ]
# )

# upstream_boundaries = copy(boundaries)%>%
# lapply(
#     function(x)x=x%>%dplyr::select(c('seqnames',contains('upstream'),'species','tadID'))%>%
#     setnames(old=c('upstream_start','upstream_end'),new=range_keys[-1])%>%
#     mutate('boundary'='upstream')
# )
# downstream_boundaries = copy(boundaries)%>%
# lapply(
#     function(x)x=x%>%dplyr::select(c('seqnames',contains('downstream'),'species','tadID'))%>%
#     setnames(old=c('downstream_start','downstream_end'),new=range_keys[c(2,3)])%>%
#     mutate('boundary'='downstream')
# )

# tad_boundaries = purrr::map2(upstream_boundaries,downstream_boundaries,rbind)
# lapply(tad_boundaries,function(x)setkeyv(x,range_keys))

# ## look at overlap @ TAD boundaries
# ctcf_boundary_overlap <- lapply(tad_boundaries,function(x){
#     overlap <- foverlaps(ctcf_peaks[[2]],x,type='any')%>%na.omit()
#     overlap <- overlap[
#         ,distance_ctcf_tad:= ifelse(boundary=='upstream',i.end-end,i.start-start)
#         ][
#             ,.SD[which.min(abs(distance_ctcf_tad))], by=.(peakID,species)
#             ][
#                 ,c('start','end','species','boundary'):=NULL
#                 ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
#                 }
# )
# ## get some % CTCF peaks within tab boundaries 
# totnumb_ctcfs_associated <- copy(ctcf_boundary_overlap)%>%lapply(function(x){x=x[,numb_ctcfs:=.N,by=.(DA)][,c('numb_ctcfs','DA')]%>%unique()%>%setorderv('DA',1)%>%split(by='DA')%>%lapply(function(x)x$numb_ctcfs)})
# prop_ctcfs_associated <- lapply(totnumb_ctcfs_associated,function(x){
#     purrr::map2(totnumb_ctcf_mypeaks,x,function(tot,subset){
#     prop<-round((subset/tot)*100,2)
#     return(prop)
#     })
# })
# ## number and proportion of CTCF associated to either common or ssp TAD boundaries 
# ## -------------------------
# # $common$da
# # [1] 111
# # $common$non_da
# # [1] 3645
# # $ssp$da
# # [1] 141
# # $ssp$non_da
# # [1] 3971
# ## -------------------------
# # $common$da
# # [1] 21.89
# # $common$non_da
# # [1] 27.18
# # $ssp$da
# # [1] 27.81
# # $ssp$non_da
# # [1] 29.61
# ## -------------------------

# ## test enrichment for: 
# ## 1) non-DA CTCF-peaks @ common TAD boundaries
# matrix_nonda_common <- as.matrix(
#         data.table(
#             col1=c(
#                 totnumb_ctcfs_associated[[1]][[2]],
#                 totnumb_ctcfs_associated[[1]][[1]]
#                 ),
#             col2=c(
#                 totnumb_ctcf_mypeaks[[2]]-totnumb_ctcfs_associated[[1]][[2]],
#                 totnumb_ctcf_mypeaks[[1]]-totnumb_ctcfs_associated[[1]][[1]]
#                 )
#             )
# )
# fisher_nonda_common <- fisher.test(matrix_nonda_common)
# fisher_nonda_common
# # 	Fisher's Exact Test for Count Data
# # data:  matrix
# # p-value = 0.008016
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  1.072238 1.664520
# # sample estimates:
# # odds ratio 
# #   1.331476 
# ##-----------------------------------------------
# ## 2) DA CTCF-peaks @ ssp TAD boundaries
# matrix_da_ssp <- as.matrix(
#         data.table(
#             col1=c(
#                 totnumb_ctcfs_associated[[2]][[1]],
#                 totnumb_ctcfs_associated[[2]][[2]]
#                 ),
#             col2=c(
#                 totnumb_ctcf_mypeaks[[1]]-totnumb_ctcfs_associated[[2]][[1]],
#                 totnumb_ctcf_mypeaks[[2]]-totnumb_ctcfs_associated[[2]][[2]]
#                 )
#             )
# )
# fisher_da_ssp <- fisher.test(matrix_da_ssp)
# fisher_da_ssp
# # 	Fisher's Exact Test for Count Data
# # data:  matrix
# # p-value = 0.3995
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  0.7459866 1.1193649
# # sample estimates:
# # odds ratio 
# #  0.9158273 
# ##-----------------------------------------------

# ## keep only da or non-da ctcf respectively associated to ssp and common TAD boundaries
# non_da_ctcf_common <- copy(ctcf_boundary_overlap[[1]])[DA=='non_da']
# da_ctcf_ssp <- copy(ctcf_boundary_overlap[[2]])[DA=='da']

# refined_ctcf_list <- list(non_da_ctcf_common,da_ctcf_ssp)
# names(refined_ctcf_list) = names(ctcf_boundary_overlap)

# density_myctcf_boundary_overlap <- copy(refined_ctcf_list)%>%lapply(function(x){
#     x<-x[
#         ,c(..range_keys,'distance_ctcf_tad','peakID','prediction','DA')
#             ][
#                 ,binned_abs_dist:=plyr::round_any(abs(distance_ctcf_tad), 200)
#                 ][
#                     ,numb_peaks_binned_dist:=.N,by=.(binned_abs_dist)
#                     ]
# })
# lapply(density_myctcf_boundary_overlap,function(x)cor.test(x$numb_peaks_binned_dist,x$binned_abs_dist,method='spearman'))
# # $common
# # 	Spearman's rank correlation rho
# # data:  x$numb_peaks_binned_dist and x$binned_abs_dist
# # S = 1.5155e+10, p-value < 2.2e-16
# # alternative hypothesis: true rho is not equal to 0
# # sample estimates:
# #        rho 
# # -0.8776644 

# # $ssp
# # 	Spearman's rank correlation rho
# # data:  x$numb_peaks_binned_dist and x$binned_abs_dist
# # S = 688545, p-value = 2.946e-09
# # alternative hypothesis: true rho is not equal to 0
# # sample estimates:
# #       rho 
# # -0.473833 


# ## plot correlation CTCF densities and distance from tad
# pdf(paste(outplot_dir,'homer_ctcf_corr_peaks_tad_distance.pdf',sep=''),width= 7,height = 7)
# df <- copy(density_myctcf_boundary_overlap)%>%rbindlist()
# ggplot(df,aes(x=binned_abs_dist,y=numb_peaks_binned_dist,col=DA))+
# geom_point(size=2)+
# scale_color_manual(values=da_palette)+
# xlab('Distance from TAD @ 200bp resolution')+ylab('Number of CTCF binding sites')+
# facet_wrap(DA~.,ncol=2)+
# geom_smooth(method=loess, se=FALSE, linetype="dashed",color="black")+
# theme_classic()+
# theme(
#     legend.position = 'bottom',
#      strip.background = element_blank()
# )
# dev.off()

# ##-------------------------------------------------------------------------
# ## perform comparison between CTCF associated and non to TAD boundaries
# ##-------------------------------------------------------------------------
# ctcf_peaks_tad_info <- copy(ctcf_peaks[[2]])[,overlap:=ifelse(peakID %in% rbindlist(copy(ctcf_boundary_overlap))$peakID,'tad_bound','not_associated')]

# ## get position ctcf motifs
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(stringr)

# ctcf_motif_positions <- copy(ctcf_peaks_tad_info)[
#     ,c(..range_keys,'peakID','sequence')
#     ][
#         ,full_sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
# ]
# ctcf_motif_positions <- ctcf_motif_positions%>%split(by='sequence')%>%lapply(function(x){
#     ctcf_seq <-x$sequence
#     full_seq <- x$full_sequence
#     locations <- str_locate(full_seq,ctcf_seq)%>%as.data.table()
#     return(locations)
# })%>%rbindlist()
# colnames(ctcf_motif_positions) =c('offset_start','offset_end')

# myctcf_peaks_tad_info <- cbind(copy(ctcf_peaks_tad_info),ctcf_motif_positions)
# myctcf_peaks_tad_info <- myctcf_peaks_tad_info[
#     ,start:=start+offset_start
#     ][
#         ,end:=start+20
#         ][
#             ,c('offset_end','offset_start'):=NULL
# ]

# ## Get conservation score 
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

# ctcf_peaks_conScore <- copy(myctcf_peaks_tad_info)[,.SD[which.max(abs(end))], by=.(peakID)]%>%get_phastCons_score()

# # plot function
# plot_ctcf_tad = function(x,column,ylab,xlab){
# tadass_palette <- c('#577399','#BDD5EA')
#     df=copy(x)[
#         ,column_to_plot:=column
#     ]
#     p = ggplot(df, aes(x=overlap, y=column_to_plot,fill=overlap)) +
#     geom_violin(trim=T,scale = "width")+
#     geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
#     ylab(ylab)+xlab(xlab)+
#     scale_fill_manual(values=tadass_palette)+
#     facet_wrap(DA~.,ncol=3)+
#     stat_compare_means(
#         method = "wilcox.test",
#         label.y = max(df$column_to_plot+0.1),
#         size=5
#         )+
#     theme_classic()+
#     theme(
#         legend.position = 'none',
#         strip.background = element_blank(),
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
#     return(p)
# }

# ## plot cons scores
# pdf(paste(outplot_dir,'homer_ctcf_peaks_PhastCons_score.pdf',sep=''),width= 7,height = 5)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$avg_phastcons,ylab='PhastCons score',xlab=' ') 
# dev.off()

# ## plot Deepbind CTCF affinities
# pdf(paste(outplot_dir,'homer_ctcf_peaks_Deepbind_affinities.pdf',sep=''),width= 7,height = 5)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
# dev.off()





