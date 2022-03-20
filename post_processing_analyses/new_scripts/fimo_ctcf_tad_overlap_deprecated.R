
# library(dplyr)
# library(data.table)
# library(magrittr)
# library(GenomicRanges)
# library(ggthemes)
# library(ggplot2)
# library(ggthemes)
# library(ggpubr)
# library(phastCons7way.UCSC.hg38)
# library(GenomicScores)

# options(width=150)
# setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

# scripts_dir = './scripts/'
# source(paste(scripts_dir,'utils.R',sep=''))

# # ctcf_dir <- './output/homer/homer_output/ctcf/'
# outplot_dir <- create_dir(plot_dir,'ctcf')
# peakDir = './output/DA/peaks/'
# tads_dir =  './output/TADs/'
# files_dir <- './output/files/'

# ##-----------------------
# ## read deepBind results
# ##-----------------------
# ## read deepbind input with all the info 
# deepbind_output_file <- dir(paste(files_dir,'deepbind_predictions',sep=''),recursive=F,full.names=T,pattern='fimo')
# ctcf_peak_file <- dir(paste(files_dir,'fimo_ctcf',sep=''),recursive=F,full.names=T,pattern='fimo')

# ctcf_peaks <- fread(ctcf_peak_file,header=T,sep='\t')
# ctcf_affinities <- fread(deepbind_output_file,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(2,3),new=c('sequence','prediction'))

# ctcf_peaks <- ctcf_peaks[
#     ctcf_affinities,on=c('peakID','sequence'),nomatch=0
#     ][
#         ,peak_species:=ifelse(peak_species=='non_da',"common",ifelse(peak_species=='da_chimp','chimp','human'))
# ]
# setkeyv(ctcf_peaks,range_keys)

# ctcf_peaks <- setorderv(ctcf_peaks,'peak_species')%>%split(by='peak_species')

# ##-------------------------------------------------
# ## look for overlaps between ctcf at TAD boundaries
# ##--------------------------------------------------
# orth_tads <- read_tads('human_chimp_orth_tads.txt')
# chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL][species=='Chimp']
# human_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

# chimp_tads_hg38 <- copy(orth_tads)[species=='Chimp']
# colnames(chimp_tads_hg38)[4:6] = range_keys
# chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'species'),nomatch=0][
#   ,c(range_keys):=NULL
# ]%>%setnames(old=c(1:3),new=range_keys)

# common_tads_hg38 <- copy(human_tads_hg38)[species=='Common']

# all_tads <- list(chimp_tads_hg38,common_tads_hg38,human_tads_hg38[species!='Common'])

# all_tads <- lapply(all_tads,function(x){
#   x <- setorderv(x,range_keys,1)
#   x<-x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]
#   setkeyv(x,range_keys)
#   return(x)
#   }
# ) 
# names(all_tads) = species_names

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

# ## look at ctcf overlap at tad boundaries
# ## add logFC info to ctcf peaks 
# da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
# da_results <- fread(da_file,sep='\t',header=T,select=c('logFC','peakID'))

# common_tads = copy(tad_boundaries[[2]])
# ctcf_peaks_da = rbindlist(copy(ctcf_peaks))[da_results,on='peakID',nomatch=0]

# ctcf_boundary_overlap <-foverlaps(ctcf_peaks_da,common_tads,type='within')%>%na.omit()
# ctcf_boundary_overlap <-ctcf_boundary_overlap[
#     ,distance_ctcf_tad:= ifelse(boundary=='upstream',i.end-end,i.start-start)
#     ][
#         , .SD[which.min(abs(distance_ctcf_tad))], by=.(peakID,species)
#         ][
#             ,c('start','end','species','tadID','boundary'):=NULL
# ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))


# ## check densitity CTCF associated per 10kb binned distance 
# ctcf_boundary_overlap <- ctcf_boundary_overlap[
#     ,binned_abs_dist:=plyr::round_any(abs(distance_ctcf_tad), 10000)
#     ][
#         ,numb_peaks:=.N,by=.(DA)
#         ][
#             ,numb_peaks_binned_dist:=.N,by=.(DA,binned_abs_dist)
#             ][
#                 ,prop:=numb_peaks_binned_dist/numb_peaks
# ]

# ## get some % CTCF peaks within tab boundaries 
# numb_ctcfs_associated <- copy(ctcf_boundary_overlap)[,numb_ctcfs:=.N,by=.(DA)][,c('numb_ctcfs','DA')]%>%unique()%>%setorderv('DA',1)%>%split(by='DA')
# tot_ctcfs <- copy(ctcf_peaks_da)[,c('DA')]%>%split(by='DA')
# prop_ctcfs_associated <- purrr::map2(tot_ctcfs,numb_ctcfs_associated,function(x,y){
#     prop<-round((y$numb_ctcfs/nrow(x))*100,2)
#     return(prop)
# })

# ## enrichment non-DA at common TAD boundaries?
# ## correct 2x2 table structure
# matrix <- as.matrix(data.table(col1=c(
#     numb_ctcfs_associated[[2]]$numb_ctcfs,
#     numb_ctcfs_associated[[1]]$numb_ctcfs
#     ),
#     col2=c(
#     nrow(tot_ctcfs[[2]])-numb_ctcfs_associated[[2]]$numb_ctcfs,
#     nrow(tot_ctcfs[[1]])-numb_ctcfs_associated[[1]]$numb_ctcfs
#     )
# ))
# fisher.test(matrix)
# # p-value = 0.22
# # odds ratio 1.16322 

# ## plot ctcf densities
# pdf(paste(outplot_dir,'fimo_ctcf_density_binned_distance.pdf',sep=''),width= 10,height = 8)
# ggplot(ctcf_boundary_overlap,aes(x=binned_abs_dist,y=prop,fill=DA))+
# geom_bar(stat='identity',position=position_dodge())+
# scale_fill_manual(values=da_palette)+
# xlab('TAD boundary size 10kb resolution')+ylab('Proportion CTCF overlapping TAD boundary')+
# theme_bw()+
# theme()
# dev.off()

# ## compare logFC of DA peaks with CTCF associated to TAD vs those not 
# ctcf_peaks_da <- ctcf_peaks_da[
#             ,tad_ass:=ifelse(peakID %in% ctcf_boundary_overlap$peakID,'ass','no_ass')
# ]

# pdf(paste(outplot_dir,'fimo_ctcf_da_peak_logFC_tad_overlap.pdf',sep=''),width= 10,height = 8)
# df <- copy(ctcf_peaks_da)[DA=='da']

# ggplot(df,aes(x=tad_ass,y=abs(logFC),fill=tad_ass))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
# xlab('TAD boundary associations')+ylab('absolute logFC')+
# geom_hline(yintercept=median(abs(df[tad_ass=='ass']$logFC)),linetype='dashed')+
# stat_compare_means(
#     method = "wilcox.test",
#     size=5,
#     ref.group = 'no_ass'
#     )+
# theme_bw()+
# theme(
#     legend.position ='none',
#     axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#     )
# dev.off()

# ##------------------------------------------------------------------------------------------------
# ## perform comparison between CTCF associated and non to TAD boundaries, merging DA and non-DA 
# ##------------------------------------------------------------------------------------------------
# ## Get conservation score 
# phastCons <- phastCons7way.UCSC.hg38

# ctcf_peaks_conScore <- copy(ctcf_peaks_da)%>%makeGRangesFromDataFrame(keep.extra.columns=T)
# ctcf_peaks_conScore <- gscores(phastCons,ctcf_peaks_conScore)%>%as.data.table()

# ctcf_peaks_conScore <- ctcf_peaks_conScore[,default:=nafill(default,fill=0)]

# ## plot function
# plot_ctcf_tad = function(x,column,ylab,xlab){
#     df=copy(x)[
#         ,column_to_plot:=column
#     ]
#     p = ggplot(df, aes(x=tad_ass, y=column_to_plot,fill=tad_ass)) +
#     geom_violin(trim=T,scale = "width")+
#     geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
#     ylab(ylab)+xlab(xlab)+
#     geom_hline(yintercept=median(abs(df[tad_ass=='ass']$column_to_plot)),linetype='dashed')+
#     stat_compare_means(
#         method = "wilcox.test",
#         label.y = max(df$column_to_plot+0.1),
#         size=5
#         )+
#     theme_classic()+
#     theme(
#         legend.position = 'none',
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
#     return(p)
# }

# ## plot cons scores
# pdf(paste(outplot_dir,'fimo_ctcf_peaks_PhastCons_score.pdf',sep=''),width= 7,height = 7)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$default,ylab='PhastCons score',xlab=' ') 
# dev.off()

# ## plot Deepbind CTCF affinities
# pdf(paste(outplot_dir,'fimo_ctcf_peaks_Deepbind_affinities.pdf',sep=''),width= 7,height = 7)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
# dev.off()

# # ## plot Motifscores
# # pdf(paste(plot_dir,'ctcf_peaks_Motifscores.pdf',sep=''),width= 10,height = 8)
# # plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$MotifScore,ylab='CTCF HOMER motif score',xlab=' ') 
# # dev.off()





