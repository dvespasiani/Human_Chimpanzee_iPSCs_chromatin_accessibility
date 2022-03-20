# ## use this script to analyse the predicted CTCF binding profile 
# ## in particular compare the CTCF binding affinities between peaks containing the ctcf bs and all the genome-wide bs fo ctcf
# ## then do the analysis by comparing the affinities btwn tad associated and non-tad associated ctcf
# ## then look if there is an enrichment of CTCF bs near DC TADs 

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

# ctcf_dir <- './output/homer/homer_output/ctcf/'
# outplot_dir <- create_dir(plot_dir,'ctcf')
# peakDir = './output/DA/peaks/'
# tads_dir =  './output/TADs/'
# tmp_files_dir = './output/temp_files/'

# ##-----------------------
# ## read deepBind results
# ##-----------------------
# ctcf_affinities <- list.files(
#     paste(ctcf_dir,'deepbind_predictions/',sep=''),recursive=F,full.names=T,pattern='affinity')%>%
#     lapply(
#         function(x){
#             x <- fread(x,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(2,3),new=c('sequence','prediction'))
#         }

# ) 
# names(ctcf_affinities) =  species_names

# ## because homer returns multiple CTCFs for a single peak which are offset by few bp from each other 
# ## take the sequence with the highest absolute binding affinity for each single peak

# ## read ctcf peaks and get coordinates ctcf bs not peak
# ctcf_peak_coords <- list.files(paste(ctcf_dir,'ctcf_peaks/',sep=''),recursive=F,full.names=T)%>%lapply(function(x)fread(x,sep='\t',header=T)%>%unique())
# names(ctcf_peak_coords) = species_names

# peaks_with_ctcf = purrr::map2(ctcf_affinities,ctcf_peak_coords,function(x,y){
#     z<-inner_join(x,y,by=c('peakID','sequence'))
#     z <- z[
#         ,.SD[which.max(abs(prediction))], by=.(peakID)
#         ][
#             ,peak_center:= round((end-start)/2,0)
#             ][
#                 ,ctcf_peak_start:=start+(peak_center+Offset) # this returns the coord of the start of ctcf peak
#                 ][
#                     ,ctcf_peak_end:=ctcf_peak_start+nchar(sequence)  # this returns the coord of the end of ctcf peak (based on length of motif sequence)
#                     ][
#                         ,c('peak_center','sequence'):=NULL
#                         ]
#     return(z)
#     }
# )

# # ##-----------------------------------------------
# # ## Compare the predicted binding affinities 
# # ## and Homer MotifScores differ between species
# # ##-----------------------------------------------
# # all_ctcf_peaks = copy(ctcf_peaks)
# # all_ctcf_peaks = Map(mutate,all_ctcf_peaks,species=species_names)%>%rbindlist()

# # stat_test_ctcf_affinities = 
# # compare_means(
# #   prediction~species,
# #   all_ctcf_peaks,
# #   method = "wilcox.test",
# #   paired = FALSE,
# #   ref.group = 'common',
# #   p.adjust.method = "fdr"
# # ) %>% as.data.table()%>%setnames(old='group2',new='species')

# # # stat_test_ctcf_affinities =Map(mutate,stat_test_ctcf_affinities,file=species_names[-2])%>%rbindlist()

# # stat_test_MotifScore = 
# # compare_means(
# #   MotifScore~species,
# #   all_ctcf_peaks,
# #   method = "wilcox.test",
# #   paired = FALSE,
# #   ref.group = 'common',
# #   p.adjust.method = "fdr"
# # ) %>% as.data.table()%>%setnames(old='group2',new='species')

# # # stat_test_MotifScore =Map(mutate,stat_test_MotifScore,file=species_names[-2])%>%rbindlist()

# # # all_ctcf_peaks = Map(mutate,all_ctcf_peaks,file=species_names[-2])%>%rbindlist()

# # ##plot
# # pdf(paste(plot_dir,'ctcf_binding_affinities.pdf',sep=''),width=10,height = 7)
# # ggplot(all_ctcf_peaks, aes(x=species, y=prediction,fill=species)) +
# # geom_boxplot(notch=T) +ylab('binding affinity prediction score')+
# # geom_text(data=stat_test_ctcf_affinities, aes(x=species, y=max(all_ctcf_peaks$prediction+0.01), label=p.signif), col='black', size=7)+
# # theme(
# #     axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
# #     )
# # dev.off()

# # pdf(paste(plot_dir,'ctcf_motifscores.pdf',sep=''),width=10,height = 7)
# # ggplot(all_ctcf_peaks, aes(x=species, y=MotifScore,fill=species)) +
# # geom_boxplot(notch=T) +ylab('Homer predicted MotifScore')+
# # geom_text(data=stat_test_ctcf_affinities, aes(x=species, y=max(all_ctcf_peaks$MotifScore+0.01), label=p.signif), col='black', size=7)+
# # theme(
# #     axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
# #     )
# # dev.off()

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
# boundary =  300000
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
# ctcf_boundary_overlap = purrr::map2(peaks_with_ctcf,tad_boundaries,function(x,y)foverlaps(x,y,type='within')%>%na.omit())
# ctcf_boundary_overlap = lapply(ctcf_boundary_overlap,
#     function(x)x=x[
#         ,distance_ctcf_tad:= ifelse(boundary=='upstream',ctcf_peak_end-end,ctcf_peak_start-start)
#     ]
# )

# ## plot distribution distances between CTCF and TAD boundaries
# ## if peak is associated with multiple tads keep closest one
# distance_ctcf_boundaries <- copy(ctcf_boundary_overlap)%>%
# lapply(function(x)x=x[
#     ,c('peakID','prediction','species','distance_ctcf_tad')
#     ][
#         , .SD[which.min(abs(distance_ctcf_tad))], by=.(peakID,species)
#         ][
#             ,log10dist:=log10(abs(distance_ctcf_tad))
#             ][
#                 ,log10dist:=ifelse(distance_ctcf_tad<0, -log10dist,log10dist)
#                 ][
#                     ,associated:=ifelse(abs(distance_ctcf_tad)<=50000,'tad_associated','non_associated')
#                     ]
# )
# distance_ctcf_boundaries <- Map(mutate,distance_ctcf_boundaries,peak_species=species_names)%>%rbindlist()

# pdf(paste(outplot_dir,'ctcf_tad_distances.pdf',sep=''),width=10,height = 7)
# ggplot(distance_ctcf_boundaries, aes(x=distance_ctcf_tad,fill=peak_species)) +
# geom_density(alpha=0.5) +ylab('Density CTCF peaks')+ xlab('CTCF-TAD Distance')+
# theme(
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
# dev.off()

# ## get some % CTCF peaks within tab boundaries 
# numb_ctcf_peaks = copy(peaks_with_ctcf)%>%lapply(function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
# numb_ctcf_peaks_tad_overlap_1mb = copy(ctcf_boundary_overlap)%>%lapply(function(x)x[,c('seqnames','i.start','i.end')]%>%unique()%>%nrow())
# numb_ctcf_peaks_tad_overlap_50kb = copy(ctcf_boundary_overlap)%>%lapply(function(x)x[abs(distance_ctcf_tad)<=50000][,c('seqnames','i.start','i.end')]%>%unique()%>%nrow())

# prop_ctcf_peaks_tad_overlap_1mb = purrr::map2(numb_ctcf_peaks_tad_overlap_1mb,numb_ctcf_peaks,function(x,y)x/y*100)
# prop_ctcf_peaks_tad_overlap_50kb = purrr::map2(numb_ctcf_peaks_tad_overlap_50kb,numb_ctcf_peaks,function(x,y)x/y*100)

# ##------------------------------------------------
# ## Get conservation score 
# ## of common/species specific CTCF peaks near TADs 
# ##------------------------------------------------
# ctcf_peaks_conScore = copy(peaks_with_ctcf)%>%lapply(function(x)x[,c('seqnames','ctcf_peak_start','ctcf_peak_end','peakID')]%>%unique())

# ## get average conservation score
# phastCons = phastCons7way.UCSC.hg38
# ctcf_peaks_conScore = lapply(
#     ctcf_peaks_conScore,function(x){ 
#     x<-x%>%makeGRangesFromDataFrame(keep.extra.columns=T)
#     x<-gscores(phastCons, x)%>%as.data.table()
#     return(x)
#     }
# )

# ctcf_boundary_overlap_50kb <- lapply(ctcf_boundary_overlap,function(x)x=x[abs(distance_ctcf_tad)<=50000])

# ctcf_peaks_conScore <- purrr::map2(
#     ctcf_peaks_conScore,ctcf_boundary_overlap_50kb,function(x,y)
#     x=x[
#         ,tad_overlap:=ifelse(peakID %in% y$peakID,'yes','no')
#         ][
#             , default := nafill(default,fill=0)
#             ]
# )

# ctcf_peaks_conScore <- Map(mutate,ctcf_peaks_conScore,peak_species=names(ctcf_peaks_conScore))%>%rbindlist()

# ## plot function
# plot_ctcf_tad = function(x,column,ylab,xlab){
#     df=copy(x)[
#         ,column_to_plot:=column
#     ]
#     p = ggplot(df, aes(x=tad_overlap, y=column_to_plot,fill=tad_overlap)) +
#     geom_violin(trim=T,scale = "width")+
#     geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+    ylab(ylab)+xlab(xlab)+
#     facet_wrap(peak_species~.,ncol=3)+
#     stat_compare_means(
#         method = "wilcox.test",
#         # label = "p.signif",
#         label.y = max(df$column_to_plot+0.1),
#         size=5
#         )+
#     theme(
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
#     return(p)
# }

# ## plot cons scores
# pdf(paste(outplot_dir,'ctcf_peaks_PhastCons_score.pdf',sep=''),width= 10,height = 8)
# plot_ctcf_tad(ctcf_peaks_conScore,ctcf_peaks_conScore$default,ylab='PhastCons score',xlab=' ') 
# dev.off()

# ## now plot motif binding affinities for CTCFs within 50kb
# ctcf_peak_affinities <- purrr::map2(
#     peaks_with_ctcf,ctcf_boundary_overlap_50kb,function(x,y)
#     x=x[
#         ,tad_overlap:=ifelse(peakID %in% y$peakID,'yes','no')
#         ]
# )
# ctcf_peak_affinities <- Map(mutate,ctcf_peak_affinities,species=names(ctcf_peak_affinities))%>%rbindlist()

# ## plot Deepbind CTCF affinities
# pdf(paste(outplot_dir,'ctcf_peaks_Deepbind_affinities.pdf',sep=''),width= 10,height = 8)
# plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
# dev.off()

# # ## plot Motifscores
# # pdf(paste(plot_dir,'ctcf_peaks_Motifscores.pdf',sep=''),width= 10,height = 8)
# # plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$MotifScore,ylab='CTCF HOMER motif score',xlab=' ') 
# # dev.off()

