# ## script used to analyses the distribution of DA peaks nearby lncRNA TSSs
# library(dplyr)
# library(data.table)
# library(magrittr)
# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome.Ptroglodytes.UCSC.panTro5)
# library(ggthemes)
# library(ggplot2)
# library(ggpubr)

# ## read all and DA peaks
# ## and see if the DA peaks have different GC/CpG content 

# options(width=150)
# setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

# scripts_dir = './post_processing_analyses/scripts/'
# source(paste(scripts_dir,'utils.R',sep=''))

# peakDir = './post_processing_analyses/output/DA/peaks/'
# plot_dir = './post_processing_analyses/output/plots/GC/'

# human_peaks_hg38 = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
# chimp_peaks_pantro5 = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
# common_peaks = read_da_peaks('common_regions.txt',c('significant','non_significant'))

# ##--------------------------------------------
# ## First calculate CpG content of these peaks
# ##--------------------------------------------
# ## expectations are that common peaks should have higher CpG content than species-specific peaks
# ## if these are enriched for promoter-like elements
# get_cpg =function(peaks,species){
#     cpg=copy(peaks)%>%makeGRangesFromDataFrame()%>%Repitools::cpgDensityCalc(organism=species)
#     peaks = peaks[,cpg_content:=cpg]
#     return(peaks)
# }

# chimp_peaks_pantro5=get_cpg(chimp_peaks_pantro5,Ptroglodytes)
# common_peaks=get_cpg(common_peaks,Hsapiens)
# human_peaks_hg38=get_cpg(human_peaks_hg38,Hsapiens)

# common_peaks = common_peaks[!cpg_content>6000] ## single outlier
# all_peaks = list(chimp_peaks_pantro5,common_peaks,human_peaks_hg38)

# ## plot CpG densities
# all_peaks_cpg = copy(all_peaks)%>%lapply(function(x)x=x[,c('cpg_content','species')])%>%rbindlist()

# pdf(paste0(plot_dir,'peaks_CpG_content.pdf',sep=''),width = 7, height = 7)
# ggplot(all_peaks_cpg,aes(x=species,y=cpg_content,fill=species))+
# geom_boxplot(notch=T)+
# xlab('')+ylab('Raw count CpG')+
# stat_compare_means(
#   method = "wilcox.test",
#   label.y = (max(all_peaks_cpg$cpg_content)+0.5),
#   size=5,
#   ref.group='common'
#   )+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()

# ##-------------------------
# ## Now look at GC content
# ##-------------------------
# gc_content = function(x,genome){
#     df = copy(x)%>%as.data.table()
#     df=df[
#         ,sequence:=as.character(getSeq(genome, seqnames,start, end))
#         ][
#             !sequence%like%'N' ## do this because u need to remove sequences with Ns 
#             ][
#             ,gc_content:=letterFrequency(DNAStringSet(sequence),"GC")/width(DNAStringSet(sequence))
#             ][
#                 ,sequence:=NULL
#                 ]

#     return(df)
# }

# chimp_peaks_pantro5_gc=gc_content(chimp_peaks_pantro5,Ptroglodytes)
# common_peaks_hg38_gc=gc_content(common_peaks,Hsapiens)
# human_peaks_hg38_gc=gc_content(human_peaks_hg38,Hsapiens)

# ## plot CpG densities
# all_peaks_gc = rbind(
#     chimp_peaks_pantro5_gc[,c('gc_content','species')],
#     common_peaks_hg38_gc[,c('gc_content','species')],
#     human_peaks_hg38_gc[,c('gc_content','species')]
# )

# pdf(paste0(plot_dir,'peaks_GC_content.pdf',sep=''),width = 7, height = 7)
# ggplot(all_peaks_gc,aes(x=species,y=gc_content,fill=species))+
# geom_boxplot(notch=T)+
# xlab('')+ylab('Fraction GC over peak width')+
# stat_compare_means(
#   method = "wilcox.test",
#   label.y = (max(all_peaks_gc$gc_content)+0.05),
#   size=5,
#   ref.group='common'
#   )+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()



# # ## to search for overlaps with lncRNAs liftover everything to hg38
# # ## then for the peaks in common get the common regions
# # chain_path <- './data/LiftOver_chains/'
# # chimp_peaks_hg38 = convert_coord(chimp_peaks_pantro5,'panTro5ToHg38.over.chain')

# # chimp_common_peaks_pantro5 = copy(common_peaks)%>%dplyr::select(c(contains('chimp'),'significant'))
# # colnames(chimp_common_peaks_pantro5)[1:4] = c('peakID',range_keys)

# # chimp_common_peaks_hg38 = convert_coord(chimp_common_peaks_pantro5,'panTro5ToHg38.over.chain')
# # setkeyv(chimp_common_peaks_hg38,range_keys)

# # human_common_peaks_hg38 = copy(common_peaks)%>%dplyr::select(c(contains('human'),'significant'))
# # colnames(human_common_peaks_hg38)[1:4] = c('peakID',range_keys)
# # setkeyv(human_common_peaks_hg38,range_keys)

# # ## now get the common regions
# # common_regions_hg38 = foverlaps(human_common_peaks_hg38,chimp_common_peaks_hg38,type='any')
# # common_regions_hg38 = common_regions_hg38[
# #     ,start:=ifelse(start <= i.start,start,i.start)
# #     ][
# #     ,end:=ifelse(end >= i.end,end,i.end)
# #     ][
# #         ,c(paste('i',range_keys[-1],sep='.')):=NULL
# # ]


# # all_peaks = list(common_regions_hg38,chimp_peaks_hg38,human_peaks_hg38)%>%
# #     lapply(function(x)
# #         x=x[
# #             ,c(..range_keys,'peakID','significant')
# #             ]
# # )
# # names(all_peaks) = c('common_regions_hg38','chimp_peaks_hg38','human_peaks_hg38')

# # lapply(all_peaks,function(x)setkeyv(x,range_keys))


# # ## first see if there is an enrichment nearby lncRNA TSSs
# # ## read genecode lncRNAs (hg38)
# # ## GENECODE/NONCODE (for this latter perhaps take 2kb up/downstream tss and see if there are any peaks nearby)
# # ## then also look at the pattern of binding around these TSS if there is an interesting significant excess of DA regions.
# # ## do this on human perspective, i.e. convert pantro coordinates to humans (and vice-versa) and check GC content

# # lncRNAs = fread('./post_processing_analyses/data/genecodev36_lncRNA_hg38',sep='\t',header=T)

# # lncRNA_upstream_reg_region =  copy(lncRNAs)[
# #     ,end:=start
# #     ][
# #         ,start:=start-5000
# #         ][
# #             ,location:='upstream'
# # ]
# # lncRNA_downstream_reg_region =  copy(lncRNAs)[
# #     ,start:=end
# #     ][
# #         ,end:=end+5000
# #         ][
# #             ,location:='downstream'
# # ]

# # lncRNA_reg_region = rbind(lncRNA_upstream_reg_region,lncRNA_downstream_reg_region)%>%setorderv('id')

# # setkeyv(lncRNA_reg_region,range_keys)

# # ## overlaps only with hg38 coordinates and then lift-back to get pantro5 coordinates
# # overlap =function(x,y){
# #     df = foverlaps(x,y,type='any')[
# #     ,c(range_keys[-1]):=NULL]%>%na.omit()%>%setorderv('id')%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
# #     return(df)
# # }
# # all_peaks_lncRNAoverlap = lapply(all_peaks,function(x)overlap(x,lncRNA_reg_region))

# # # lapply(all_peaks_lncRNAoverlap,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
# # # lapply(all_peaks_lncRNAoverlap,function(x)x[,c(2)]%>%unique()%>%nrow())

# # ##-----------------------------------------------
# # ## Calculate odds ratio between da/non-da peaks
# # ## near or not lncRNA TSSs 
# # ## for all common and species-specific peaks
# # ##-----------------------------------------------
# # count_peaks = function(x,type){
# #     count =copy(x)
# #     count=count[
# #         ,numb_peaks:=.N,by=.(significant)
# #         ][
# #             ,peak_type:=type
# #         ][
# #         ,c("peak_type","numb_peaks",'significant')
# #         ]%>%unique()
# #     return(count)
# # }

# # count_all_peaks = lapply(all_peaks,function(y)count_peaks(y,'all_peaks'))


# # count_all_lncrna_peaks = lapply(all_peaks_lncRNAoverlap,function(y)count_peaks(y,'lncrna_peaks'))
# # # count_all_lncrna_peaks = Map(mutate,count_all_lncrna_peaks,'file'=names(count_all_lncrna_peaks))

# # ## test if species-specific peaks are more enriched near lncRNAs than common peaks
# # ## use fisher test
# # fisher_test_matrix = purrr::map2(count_all_lncrna_peaks,count_all_peaks,merge,by=c('significant'))

# # # human_specific = rbindlist(fisher_test_matrix[c(1,4)])
# # # chimp_specific = rbindlist(fisher_test_matrix[c(1,2)])
# # # chimp_hg38_specific = rbindlist(fisher_test_matrix[c(1,3)])

# # fisher_test_matrix = fisher_test_matrix%>%
# #     lapply(function(x)x=x[
# #         ,non_in_lncrna_peaks:=numb_peaks.y-numb_peaks.x
# #         ][
# #             ,c(1,2,4,5):=NULL
# #             ]%>%setnames(old='numb_peaks.x',new='in_lncrna_peaks')%>%as.matrix()
# # )
# # ## ps matrix rownames are: c(ns,s)

# # ## there is significant enrichment of species-specific peaks near lncRNA TSSs
# # fisher_test_results = lapply(fisher_test_matrix,function(x)fisher.test(x))

# # ##-----------------------------------------------------------------
# # ## compare GC content of specie-specific vs common DA peaks 
# # ## near and not near lncRNA TSSs
# # ##-----------------------------------------------------------------

# # gc_content = function(x,genome){
# #     df = copy(x)%>%as.data.table()
# #     df=df[
# #         ,sequence:=as.character(getSeq(genome, seqnames,start, end))
# #         ][
# #             !sequence%like%'N' ## do this because u need to remove sequences with Ns 
# #             ][
# #             ,gc_content:=letterFrequency(DNAStringSet(sequence),"GC")/width(DNAStringSet(sequence))
# #             ][
# #                 ,sequence:=NULL
# #                 ]

# #     return(df)
# # }


# # ## non lncRNA peaks
# # non_lncrna_peaks = purrr::map2(all_peaks,all_peaks_lncRNAoverlap,anti_join,by=c(range_keys))
# # non_lncrna_peaks_hg38_gc = lapply(non_lncrna_peaks,function(y)y%>%gc_content(Hsapiens))
# # # non_lncrna_peaks_pantro5_gc = gc_content(non_lncrna_peaks[[2]], Ptroglodytes)

# # # non_lncrna_peaks_gc = list(non_lncrna_peaks_hg38_gc[[1]],non_lncrna_peaks_pantro5_gc,non_lncrna_peaks_hg38_gc[[2]],non_lncrna_peaks_hg38_gc[[3]])
# # non_lncrna_peaks_gc = Map(mutate,non_lncrna_peaks_hg38_gc,'group'=names(all_peaks))%>%rbindlist()
# # non_lncrna_peaks_gc = non_lncrna_peaks_gc[,type:='non_lncrna_peak']

# # ## lncrna peaks
# # lncrna_peaks_hg38_gc = lapply(all_peaks_lncRNAoverlap,function(y)y%>%gc_content(Hsapiens))
# # # lncrna_peaks_pantro5_gc = gc_content(all_peaks_lncRNAoverlap[[2]], Ptroglodytes)

# # # lncrna_peaks_gc = list(lncrna_peaks_hg38_gc[[1]],lncrna_peaks_pantro5_gc,lncrna_peaks_hg38_gc[[2]],lncrna_peaks_hg38_gc[[3]])
# # lncrna_peaks_gc = Map(mutate,lncrna_peaks_hg38_gc,'group'=names(all_peaks))%>%rbindlist()
# # lncrna_peaks_gc = lncrna_peaks_gc[,type:='lncrna_peak']


# # final_all_peaks_gc = rbind(non_lncrna_peaks_gc,lncrna_peaks_gc[,c(..range_keys,'peakID','significant','gc_content','group','type')])
# # final_all_peaks_gc = final_all_peaks_gc[
# #     ,group:=paste(group,type,sep='.')
# #     ][
# #         significant=='significant'
# # ]



# # test_comparisons = list(
# #     c("common_regions_hg38.lncrna_peak", "common_regions_hg38.non_lncrna_peak"),
# #    c("human_peaks_hg38.lncrna_peak" , "human_peaks_hg38.non_lncrna_peak"),
# #    c("chimp_peaks_hg38.lncrna_peak" , "chimp_peaks_hg38.non_lncrna_peak")
# # ) 

# # # plot
# # pdf(paste(plot_dir,'lncrna_vs_nonlncrna_peaks_gc_content.pdf',sep=''),width=7,height = 7)
# # ggplot(final_all_peaks_gc, aes(x=group, y=gc_content,fill=type)) +
# #   geom_boxplot(notch=T) +
# #   geom_signif(
# #     test = "t.test",
# #     comparisons = test_comparisons,
# #     map_signif_level = TRUE,
# #     textsize = 6)+
# #      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
# # dev.off()


# # # ## overall gc content species-specific vs common regions
# # # all_peaks_hg38_gc = lapply(all_peaks[c(1,3:4)],function(y)y%>%gc_content(Hsapiens))
# # # all_peaks_pantro5_gc = gc_content(all_peaks[[2]], Ptroglodytes)

# # # all_peaks_gc = list(all_peaks_hg38_gc[[1]],all_peaks_pantro5_gc,all_peaks_hg38_gc[[2]],all_peaks_hg38_gc[[3]])
# # # all_peaks_gc = Map(mutate,all_peaks_gc,'group'=names(all_peaks))%>%rbindlist()

# # # all_peaks_comparisons = list(
# # #     c("common", "human"),
# # #     c("chimp_hg38", "human"), 
# # #     c("chimp_pantro5", "human"),
# # #     c("common", "chimp_hg38"),
# # #     c("common", "chimp_pantro5")
# # # ) 

# # # # plot
# # # pdf(paste(plot_dir,'all_peaks_gc_content.pdf',sep=''),width=7,height = 7)
# # # ggplot(all_peaks_gc, aes(x=group, y=gc_content,fill=group)) +
# # #   geom_boxplot(notch=T) +
# # #   geom_signif(
# # #     test = "t.test",
# # #     comparisons = all_peaks_comparisons,
# # #     map_signif_level = TRUE,
# # #     textsize = 6)+
# # #      theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
# # # dev.off()


# # ## look proportion promoter vs enhancers at these peaks
# # library(annotatr)
# # # Build the annotations
# # annots = c('hg38_basicgenes', 'hg38_genes_intergenic')
# # annotations = build_annotations(genome = 'hg38', annotations = annots)

# # ## annotate significant DA peaks near lcnRNAs
# # lncrna_peaks_hg38_annot = lapply(lncrna_peaks_hg38_gc,function(x)x=x[
# #     significant=='significant'
# #     ]%>%makeGRangesFromDataFrame()
# # )
# # lncrna_peaks_hg38_annot = lapply(lncrna_peaks_hg38_annot,function(x)
# #     x = annotate_regions(
# #     regions = x,
# #     annotations = annotations,
# #     ignore.strand = TRUE,
# #     quiet = FALSE)%>%as.data.table()
# # )

# # lncrna_peaks_hg38_annot = lapply(lncrna_peaks_hg38_annot,function(x)
# # x=x[
# #     ,c(..range_keys,'annot.type')
# #     ][
# #         ,tot_peaks:=.N 
# #         ][
# #           ,tot_peaks_per_annot:=.N,by=.(annot.type)
# #           ][
# #             ,fraction:=tot_peaks_per_annot/tot_peaks
# #             ][
# #                 ,genomic_element:=plyr::revalue(annot.type,c(
# #                     'hg38_genes_exons'='Exons',
# #                     'hg38_genes_introns'='Introns',
# #                     'hg38_genes_intergenic'='Intergenic',
# #                     'hg38_genes_promoters'='Promoters',
# #                     'hg38_genes_1to5kb'='Enhancers',
# #                     'hg38_genes_3UTRs'='3UTRs',
# #                     'hg38_genes_5UTRs'='5UTRs'))
# #                     ][
# #                         ,annot.type:=NULL
# #                         ]

# # )
# # lncrna_peaks_hg38_annot = Map(mutate,lncrna_peaks_hg38_annot,'group'=names(lncrna_peaks_hg38_annot))%>%rbindlist()


# # pdf(paste(plot_dir,'annotation_lncrna_peaks.pdf',sep=''),width=7,height = 7)
# # ggplot(lncrna_peaks_hg38_annot[,c('genomic_element','group','fraction')]%>%unique(), aes(x=reorder(genomic_element,-fraction),
# #                               y=fraction,fill=group))+
# #   geom_bar(stat="identity",position ='dodge')+
# #   ylab('\n Proportion peaks near lncRNA TSSs \n')+
# #   xlab('\n  \n')+
# #   theme(panel.background =element_rect(fill = 'white', colour = 'white'),
# #         panel.grid.minor = element_blank(),
# #         panel.grid.major = element_blank(),
# #         legend.text = element_text(),
# #         legend.title = element_text(),
# #         legend.margin = margin(c(0.5, 2, 8, 25)),
# #         legend.spacing.x = unit(0.5, 'cm'),
# #         axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1),
# #         axis.text.y = element_text(),
# #         axis.title.y = element_text(hjust=0.5),
# #         axis.text=element_text(),
# #         axis.line = element_line(color = "black",size = 0.5, linetype = "solid"))
# # dev.off()


