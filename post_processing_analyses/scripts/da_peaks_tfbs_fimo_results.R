## script used to analyse homer tfbs known motifs results and get enrichment
## script to calculate tf enrichment after fimo meme analysis
## enrichment for tfs across chrom state in da vs non da peaks (make heatmap)
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)


options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))


da_peaks_dir = './output/DA/peaks/'
fimo_dir = './output/fimo_meme'
plot_dir = './output/plots/tfbs/'

range_keys = c('seqnames','start','end')

read_tfbs_results = function(dir,pattern){
  df = list.files(dir,recursive=T,full.names=T,pattern=pattern)%>%
    lapply(function(x)x=fread(x,sep='\t',header=T))
  names(df) =  gsub("\\/.*","",list.files(dir,recursive=T,full.names=F,pattern=pattern))
  return(df)
}


fimo_results = read_tfbs_results(fimo_dir,'fimo.tsv')

fimo_results = lapply(fimo_results,function(x)setnames(x,old=c('p-value','q-value'),new=c('p_value','q_value'))[
    q_value<0.01
    ][
      ,log10_adjP:=-log10(q_value)
    ][
      ,c('q_value','p_value'):=NULL
      ][
        , .SD[which.max(log10_adjP)], by=sequence_name
        ]
)

## add cluster info
motif_cluster = fread('/data/projects/punim0586/dvespasiani/Files/Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

fimo_cluster = copy(fimo_results)
fimo_cluster = lapply(fimo_cluster,function(x)x=x[
  motif_cluster,on='motif_id',nomatch=0
  ][
    ,c('motif_id','log10_adjP','Name')
    ][
    ,motif_id:= gsub("\\_.*","",motif_id)
    ][
      ,.SD[which.max(log10_adjP)], by=motif_id ## for each tf keep only the peak with max log10 (to remove redundancy for heatmap)
      ]%>%unique()
)
## keep only relevant files
# fimo_cluster = fimo_cluster[c(1,4,5)]
fimo_cluster = Map(mutate,fimo_cluster,'file'=names(fimo_cluster))%>%rbindlist()

fimo_cluster = fimo_cluster[
    ,tf_cluster:=paste(motif_id,Name,sep=':')
    ][
      ,c('motif_id','Name'):=NULL
]
## make matrix
matrix_fimo =copy(fimo_cluster)
matrix_fimo = dcast(fimo_cluster, tf_cluster~file,value.var='log10_adjP')
matrix_fimo_rownames = matrix_fimo$tf_cluster
matrix_fimo = matrix_fimo[
  ,tf_cluster:=NULL
] 
matrix_fimo[is.na(matrix_fimo)]=0
matrix_fimo = as.matrix(matrix_fimo)
rownames(matrix_fimo) = matrix_fimo_rownames

## to keep only most enriched clusters remove rows with rowsum < x
matrix_fimo=matrix_fimo[rowSums(matrix_fimo)>=9, ]

## get random colors for motif cluster
library(RColorBrewer)
n <- unique(gsub(".*:","",rownames(matrix_fimo)))%>%length()
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

## need to find better solution than this to reduce the number fo colors/clusters to show
col_vector = c(col_vector,"#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
"#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
"#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
"#8A7C64", "#599861")%>%unique()

cluster_colors =sample(col_vector,n)

names(cluster_colors) = unique(gsub(".*:","",rownames(matrix_fimo)))


## plot heatmamp
fimo_heatmap = Heatmap(
      matrix_fimo, 
      border = T, 
      col= viridis(1000),
      row_dend_reorder = T,
      show_heatmap_legend = T,
      heatmap_legend_param = list(title = '-log10 fdr adj P'),
      show_row_names = T,
      row_labels = gsub(".*:","", rownames(matrix_fimo)),
      # row_order= order(gsub(".*:","", rownames(matrix_fimo))), 
      row_title =" ",
      show_column_names = T,
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 70,
      column_title =' ',
      right_annotation = HeatmapAnnotation(
                  which='row',
                  width = unit(1.5,'cm'),
                  cluster = anno_simple(
                      gsub(".*:","",rownames(matrix_fimo)),
                      col= cluster_colors
                      ),
                      show_annotation_name = F)
)


pdf(paste(plot_dir,'fimo_heatmap_pvals.pdf',sep=''),width= 10,height = 20)
fimo_heatmap
dev.off()



# ## script to calculate tf enrichment after fimo meme analysis
# ## enrichment for tfs across chrom state in da vs non da peaks (make heatmap)
# library(dplyr)
# library(data.table)
# library(magrittr)
# library(GenomicRanges)
# library(ggthemes);library(ggplot2);library(ggpubr)
# library(viridis);library(viridisLite)
# library(ggrepel)

# # hg38_fimo = snakemake@input[[]]
# # pantro5_fimo = snakemake@input[[]]
# # motif_cluster_file = snakemake@input[[]]
# # peaks = snakemake@input[[]]

# ## output
# plot_dir = '/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/plots/'

# merging_keys =  c("seqnames", "start",'end')

# read_fimo = function(x){
#     df= fread(x,sep='\t',header=T)
#     df=df[,motif_alt_id:=NULL]%>%
#     na.omit()%>%
#     setnames(old=c('start','stop','p-value','q-value'),new=c('motif_start','motif_end','p_value','q_value'))
#     df=df[q_value<=0.01]
#     df=df[,c(merging_keys,'species') := tstrsplit(sequence_name, "_", fixed=TRUE)]
#     df=df[,start:=as.numeric(start)][,end:=as.numeric(end)]
#     return(df)
#     }

# # hg38_fimo_results = read_fimo(hg38_fimo)
# # pantro5_fimo_results = read_fimo(pantro5_fimo)
# # all_tmm_peaks = fread(peaks, sep='\t',header=T,select=c('seqnames','start','end','width','FDR','direction'))
# # motif_cluster = fread(motif_cluster_file,sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

# hg38_fimo_results = read_fimo('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/fimo_meme/hg38/fimo.tsv')
# pantro5_fimo_results = read_fimo('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/fimo_meme/pantro5/fimo.tsv')

# all_fimo_results = rbind(hg38_fimo_results,pantro5_fimo_results)

# all_tmm_peaks = fread('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/DA/peaks/tmm_human_vs_chimp_csaw_DA_all_peaks.txt',
#     sep='\t',header=T)[,peak_type:=ifelse(FDR<=0.01,'da','non_da')]

# ## bimodal distribution fold change
# pdf(paste(plot_dir,'distribution_logFC.pdf',sep=''),width = 7,height = 7)
# ggplot(all_tmm_peaks,aes(x=logFC))+
#     geom_density()
# dev.off()

# ## add cluster info
# motif_cluster = fread('/data/projects/punim0586/dvespasiani/Files/Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

# all_fimo_results = all_fimo_results[motif_cluster,on='motif_id',nomatch=0]

# ## rank DA peaks by logFC
# da_peaks = copy(all_tmm_peaks)[peak_type=='da']
# da_peaks = da_peaks[
#     ,rank_fc:=logFC/max(abs(logFC))
#     ]

# da_peaks_tf = all_fimo_results[da_peaks,on=c(merging_keys),nomatch=0][
#     ,c(..merging_keys,'Cluster','Name','logFC','direction','FDR','peak_type','rank_fc')
# ]%>%unique()
# # [
# #     ,motif_id:=gsub("\\_.*","",motif_id)
# # ]%>%unique()



# ##--------------------------------------------------
# ## calculate enrichment per cluster 
# ## for tf in DA peaks vs tf in non-DA peaks
# ## for each species
# ##--------------------------------------------------
# ## first count number tf cluster in da/non-da peaks per species

# manipulate_file = function(a){

# create_matrix = function(x,organism){
#     df=copy(all_fimo_results)
#     df=df[species%in%organism][,c('seqnames','start','end','Cluster','Name','peak_type','species')] %>% unique()

#     df_da = copy(df)
#     df_da = df_da[peak_type=='da'][,tfcluster_da_peaks_per_species:=.N,by=.(Cluster,species)][,c(merging_keys):=NULL]%>%unique()
#     colnames(df_da)[5]=paste(colnames(df_da)[5],organism,sep='_')

#     df_nonda = copy(df)
#     df_nonda = df_nonda[peak_type=='non_da'][,tfcluster_non_da_peaks:=.N,by=.(Cluster)][,c(merging_keys):=NULL]%>%unique()
#     colnames(df_nonda)[5]=paste(colnames(df_nonda)[5],organism,sep='_')

#     df_merged = merge(df_da, df_nonda,by=c('Cluster','Name','species'), all=TRUE) 
#     df_merged$peak_type.y[is.na(df_merged$peak_type.y)] = 'non_da'
#     df_merged[is.na(df_merged)] = 1
#     df_merged = df_merged[,c('peak_type.x','peak_type.y','species'):=NULL]
#     return(df_merged)
# }

# chimp_tf_inpeaks = create_matrix(all_fimo_results,'chimp')
# human_tf_inpeaks = create_matrix(all_fimo_results,'human')


# matrix = chimp_tf_inpeaks[human_tf_inpeaks,on=c('Cluster','Name'),nomatch=0]

# fisher_pvalues=function(a){
#   pvals=copy(a)
#   pvals=pvals[,c(1:2):=NULL]
#   pvals=apply(pvals, 1, 
#               function(x) {
#                 tbl <- matrix(as.numeric(x[1:4]), ncol=2, byrow=T)
#                 fisher.test(tbl)$p.value
#               })
  
#   pvalues_table=data.table(pval=pvals)
#   pvalues_table=pvalues_table[
#     ,adj_p:=p.adjust(pval,method = 'fdr')][
#       ,log10_p_adjust:=-log10(adj_p)
#       ][
#         ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
#                                    ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
#                                           ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
#                                                  ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
#         ]
  
#   return(pvalues_table)
# }

# calculate_enrichment=function(x){

#   df = copy(x)

  
#   df_pvalues=copy(df)
#   df_pvalues=fisher_pvalues(df_pvalues)
  
#   df_final = cbind(df,df_pvalues)
#   df_final=df_final[
#     ,ratio_human_cluster_in_da_peaks := tfcluster_da_peaks_per_species_human/tfcluster_da_peaks_per_species_chimp][
#     ,mean_ratio:=mean(ratio_human_cluster_in_da_peaks)
#     ][
#       ,log2_fold_enrichment:=log2(ratio_human_cluster_in_da_peaks/mean_ratio)
#       ][
#         ,log10_numb_snps_tf:=log10(ratio_human_cluster_in_da_peaks) 
#         ][
#           ,c('Cluster','Name','log2_fold_enrichment','pval','adj_p',
#              'significant_score','log10_p_adjust','log10_numb_snps_tf')
#           ]%>% setorderv('log10_p_adjust',-1)
  
#   return(df_final)  
# }

# tf_enrichment = calculate_enrichment(matrix)

# ## plot
# tf_enrich_plot=function(x){
#   df=copy(x)[,'Log10 total number aSNPs per TF':=log10_numb_snps_tf]
#   gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per TF`),option="inferno",discrete = F)
#   text=ifelse(!df$significant_score%in%' ',df$Name,'')
  
#   ggplot(df,aes(x=log2_fold_enrichment,log10_p_adjust,label = text,col=log10_numb_snps_tf))+
#     geom_point(size=2)+
#     geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
#     geom_text_repel(size = 5,color='black',
#                     box.padding = unit(0.5, "lines"),
#                     point.padding = unit(0.5, "lines")
#     )+
#     gradient+
#     xlab('\n Log2 fold enrichment \n')+
#     ylab('-Log10 (adj-P)')+
#     xlim(-2,2)+
#        theme(
#       legend.position = "bottom",
#       legend.key = element_rect(fill = "white", colour = "black"),
#      axis.line = element_blank()
#      )
# }


# pdf(plot_tf_enrichment,width = 7,height = 7)
# tf_enrich_plot(tf_enrichment)
# dev.off()

