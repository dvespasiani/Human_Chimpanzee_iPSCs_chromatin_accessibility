## script used to analyse homer tfbs known motifs results and get enrichment
library(dplyr)
library(data.table)
library(magrittr)
library(ComplexHeatmap)
library(circlize)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(viridis)
library(viridisLite)
library(RColorBrewer)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

da_peaks_dir = './output/DA/peaks/'
homer_dir = './output/homer/homer_output/'
plot_dir = './output/plots/tfbs/'

##-------------
## read motifs
##-------------
read_known_motifs = function(dir){
  known_motifs = fread(paste(homer_dir,dir,'knownResults.txt',sep=''),sep='\t',header=T)
  known_motifs = known_motifs[,motif_name:=gsub("\\/.*","",`Motif Name`)]
  families =stringr::str_extract(string = known_motifs$motif_name,pattern = "(?<=\\().*(?=\\))")
  
  known_motifs= known_motifs[
    ,motif_family:=families
    ][
      ,motif_name:=gsub("\\(.*","",motif_name)
      ][
        ,prop_motif_in_target:=as.numeric(gsub("\\%.*","",`% of Target Sequences with Motif`))
        ][
          ,c('Motif Name','Consensus','% of Target Sequences with Motif'):=NULL
          ]
  return(known_motifs)
}

chimp_known = read_known_motifs('chimp_specific_pantro5/')
common_known = read_known_motifs('common_regions_hg38/')
human_known = read_known_motifs('human_specific_hg38/')

known_motif = list(chimp_known,common_known,human_known)
known_motif = Map(mutate,known_motif,species=species_names)

number_initial_motifs = copy(known_motif)%>%lapply(function(x)as.data.table(x)[,'motif_name']%>%unique()%>%nrow())
##----------------
## filter motifs
##----------------
## FDR adjust pvalues
## remove motifs with fdr adj.pval > 0.01 
## remove motifs present in less than 5% of target sequences

filter_motifs= function(motif){
  filtered=copy(motif)
  # %>%as.data.table()
  # filtered$log_p_value[(filtered$log_p_value == Inf)] = max(filtered$log_p_value[is.finite(filtered$log_p_value)]+200)## this to counteact Inf
  filtered=filtered[
    ,log_pval:=ifelse(`Log P-value`< -729,-745,`Log P-value`) ## done this to control R from converting <<< 0 values to 0
    ][
      ,log10_fdr_adjP:=-log10(p.adjust(exp(log_pval),'fdr'))
      ][
        log10_fdr_adjP > 2
        ][
          prop_motif_in_target>5
          ][
            !motif_name %like% 'Unknown'
            ][
              ,motif_name:=ifelse(motif_name %like% 'BORIS','CTCFL',motif_name)
              ]

  # filtered=filtered[
  #   p.adj< 0.01
  #   ][
  #     ,tgt_pct:=round(tgt_pct,2)
  #     ][
  #       tgt_pct>0.05
  #       ][
  #         !motif_name %like% 'Unknown'
  #         ][
  #           ,motif_name:=ifelse(motif_name %like% 'BORIS','CTCFL',motif_name)
  #           ]
  # filtered=tibble(filtered)
  # names(filtered$motif_pwm)=filtered$motif_name
  return(filtered)

}

filtered_known_motif = lapply(known_motif,function(x)filter_motifs(x))

number_filtered_motifs = copy(filtered_known_motif)%>%lapply(function(x)as.data.table(x)[,'motif_name']%>%unique()%>%nrow())

##--------------
## plot results
##--------------
## important: GFY is a general factor Y, which refers to a motif for which HOMER doesnt know which if the factor that binds it (see http://homer.ucsd.edu/homer/motif/motifDatabase.html)
## however, these motifs are real and some TFs have been described binding them (see http://genesdev.cshlp.org/content/24/14/1479.long)
## what I'll do is to report the family for this as 'unknown' and in the text I will mention this and cite the above paper.
motifs_table = copy(filtered_known_motif)%>%lapply(function(x)x=x[,c('motif_name','motif_family','species','log10_fdr_adjP')])%>%rbindlist()
motifs_table = motifs_table[,motif_family:=ifelse(motif_family=='?','Unknown',motif_family)][,motif_name_family:=paste(motif_name,motif_family,sep=':')][,c('motif_name','motif_family'):=NULL]

## create matrix for heatmap
matrix = copy(motifs_table)
matrix = dcast(matrix, motif_name_family~species,value.var='log10_fdr_adjP')
matrix[is.na(matrix)]=0
matrix_rownames = matrix$motif_name_family
matrix = matrix[,motif_name_family:=NULL]%>%as.matrix()
rownames(matrix) = matrix_rownames

filtered_matrix = copy(matrix)
filtered_matrix=filtered_matrix[rowSums(filtered_matrix)>17, ]
filtered_matrix[!rownames(filtered_matrix) %like% 'Zic2', ]

## get random colors for motif families 
n_tf_families = length(unique(gsub(".*:","",rownames(filtered_matrix))))

qual_col_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_palette$maxcolors, rownames(qual_col_palette)))

tf_family_colors = sample(col_vector,n_tf_families)
names(tf_family_colors) = unique(gsub(".*:","",rownames(filtered_matrix)))


## plot heatmamp
homer_heatmap = function(x){
  Heatmap(
      x, 
      border = T, 
      col= viridis(1000),
      row_dend_reorder = T,
      # row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
      show_heatmap_legend = T,
      heatmap_legend_param = list(title = '-log10 fdr adj P'),
      show_row_names = T,
      row_labels = gsub("\\:.*", "", rownames(x)),
      row_title =" ",
      show_column_names = T,
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 70,
      column_title =' ',
      right_annotation = HeatmapAnnotation(
                  which='row',
                  width = unit(1.5,'cm'),
                  cluster = anno_simple(
                      gsub(".*:","",rownames(x)),
                      col= tf_family_colors
                      ),
                      show_annotation_name = F)
                      )
}

pdf(paste(plot_dir,'homer_known_heatmap.pdf',sep=''),width= 10,height = 20)
homer_heatmap(filtered_matrix)
dev.off()


##----------------------------
## create motif logo object
##----------------------------
## for this read again HOMER results but using marge 
## this allows to retrieve the PWMs for plotting
library(marge)
library(motifStack)

chimp_known_pwms = read_known_results(paste(homer_dir,'chimp_specific_pantro5/',sep=''))
common_known_pwms = read_known_results(paste(homer_dir,'common_regions_hg38/',sep=''))
human_known_pwms = read_known_results(paste(homer_dir,'human_specific_hg38/',sep=''))

my_motifs = list(chimp_known_pwms,common_known_pwms,human_known_pwms)%>%lapply(function(y)y=dplyr::select(y,c('motif_name','motif_pwm'))%>%as.data.table())
my_motifs = lapply(my_motifs,function(x)x=x[,motif_name:=ifelse(motif_name %like% 'BORIS','CTCFL',motif_name)][motif_name %in% gsub('\\:.*','',rownames(filtered_matrix))])%>%rbindlist()%>%tibble()%>%unique()
names(my_motifs$motif_pwm)=my_motifs$motif_name

motif_names=as.list(my_motifs$motif_name)
motifs_pwm=copy(my_motifs$motif_pwm)

get_motif_logo=function(pwm,names){
  motifs=copy(pwm)%>%as.matrix()%>%t()
  colnames(motifs)=paste('V',1:dim(motifs)[2],sep='')

  motif_object= new("pcm", mat=motifs, name=names)
  return(motif_object)
}

my_motif_object= purrr::map2(motifs_pwm,motif_names,function(x,y)get_motif_logo(x,y))

pdf(paste(plot_dir,'logo_homer_top_known_motifs.pdf',sep=''),width= 15,height = 5)
lapply(my_motif_object,function(x)plot(x))
dev.off()

## now make a barplot with proportion motif families
motif_families=data.table(motif_name_family=rownames(filtered_matrix))
motif_families = motif_families[
  ,c("motif_name", "motif_family") := tstrsplit(motif_name_family, ":", fixed=TRUE)
  ][
    ,'motif_name_family':=NULL
    ][
      ,numb_motifs_by_family:=.N,by=.(motif_family)
      ][
        ,numb_motifs:=.N
        ][
          ,prop_family:=numb_motifs_by_family/numb_motifs
          ][
            ,mock_col:='tf_family'
          ]
prop_motif_families = copy(motif_families)[,c('prop_family','motif_family','mock_col')]%>%unique()

pdf(paste(plot_dir,'prop_families_homer_top_known_motifs.pdf',sep=''),width= 5,height = 10)
ggplot(prop_motif_families,aes(x=mock_col,y=prop_family,fill=motif_family))+
geom_bar(position="fill", stat="identity")+
scale_fill_manual(values = tf_family_colors)+ theme_void()
dev.off()


## DEPRECATED
# read_tfbs_results = function(dir,pattern){
#   df = list.files(dir,recursive=T,full.names=T,pattern=pattern)%>%
#     lapply(function(x)x=fread(x,sep='\t',header=T))
#   names(df) =  gsub("\\/.*","",list.files(dir,recursive=T,full.names=F,pattern=pattern))
#   return(df)
# }

# ## read both de novo and known tfbs
# homer_denovo_results = read_tfbs_results(homer_dir,'homer_denovo_motifs.txt')

# ## numb motifs
# lapply(homer_denovo_results,function(x)x=x[!file_name %like% 'similar|RV'][,c('file_name')]%>%unique()%>%nrow())

# ##------------------------------------------
# ## start manipulating the de novo results
# ##------------------------------------------
# # homer_denovo_results = 
# # lapply(
# #   homer_denovo_results,function(x)x=x[
# #     ,motif_occurrences:=gsub("\\,.*","",motif_occurrences)
# #     ][
# #       ,motif_occurrences:=sub("\\).*", "", sub(".*\\(", "", motif_occurrences)) 
# #       ][
# #         ,motif_occurrences := as.numeric(gsub("\\%.*","",motif_occurrences))
# #         ][
# #           motif_occurrences>5 ## remove tfs present in < 5% of test set as suggested by homer
# #           ][
# #             ,log_pval:=ifelse(log_pval< -729,-745,log_pval) ## done this to control R from converting <<< 0 values to 0
# #             ][
# #               ,log10_fdr_adjP:=-log10(p.adjust(exp(log_pval),'fdr'))
# #               ]
# # )
# homer_denovo_results = 
# lapply(
#   homer_denovo_results,function(x)x=x[
#     ,motif_occurrences:=gsub("\\,.*","",motif_occurrences)
#     ][
#       ,motif_occurrences:=sub("\\).*", "", sub(".*\\(", "", motif_occurrences)) 
#       ][
#         ,motif_occurrences := as.numeric(gsub("\\%.*","",motif_occurrences))
#         ][
#           motif_occurrences>5 ## remove tfs present in < 5% of test set as suggested by homer
#           ][
#             ,log_pval:=ifelse(log_pval< -729,-745,log_pval) ## done this to control R from converting <<< 0 values to 0
#             ][
#               ,log10_fdr_adjP:=-log10(p.adjust(exp(log_pval),'fdr'))
#               ]
# )
# homer_denovo_results  = Map(mutate,homer_denovo_results,"file"=(names(homer_denovo_results)))

# homer_denovo_results = 
# lapply(homer_denovo_results,function(x)x=x[
#   !file_name %like% 'similar|RV'
#   ][motif_name %like% 'BestGuess'][,file_name:=NULL]%>%unique()
# )
# homer_denovo_results = lapply(
#   homer_denovo_results,function(x)x=x[
#     ,motif_id:=toupper(gsub(".*:","",gsub("\\/.*","",motif_name)))
#     ][
#       ,motif_id:=gsub("\\s*\\([^\\)]+\\)","",motif_id) ## remove stuff within parentheses
#       ][
#         ,motif_id:=ifelse(motif_id=='BORIS','CTCFL',motif_id) 
#         ]
# )%>%rbindlist()
# #    [
# #     ,motif_id:=toupper(gsub(".*:","",gsub("\\/.*","",motif_name)))
# #     ][
# #       ,motif_id:=gsub("\\s*\\([^\\)]+\\)","",motif_id) ## remove stuff within parentheses
# #       ][
# #         ,motif_id:=ifelse(motif_id=='BORIS','CTCFL',motif_id) 
# #         ][
# #           !motif_id %like% 'UNKNOWN' ## remove these
# #           ]
# # )

# number_final_denovo_motifs = copy(homer_denovo_results)[,c('file','motif_id')]%>%unique()
# number_final_denovo_motifs = number_final_denovo_motifs[,number_final_denovo_motifs:=.N,by=.(file)][,c('file','number_final_denovo_motifs')]%>%unique()

# ##------------------------------------------
# ## now do the same with the known results
# ##------------------------------------------
# # homer_known_results = read_tfbs_results(homer_dir,'knownResults.txt')
# # homer_known_results = lapply(homer_known_results,function(x)x=setnames(x,old=c('Motif Name','P-value','% of Target Sequences with Motif','Log P-value'),
# # new=c('motif_name','pvalue','prop_target_seq_with_motif','log_pval'))%>%setorderv('log_pval',1))

# # ## numb motifs
# # lapply(homer_known_results,function(x)x=x[,motif_id:=toupper(gsub("\\(.*","",motif_name))][,c('motif_id')]%>%unique()%>%nrow())

# # homer_known_results = lapply(homer_known_results,function(x)x=x[
# #     ,prop_target_seq_with_motif:=as.numeric(gsub("\\%.*","",prop_target_seq_with_motif))
# #     ][
# #       prop_target_seq_with_motif>5
# #       ][
# #         ,log_pval:=ifelse(log_pval< -729,-745,log_pval) ## done this to control R from converting <<< 0 values to 0
# #         ][
# #           ,log10_fdr_adjP:=-log10(p.adjust(exp(log_pval),'fdr'))
# #           ][
# #             ,motif_id:=toupper(gsub("\\(.*","",motif_name))
# #             ][
# #                 ,motif_id:=ifelse(motif_id =='BORIS','CTCFL',motif_id) 
# #                 ][
# #                   !motif_id%like%'UNKNOWN' ## remove these
# #                 ][
# #                   ,c('motif_id','log10_fdr_adjP')
# #                 ]%>%unique()
# # )

# # homer_known_results = Map(mutate,homer_known_results,'file'=names(homer_known_results))%>%rbindlist()

# # number_final_known_motifs = copy(homer_known_results)[,c('file','motif_id')]%>%unique()
# # number_final_known_motifs =number_final_known_motifs[,number_final_known_motifs:=.N,by=.(file)][,c('file','number_final_known_motifs')]%>%unique()

# ## add motif cluster info
# motif_cluster = fread('/data/projects/punim0586/dvespasiani/Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')
# motif_cluster = motif_cluster[,motif_id:=toupper(motif_id)][,motif_id:=gsub("\\_.*","",motif_id)]

# ## keep only those TF that can be assigned to a cluster
# homer_denovo_results = merge(homer_denovo_results,motif_cluster, by = "motif_id")

# ## same for known tfbs
# homer_known_results = merge(homer_known_results,motif_cluster, by = "motif_id")

# ## make matrix
# make_matrix = function(x){
#   matrix=copy(x)[,c('log10_fdr_adjP','file','motif_id','Name')]%>%unique()
#   all_tfs= copy(matrix)[,c('motif_id','Name')]%>%unique()
  
#   matrix = split(matrix,by=c('file'))%>%
#   lapply(function(x)full_join(x,all_tfs,by=c('motif_id','Name'))[
#     ,log10_fdr_adjP:=ifelse(is.na(log10_fdr_adjP)==T,0,log10_fdr_adjP)
#     ][
#       ,file:=NULL
#       ][
#         , .SD[which.max(log10_fdr_adjP)], by=motif_id
#         ]
#   )
#   print(lapply(matrix,function(x)x[log10_fdr_adjP>0]%>%nrow())) ## prints number of tfs retained after this last step
#   matrix  = Map(mutate,matrix,"file"=(names(matrix)))%>%rbindlist()
#   matrix = matrix[,tf_cluster:=paste(motif_id,Name,sep=':')][,motif_id:=NULL]

#   matrix = dcast(matrix, tf_cluster~file,value.var='log10_fdr_adjP')
#   matrix_rownames = matrix$tf_cluster
  
#   matrix = matrix[,tf_cluster:=NULL]%>%as.matrix()
#   rownames(matrix) = matrix_rownames
#   return(matrix)

# }
# matrix_homer_denovo_results = make_matrix(homer_denovo_results)
# matrix_homer_known_results = make_matrix(homer_known_results)

# ## do this to simplify matrix and keep only most significant ones to plot
# matrix_homer_denovo_results = matrix_homer_denovo_results[rowSums(matrix_homer_denovo_results)>0, ]
# matrix_homer_known_results = matrix_homer_known_results[rowSums(matrix_homer_known_results)>30, ]

# ## get random colors for motif cluster
# library(RColorBrewer)
# n <- unique(
#   c(gsub(".*:","",rownames(matrix_homer_known_results)),
#    gsub(".*:","",rownames(matrix_homer_denovo_results)) 
# ))%>%length()

# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# ## need to find better solution than this to reduce the number fo colors/clusters to show
# col_vector = c(col_vector,"#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
# "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
# "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
# "#8A7C64", "#599861")%>%unique()

# cluster_colors = sample(col_vector,n)

# names(cluster_colors) = unique( 
#   c(gsub(".*:","",rownames(matrix_homer_known_results)), 
#   gsub(".*:","",rownames(matrix_homer_denovo_results))
#   )
# )

# ## plot heatmamp
# homer_heatmap = function(x){
#   Heatmap(
#       x, 
#       border = T, 
#       col= viridis(1000),
#       row_dend_reorder = T,
#       # row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
#       show_heatmap_legend = T,
#       heatmap_legend_param = list(title = '-log10 fdr adj P'),
#       show_row_names = T,
#       row_labels = gsub("\\:.*", "", rownames(x)),
#       row_title =" ",
#       show_column_names = T,
#       column_names_gp = gpar(fontsize = 10),
#       column_names_rot = 70,
#       column_title =' ',
#       right_annotation = HeatmapAnnotation(
#                   which='row',
#                   width = unit(1.5,'cm'),
#                   cluster = anno_simple(
#                       gsub(".*:","",rownames(x)),
#                       col= cluster_colors
#                       ),
#                       show_annotation_name = F)
#                       )
# }

# pdf(paste(plot_dir,'homer_denovo_heatmap.pdf',sep=''),width= 10,height = 20)
# homer_heatmap(matrix_homer_denovo_results)
# dev.off()

# pdf(paste(plot_dir,'homer_known_heatmap.pdf',sep=''),width= 10,height = 20)
# homer_heatmap(matrix_homer_known_results)
# dev.off()