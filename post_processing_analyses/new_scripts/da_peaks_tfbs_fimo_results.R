library(dplyr)
library(data.table)
library(magrittr)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(biomaRt)
library(VennDiagram)
library(viridis)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

fimo_dir <- './output/fimo_meme'
outplot_dir <- create_dir(plot_dir,'tfbs')

outfile_dir <- create_dir(outdir,'files/fimo_ctcf')

## read da peaks
da_file <- paste(da_dir,'da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))

## read fimo results
fimo_results_files <- list.files(fimo_dir,recursive=T,full.names=T,pattern='fimo.tsv')
fimo_results <- fread(fimo_results_files,sep='\t',header=T)%>%setnames(old=c('p-value','q-value','sequence_name'),new=c('p_value','q_value','peakID'))
# adj_pvals <- copy(fimo_results)%>% dplyr::pull('p_value')%>%p.adjust('bonferroni')%>%as.data.table()%>%setnames('adj_pval')

fimo_results_adj <- copy(fimo_results)[
    q_value<=0.005
    ][
      ,log10_adjP:=-log10(q_value)
      ][
        ,c('motif_alt_id','strand'):=NULL
        ][
          da_results,on='peakID',nomatch=0
          ][
            ,gene:=gsub('\\_.*','',motif_id)
]
fimo_results_adj$seqlength <- nchar(fimo_results_adj$matched_sequence)
##----------------------------------------------------------------
## filter out motifs that are not expressed and add cluster info
##----------------------------------------------------------------
motif_cluster <- fread('/data/projects/punim0586/dvespasiani/Annotation_and_other_files/Motifs_clusters/motifs_clusters.txt',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

# expression_data <- fread(
#   "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
#   sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
#   )[
#     ,de_logFC:=-logFC
#     ][
#         ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
#         ][
#           ,c('P.Value','t','B','genes','logFC'):=NULL
# ]

# mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# ensembl_gene_symbols = getBM(
#     attributes = c('hgnc_symbol','ensembl_gene_id'), ## get all orth genes
#     mart = mart,
#     filters = 'ensembl_gene_id',
#     values = expression_data$EnsemblID
# )%>%as.data.table()%>%setnames(c('gene','EnsemblID'))

# expression_data <- expression_data[ensembl_gene_symbols,on='EnsemblID',nomatch=0][,ranked_expr:=abs(AveExpr)/max(abs(AveExpr))]

# fimo_results_adj_expr <- copy(fimo_results_adj)[expression_data,on='gene',nomatch=0]
# fimo_results_adj_expr_clust <- copy(fimo_results_adj_expr)[motif_cluster,on='motif_id',nomatch=0]

# ## count numb peaks with at least 1 motif and prop of tot numb peaks 
# seqid_allpeaks<-copy(da_results)[,c('peakID','DA')]%>%unique()
# seqid_allpeaks <- seqid_allpeaks[,tot:=.N,by=.(DA)][,c('DA','tot')]%>%unique()

# seqid_motif <-copy(fimo_results_adj)[,c('peakID','DA')]%>%unique()
# seqid_motif <- seqid_motif[,n:=.N,by=.(DA)][,c('DA','n')]%>%unique()

# prop_motif <- seqid_motif[seqid_allpeaks,on="DA"][,prop:=round((n/tot)*100,2)]
# prop_motif
# #        DA     n    tot  prop
# # 1:     da  2530   5262 48.08
# # 2: non_da 51298 110337 46.49

# seqid_motif_w_expr <-copy(fimo_results_adj_expr_clust)[,c('peakID','DA')]%>%unique()
# seqid_motif_w_expr <- seqid_motif_w_expr[,n:=.N,by=.(DA)][,c('DA','n')]%>%unique()

# prop_motif_w_expr <- seqid_motif_w_expr[seqid_allpeaks,on="DA"][,prop:=round((n/tot)*100,2)]
# prop_motif_w_expr
# #        DA     n    tot  prop
# # 1:     da  1339   5262 25.45
# # 2: non_da 30970 110337 28.07

# ##--------
# ## QCs
# ##--------
# ## number of motifs per peak (and avg)
# numb_motif_peak <- copy(fimo_results_adj_expr_clust)[,c('peakID','motif_id','DA')]%>%unique()
# numb_motif_peak<-numb_motif_peak[
#     ,numb_motifs_per_peak:=.N,by=.(peakID,DA)
#     ][
#       ,avg_motif_per_peak:=round(mean(numb_motifs_per_peak),2),by=.(DA)
# ]
# ## avg motifs per peak: 4.34 da; 4.13 non-da

# ## distribution number motifs per peak
# pdf(paste(outplot_dir,'qc_number_motifs_per_peak.pdf',sep=''),width = 10,height = 7)
# ggplot(numb_motif_peak,aes(x=numb_motifs_per_peak,fill=DA))+
# geom_bar()+facet_wrap(DA~.,scale='free')
# dev.off()

# ## check correlation with peak sizes
# corr_numb_motif_peak_size <- numb_motif_peak[da_results,on=c('peakID','DA'),nomatch=0][,width:=end-start]

# cor.test(corr_numb_motif_peak_size$width,corr_numb_motif_peak_size$numb_motifs_per_peak)
# # data:  corr_numb_motif_peak_size$width and corr_numb_motif_peak_size$numb_motifs_per_peak
# # t = 102.13, df = 84136, p-value < 2.2e-16
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.3260805 0.3381040
# # sample estimates:
# #       cor 
# # 0.3321057 

##--------------------
## permutation
##--------------------
## read all hocomoco motifs used in FIMO
hocomoco_motifs <- universalmotif::read_meme('/data/projects/punim0586/dvespasiani/Annotation_and_other_files/MEME_files/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')
hocomoco_motifs <- lapply(hocomoco_motifs,function(x)as.data.table(x@name))%>%rbindlist()%>%setnames('motif_id')

myclusters <- copy(fimo_results_adj[,c('peakID','motif_id')]%>%unique()%>%dplyr::pull('motif_id'))

observed_value = as.data.table(table(myclusters))%>%setnames(old=c(1:2),new=c('motif_id','score'))%>%split(by='motif_id')

list_df_permuted_scores = list()
for(i in 1:10000){
      permuted_scores=table(sample(hocomoco_motifs$motif_id, length(unique(fimo_results_adj$peakID)), replace = T))
      list_df_permuted_scores[[i]] = as.data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('motif_id','permuted_score'))
}
list_df_permuted_scores = rbindlist(list_df_permuted_scores)[motif_id%in%names(observed_value)]%>%split(by='motif_id')

genwide_motif_enrichment = purrr::map2(list_df_permuted_scores,observed_value,function(p,o){
  zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
  pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
  stat =  data.table(zscore = zscore,pval=pvalue)
  return(stat)
  }
)
genwide_motif_enrichment <- Map(mutate,genwide_motif_enrichment,motif_id=names(genwide_motif_enrichment))%>%rbindlist()%>%setorderv('zscore',1)

genwide_motif_enrichment <- genwide_motif_enrichment[
  ,logzscore:=ifelse(zscore<0,-log(abs(zscore)),log(abs(zscore)))
  ][
    ,gene:=gsub('\\_.*','',motif_id)
    ][
      expression_data, on='gene',nomatch=0
      ][
        motif_cluster,on='motif_id',nomatch=0
        ][
          ,c('Cluster','Database','Seed_motif','motif_id'):=NULL
          ][
            ,tf_symbol:=paste(gene,Name,sep=':')
]


## make heatmaps
library(ComplexHeatmap)

pvals <- copy(genwide_motif_enrichment)[,c('DE','tf_symbol')][,p_sign:=ifelse(DE=='de','*',NA)]%>%dplyr::pull('p_sign')
names(pvals) = genwide_motif_enrichment$tf_symbol

matzscore_rownames <- copy(genwide_motif_enrichment$tf_symbol)
mat_zscore <- copy(genwide_motif_enrichment)[,logzscore]%>%as.matrix()
rownames(mat_zscore) = matzscore_rownames

colfun_zscore = circlize::colorRamp2(c(min(mat_zscore), 0, max(mat_zscore)), c("#0582CA", "white", "#d62828"))

## cluster colors
qual_col_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_palette$maxcolors, rownames(qual_col_palette)))

cluster_colors = sample(col_vector,length(unique(gsub(".*:","",rownames(mat_zscore)))))
names(cluster_colors) = unique(gsub(".*:","",rownames(mat_zscore)))

heatmap_zscore <- Heatmap(
  mat_zscore,
  name='log zscore \n motif enrichment',
  col=colfun_zscore,
  right_annotation = HeatmapAnnotation(
                  which='row',
                  de_pvalue = anno_simple(pvals,width = unit(1,'cm'), pch = pvals),
                  cluster = anno_simple(
                    gsub(".*:","",rownames(mat_zscore)),
                    width = unit(1,'cm'),
                    col=cluster_colors
                    )
                  )
)

mat_expr <- copy(genwide_motif_enrichment)[,ranked_expr]%>%as.matrix()%>%abs()
rownames(mat_expr) = matzscore_rownames

colfun_expr = circlize::colorRamp2(c(min(mat_expr),median(mat_expr), max(mat_expr)), c("#FFFF3F",'#BFD200',"#2B9348"))

heatmap_expr <- Heatmap(
  mat_expr,
  name='ranked avg gene expression \n (log2 RPKM)',
  col=colfun_expr
)

pdf(paste(outplot_dir,'heatmap_expressed_enriched_motifs.pdf',sep=''),width=10,height = 15)
heatmap_expr + heatmap_zscore
dev.off()

##--------------------------------
## overlap da vs non-da
##--------------------------------

motifs_ids <- copy(fimo_results_adj_expr_clust)[,c('Name','DA')]%>%unique()%>%split(by='DA')%>%lapply(function(x)x$Name)

venn.diagram(
    x = motifs_ids,
    category.names = c('non_da','da'),
    filename = paste(outplot_dir,'set_of_motif_clusters_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col=c(da_palette[[2]],da_palette[[1]]),
    fill = c(alpha(da_palette[[2]],0.3), alpha(da_palette[[1]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-30, 30),
    cat.dist = c(0.1, 0.1),
    cat.fontfamily = "sans",
    cat.col = c(da_palette[[2]],da_palette[[1]])
)

## OR
fimo_results_adj_expr_clust_counts <- copy(fimo_results_adj_expr_clust)[
    ,numb_motifs_cluster:=.N,by=.(Name,DA)
      ][
        ,numb_motifs:=.N,by=.(DA)
        ][
          ,numb_motifs_noncluster:=numb_motifs-numb_motifs_cluster
]

counts_motifs_clusters <- copy(fimo_results_adj_expr_clust_counts)%>%split(by='DA')%>%lapply(
  function(x){x<-x[,c('Name','numb_motifs_cluster','numb_motifs_noncluster')]%>%unique()}
)

da_vs_nonda_or <- calculate_or(counts_motifs_clusters[[2]],counts_motifs_clusters[[1]],'Name')%>%setnames(old='elements',new='Name')%>%adjust_pvalues()
da_vs_nonda_or <- da_vs_nonda_or[
  counts_motifs_clusters[[2]],on='Name',nomatch=0
  ][
    ,numb_motifs_noncluster:=NULL
    ][
      ,log10_numb_snps_tf:=log10(numb_motifs_cluster)
      ][
        ,log10_padj:= -log10(p.adj)
        ][
          ,log10_padj:=ifelse(log10_padj=='Inf',0,log10_padj)
          ][
            ,replace_p:=ifelse(p.signif!=' ' & log10_padj==0,'y','n')
]

## plot
tf_enrich_plot=function(x){
  df <- copy(x)[,'Log10 total number motifs per cluster':=log10_numb_snps_tf]
  gradient <- scale_colour_viridis(aes(`Log10 total number motifs per cluster`),option="inferno",discrete = F)
  text <- ifelse(!df$p.signif%in%' ',df$Name,'')
  
  ggplot(df,aes(x=log2(odds_ratio),log10_padj,label = text,col=log10_numb_snps_tf))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
    geom_text_repel(size = 5,color='black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines"),
                    max.overlaps = 50
    )+
    gradient+
    geom_hline(yintercept=1.3,linetype='dashed',size=0.2)+
    xlab('\n Log2 DA vs non-DA OR \n')+
    ylab('-Log10 (adj-P)')+
    xlim(-max(df$odds_ratio+0.5),max(df$odds_ratio+0.5))+
    theme_classic()+ 
    theme(
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", colour = "black")
        )
}

pdf(paste(outplot_dir,'da_vs_nonda_or_fimo.pdf',sep=''),width = 7,height = 7)
tf_enrich_plot(da_vs_nonda_or)
dev.off()














# ## keep only motif with highest score per peak
# # fimo_out_topmotif <- copy(fimo_out)[,.SD[which.max(score)], by=.(peakID)]

# ## Overlap motif clusters between non-da & da peaks 
# library(VennDiagram)

# motifs_ids <- copy(fimo_out)[,c('motif_id','DA')]%>%unique()%>%split(by='DA')%>%lapply(function(x)x$motif_id)

# venn.diagram(
#     x = motifs_ids,
#     category.names = c('non_da','da'),
#     filename = paste(outplot_dir,'set_of_motif_clusters_venn.png',sep=''),
#     output = TRUE ,
#     imagetype="png" ,
#     height = 700 , 
#     width = 700 , 
#     resolution = 400,
#     lwd = 1,
#     col=c(da_palette[[2]],da_palette[[1]]),
#     fill = c(alpha(da_palette[[2]],0.3), alpha(da_palette[[1]],0.3)),
#     cex = 0.5,
#     fontfamily = "sans",
#     cat.cex = 0.3,
#     cat.default.pos = "outer",
#     cat.pos = c(-30, 30),
#     cat.dist = c(0.1, 0.1),
#     cat.fontfamily = "sans",
#     cat.col = c(da_palette[[2]],da_palette[[1]])
# )

# # ## proceed with odds ratio of tfs per cluster
# # fimo_out_cluster_top_motif <- copy(fimo_out_cluster)[
# #       ,numb_motifs_cluster:=.N,by=.(Name,DA)
# #       ][
# #         ,numb_motifs:=.N,by=.(DA)
# #         ][
# #           ,numb_motifs_noncluster:=numb_motifs-numb_motifs_cluster
# # ]


# # counts_motifs_clusters <- copy(fimo_out_cluster_top_motif)%>%split(by='DA')%>%lapply(
# #   function(x){x<-x[,c('Name','numb_motifs_cluster','numb_motifs_noncluster')]%>%unique()}
# # )

# # da_vs_nonda_or <- calculate_or(counts_motifs_clusters[[1]],counts_motifs_clusters[[2]],'Name')%>%setnames(old='elements',new='Name')%>%adjust_pvalues()
# # da_vs_nonda_or <- da_vs_nonda_or[
# #   counts_motifs_clusters[[2]],on='Name',nomatch=0
# #   ][
# #     ,numb_motifs_noncluster:=NULL
# #     ][
# #       ,log10_numb_snps_tf:=log10(numb_motifs_cluster)
# #       ][
# #         ,log10_padj:= -log10(p.adj)
# #         ][
# #           ,log10_padj:=ifelse(log10_padj=='Inf',0,log10_padj)
# #           ][
# #             ,replace_p:=ifelse(p.signif!=' ' & log10_padj==0,'y','n')
# # ]

# # ## plot
# # tf_enrich_plot=function(x){
# #   df <- copy(x)[,'Log10 total number motifs per cluster':=log10_numb_snps_tf]
# #   gradient <- scale_colour_viridis(aes(`Log10 total number motifs per cluster`),option="inferno",discrete = F)
# #   text <- ifelse(!df$p.signif%in%' ',df$Name,'')
  
# #   ggplot(df,aes(x=log2(odds_ratio),log10_padj,label = text,col=log10_numb_snps_tf))+
# #     geom_point(size=2)+
# #     geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
# #     geom_text_repel(size = 5,color='black',
# #                     box.padding = unit(0.5, "lines"),
# #                     point.padding = unit(0.5, "lines"),
# #                     max.overlaps = 50
# #     )+
# #     gradient+
# #     geom_hline(yintercept=1.3,linetype='dashed',size=0.5)+
# #     xlab('\n Log2 non da vs da OR \n')+
# #     ylab('-Log10 (adj-P)')+
# #     xlim(-8,8)+
# #     theme_bw()+ 
# #     theme(
# #         legend.position = "bottom",
# #         legend.key = element_rect(fill = "white", colour = "black"),
# #         axis.line = element_blank())
# # }

# # pdf(paste(outplot_dir,'da_vs_nonda_or_fimo.pdf',sep=''),width = 10,height = 10)
# # tf_enrich_plot(da_vs_nonda_or)
# # dev.off()
