## script used to look for overlap between DA peaks and DE genes
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(VennDiagram)

options(width=150,scipen=99)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'da_de_overlap')
outfile_dir <- create_dir(outdir,'files/GO_enrich')

target_genes_dir <- './output/files/GO_enrich'
tmp_files_dir <- './output/temp_files/'

de_palette <- c('#DDA15E','#606C38')
names(de_palette) = c('de','non_de')

# ## get DA peaks
# da_results <- read_da_results('new_da_results.txt')
# da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

# ## use chipseeker to annotate peaks within 1Mb of distance to their putative target gene(s)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# tss_region = 1000000

# peakAnno <- annotatePeak(
#   makeGRangesFromDataFrame(da_results,keep.extra.columns=T), 
#   tssRegion=c(-tss_region, tss_region),
#   TxDb=txdb, 
#   annoDb="org.Hs.eg.db"
# )%>%as.data.table()%>%setnames(old='ENSEMBL',new='EnsemblID')

# ## raw numbers and proportion of peaks with target genes
# numb_peaks_associated <- copy(peakAnno)[,c('peakID','DA','EnsemblID')]%>%unique()%>%na.omit()
# numb_peaks_associated <- numb_peaks_associated[,EnsemblID:=NULL][,n:=.N,by=.(DA)][,c('DA','n')]%>%unique()
# tot_numb_peaks <- copy(da_results)[,c('peakID','DA')][,tot_n:=.N,by=.(DA)][,c('DA','tot_n')]%>%unique()
# prop <- numb_peaks_associated[tot_numb_peaks,on='DA',nomatch=0][,prop:=n/tot_n]
# #        DA     n tot_n      prop
# # 1: non_da 52851 55731 0.9483232
# # 2:     da 16592 17701 0.9373482

# ## look @ the overlap between genes targeted by DA and non-DA peaks 
# genes <- copy(peakAnno)%>%split(by='DA')%>%lapply(function(x)x=x[,EnsemblID]%>%unique()%>%na.omit())
# venn.diagram(
#     x = genes,
#     category.names = c('da','non_da'),
#     filename = paste(outplot_dir,'overlap_genes_da_nonda.png',sep=''),
#     output = TRUE ,
#     imagetype="png" ,
#     height = 700 , 
#     width = 700 , 
#     resolution = 400,
#     lwd = 1,
#     col=da_palette,
#     fill = c(alpha(da_palette[[1]],0.3), alpha(da_palette[[2]],0.3)),
#     cex = 0.5,
#     fontfamily = "sans",
#     cat.cex = 0.3,
#     cat.default.pos = "outer",
#     cat.pos = c(-27, 27),
#     cat.dist = c(0.055, 0.055),
#     cat.fontfamily = "sans",
#     cat.col = da_palette
# )

# ## plot the distribution of the numbers of peaks per gene
# pdf(paste(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width=10,height = 7)
# df <- copy(peakAnno)[,c('EnsemblID','DA','peakID')]%>%na.omit()%>%unique()
# df <- df[,numb_peaks_per_gene:=.N,by=.(DA,EnsemblID)]
# ggplot(df,aes(x=numb_peaks_per_gene,fill=DA))+
# geom_bar(position=position_dodge())+
# scale_fill_manual(values=da_palette)+
# facet_wrap(DA~.,ncol=2)+
# theme_classic()+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black")
# )
# dev.off()

# ## mean number of peaks associated per gene
# # mean(df[DA=='da']$numb_peaks_per_gene)
# # [1] 5.25
# # mean(df[DA!='da']$numb_peaks_per_gene)
# # [1] 8.15

# ## assign closest peak to gene
# one_peak_per_gene <- copy(peakAnno)[,distTSS:=start-geneStart][,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]

# pdf(paste(outplot_dir,'dist_peak_tss.pdf',sep=''),width=10,height = 7)
# df <- copy(one_peak_per_gene)[,c('distTSS','DA','peakID')]%>%na.omit()%>%unique()
# df <- df[,log10dist:=log10(abs(distTSS)+0.1)]
# ggplot(df,aes(x=log10dist,fill=DA))+
# geom_density(alpha=0.5)+
# scale_fill_manual(values=da_palette)+
# theme_classic()+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black")
# )
# dev.off()

# ## add DE information
# ## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
# de_genes <- fread(
#   "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
#   sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
#   )[
#     ,de_logFC:=-logFC
#     ][
#         ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
#         ][
#           ,c('P.Value','t','B','genes','logFC'):=NULL
# ]

# da_peaks_target_genes <- copy(peakAnno)[
#   ,c(..range_keys,'DA','peakID','FDR','peak_species','logFC','geneStart','ENSEMBL')
#   ][
#     ,distTSS:=start-geneStart
#     ][
#       ,geneStart:=NULL
#       ][
#         de_genes,on='EnsemblID',nomatch=0
# ]

# ## raw numbers and proportion of peaks with target orthologous genes
# get_prop <- function(initial,filtered,metric){
#   i <- copy(initial)%>%dplyr::pull(metric)%>%unique()%>%na.omit()%>%length()
#   f <- copy(filtered)%>%dplyr::pull(metric)%>%unique()%>%na.omit()%>%length()
#   prop <- 100*f/i
#   return <- list(i,f,prop)
#   names(return) = c('initial_number','filtered_number','proportion')
#   return(return)
# }
# ## genes
# get_prop(peakAnno,da_peaks_target_genes,'EnsemblID')
# # $initial_number
# # [1] 18334

# # $filtered_number
# # [1] 10088

# # $proportion
# # [1] 55.02

# ## peaks
# get_prop(peakAnno,da_peaks_target_genes,'peakID')
# # $initial_number
# # [1] 73432

# # $filtered_number
# # [1] 38336

# # $proportion
# # [1] 52.20612

# ## count the number of peaks associated per gene and stratify that by de/non-de  
# da_peaks_target_genes <- da_peaks_target_genes[,numb_peaks_per_gene:=.N,by=.(EnsemblID,DE)]

# pdf(paste0(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width = 7, height = 5)
# df <- copy(da_peaks_target_genes)[,c('DE','EnsemblID','numb_peaks_per_gene')]%>%unique()
# ggplot(da_peaks_target_genes,aes(x=DE,y=numb_peaks_per_gene,fill=DE))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# xlab(' ')+ylab('Number of peaks associated per gene')+
# scale_fill_manual(values=de_palette)+
# stat_compare_means(
#   method = "wilcox.test",
#   ref.group = 'non_de',
#   size=5
#   )+
#   theme_classic()+
#   theme(legend.position='none')
# dev.off()

# ## compare proportion DA peaks for de/non-de genes
# prop_da_peaks <- copy(da_peaks_target_genes)[
#   ,binned_dist:=ifelse(abs(distTSS)<=10000,10000,"> 10000")][
#     ,c('binned_dist','logFC','de_logFC','EnsemblID','peakID','DA','DE')
# ]%>%setorderv('binned_dist',1)%>%split(by='binned_dist')%>%lapply(function(x){
#   x<- x[
#   ,numb_peak_type_per_gene:=.N,by=.(EnsemblID,DA)
#   ][
#     ,totnumb_peaks_per_gene:=.N,by=.(EnsemblID)
#     ][
#       ,prop_peaktype_per_gene:=numb_peak_type_per_gene/totnumb_peaks_per_gene
#     ][
#       DA=='da'
#       ][
#       ,c('EnsemblID','numb_peak_type_per_gene','totnumb_peaks_per_gene','prop_peaktype_per_gene','DE')
#       ]%>%unique()
#   return(x)
# })
# prop_da_peaks <- Map(mutate,prop_da_peaks,binned_dist=names(prop_da_peaks))%>%rbindlist()
# prop_da_peaks$binned_dist <- factor(prop_da_peaks$binned_dist,levels=c('10000',"> 10000"))

# prop_da_peaks_stat=copy(prop_da_peaks)%>%split(by='binned_dist')%>%
# lapply(
#     function(x)x=compare_means(
#     prop_peaktype_per_gene~DE,
#     x,
#     method='wilcox.test', 
#     ref.group='non_de'
#     )
# )
# prop_da_peaks_stat = Map(mutate,prop_da_peaks_stat,binned_dist=names(prop_da_peaks_stat))%>%rbindlist()
# prop_da_peaks_stat = prop_da_peaks_stat[,c('group1','binned_dist','p.adj','p.signif')]%>%setnames(old=1,new='DE')

# pdf(paste0(outplot_dir,'prop_da_peaks_per_gene.pdf',sep=''),width = 7, height = 5)
# dodge <- position_dodge(width = 1)
# ggplot(prop_da_peaks,aes(x=binned_dist,y=prop_peaktype_per_gene,group = interaction(binned_dist,DE)))+
# geom_violin(trim=T,scale = "width",position = dodge,aes(fill = DE))+
# geom_boxplot(position = dodge,width=0.15,notch=T,fill = "white")+
# geom_text(data = prop_da_peaks_stat, aes(x = binned_dist, y = 1.1, label = p.signif),size=5)+
# xlab('Distance from TSS')+ylab('Proportion of DA peaks associated per gene')+
# scale_fill_manual(values=de_palette)+
# theme_classic()+
# theme(legend.position='bottom')
# dev.off()

## get target genes
peaks_w_genes <- fread(paste(target_genes_dir,'/target_genes.txt',sep=''),sep='\t',header=T)%>%na.omit()
# target_genes <- target_genes[da_results,on=c(range_keys,'peakID','DA'),nomatch=0]

# ## get DE genes 
# ## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
gene_exp_data <- fread(
  "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,de_logFC:=-logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','t','B','genes','logFC'):=NULL
]

peaks_w_genes_expr <- copy(peaks_w_genes)[gene_exp_data,on='EnsemblID',nomatch=0]

## get proportion peaks and genes retained
get_prop <- function(initial,filtered,metric){
  i <- copy(initial)%>%dplyr::pull(metric)%>%unique()%>%na.omit()%>%length()
  f <- copy(filtered)%>%dplyr::pull(metric)%>%unique()%>%na.omit()%>%length()
  prop <- 100*f/i
  return <- list(i,f,prop)
  names(return) = c('initial_number','filtered_number','proportion')
  return(return)
}
get_prop(peaks_w_genes,peaks_w_genes_expr,'EnsemblID')
# $initial_number
# [1] 21084

# $filtered_number
# [1] 11091

# $proportion
# [1] 52.60387
get_prop(peaks_w_genes,peaks_w_genes_expr,'peakID')
# $initial_number
# [1] 142875

# $filtered_number
# [1] 79015

# $proportion
# [1] 55.30359

## for each regulation type keep the closest peak to each gene
closest_peaks_w_genes_expr <- copy(peaks_w_genes_expr)%>%split(by='regulation')%>%lapply(
  function(x) x<-x[ ,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]
)%>%rbindlist()


## look at enrichment DA peaks near DE genes
## enrichment de genes with da peaks
significant_da_genes <- copy(closest_peaks_w_genes_expr)[DA=='da']$EnsemblID%>%unique()
significant_de_genes <- copy(closest_peaks_w_genes_expr)[DE=='de']$EnsemblID%>%unique()

da_de_overlap <-  list(significant_da_genes,significant_de_genes)

venn.diagram(
    x = da_de_overlap,
    category.names = c("DA","DE"),
    filename = paste(outplot_dir,'da_de_overlap.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col=c(da_palette[1],de_palette[1]),
    fill = c(alpha(da_palette[1],0.3),alpha(de_palette[1],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-20, 40),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = c(da_palette[1],de_palette[1])
)

## test significance overlap
da_de_genes <- Reduce(intersect,da_de_overlap)%>%unique()
nonda_de_genes <- copy(closest_peaks_w_genes_expr)[DA!='da'][DE=='de'][,EnsemblID]%>%unique()

da_nonde_genes <- copy(closest_peaks_w_genes_expr)[DA=='da'][DE!='de'][,EnsemblID]%>%unique()
nonda_nonde_genes <- copy(closest_peaks_w_genes_expr)[DA!='da'][DE!='de'][,EnsemblID]%>%unique()

matrix <- data.frame(
  as.numeric(length(da_de_genes)),
  as.numeric(length(nonda_de_genes)), 
  as.numeric(length(da_nonde_genes)),
  as.numeric(length(nonda_nonde_genes))
)%>%unlist()%>%matrix(nrow=2,byrow=T)

enrichment_dade_overlap <-fisher.test(matrix) 

## plot cumulative distribution number DA peaks over 10kb region
cumdist_peaks_genes = copy(closest_peaks_w_genes_expr)[,binned_dist:=plyr::round_any(abs(distTSS), 10)][abs(distTSS)<=100000]

pdf(paste(outplot_dir,'cumulative_distribution_peaks_gene_distances.pdf',sep=''),width = 7,height = 7)
ggplot(cumdist_peaks_genes, aes(binned_dist,col=DA)) + stat_ecdf()+
geom_vline(xintercept=10000,linetype='dashed')+
scale_color_manual(values=da_palette)+
theme_classic()
dev.off()


## add chrom state info to peaks 
# chrom_state_dir <- '../data/iPSC_chrom_states_hg38'

# ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
# setkeyv(ipsc_chromstate,range_keys)

# closest_peaks_w_genes_expr <- foverlaps(closest_peaks_w_genes_expr,ipsc_chromstate,type='any')

# closest_peaks_w_genes_expr <- closest_peaks_w_genes_expr[
#     ,overlap:=ifelse(i.start<start,i.end-start,end-i.start)
#     ][
#         ,.SD[which.max(overlap)], by=.(peakID)
#         ][
#             ,c(range_keys[-1],'overlap'):=NULL
# ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))


## correlation atac/rna-seq logFC
logfc_da_de_genes <- copy(closest_peaks_w_genes_expr)[,da_logFC:=logFC]

corr_dade=copy(logfc_da_de_genes)[DA=='da'][DE=='de']
cor.test(corr_dade$de_logFC,corr_dade$da_logFC)
# data:  corr_dade$de_logFC and corr_dade$logFC
# t = 5.7815, df = 159, p-value = 3.808e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2801710 0.5368505
# sample estimates:
#       cor 
# 0.4167843 
corr_nonda_nonde=copy(logfc_da_de_genes)[DA!='da'][DE!='de']
cor.test(corr_nonda_nonde$de_logFC,corr_nonda_nonde$logFC)
# 	Pearson's product-moment correlation
# data:  corr_nonda_nonde$de_logFC and corr_nonda_nonde$logFC
# t = 5.1355, df = 7049, p-value = 2.89e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.03776468 0.08427457
# sample estimates:
#        cor 
# 0.06105277 

count_genes_per_quadrant = function(x){
  maxlogfc = 10
  df = copy(x)[DA=='da'][DE=='de']%>%unique()
  df = df[
    ,quadrant:=ifelse((da_logFC>0 & de_logFC<0),'1',
      ifelse((da_logFC>0 & de_logFC>0),'2',
        ifelse((da_logFC<0 & de_logFC<0),'3','4'
        )))
        ][
          ,maxDAfc:=maxlogfc
          ][
            ,minDAfc:=-maxlogfc
            ][
              ,maxDEfc:=maxlogfc
              ][
                ,minDEfc:=-maxlogfc
                ][
                  ,c("EnsemblID",'quadrant','maxDAfc','minDAfc','maxDEfc','minDEfc')
                  ]%>%unique()
  df = df[
    ,numbgenes:=.N,by=.(quadrant)
    ][
      ,EnsemblID:=NULL
      ]%>%unique()%>%setorderv('quadrant',1)

  df = df[
    ,y:=ifelse((quadrant=='1'|quadrant=='2'),maxDAfc,minDAfc)
    ][
      ,x:=ifelse((quadrant=='1'|quadrant=='3'),minDEfc,maxDEfc)
      ][,dade_overlap:='overlap']
  return(df)
}

numb_genes <- count_genes_per_quadrant(logfc_da_de_genes)

logfc_da_de_genes <-logfc_da_de_genes[
  ,dade_overlap:=as.factor(ifelse((DA=='da'& DE=='de'),'overlap','no_overlap'))
  ][
    ,c('dade_overlap','de_logFC','da_logFC')
]

pdf(paste(outplot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 7,height = 7)
ggplot(logfc_da_de_genes, aes(x=de_logFC, y=da_logFC,col=dade_overlap)) +
geom_point(alpha=0.2)+
scale_color_manual(values=c('darkgray','red'))+
geom_text(data = numb_genes, aes(label = numb_genes$numbgenes,y = numb_genes$y, x = numb_genes$x))+
geom_vline(xintercept=0,linetype='dashed')+
geom_hline(yintercept=0,linetype='dashed')+
geom_smooth(se=F,show.legend=T)+
ylab('DA logFC')+xlab('DE logFC')+
xlim(-10,10)+ylim(-10,10)+
theme_classic()+
theme(
  legend.position ='none'
)
dev.off()




# ## look @ the overlap between set of orthologous genes targeted by DA and non-DA peaks
# genes <- copy(peaks_w_orthgenes)%>%split(by='DA')%>%lapply(function(x)x=unique(na.omit(x$EnsemblID)))
# venn.diagram(
#     x = genes,
#     category.names = c('da','non_da'),
#     filename = paste(outplot_dir,'target_genes_venn.png',sep=''),
#     output = TRUE ,
#     imagetype="png" ,
#     height = 700 , 
#     width = 700 , 
#     resolution = 400,
#     lwd = 1,
#     col=da_palette,
#     fill = c(alpha(da_palette[[1]],0.3), alpha(da_palette[[2]],0.3)),
#     cex = 0.5,
#     fontfamily = "sans",
#     cat.cex = 0.3,
#     cat.default.pos = "outer",
#     cat.pos = c(-27, 27),
#     cat.dist = c(0.055, 0.055),
#     cat.fontfamily = "sans",
#     cat.col = da_palette
# )

# ## number of peaks (DA/non-DA) associated to each gene
# peaks_w_orthgenes<- peaks_w_orthgenes[,numb_da_peaks_per_gene:=.N,by=.(EnsemblID,DA)]

# ## number peaks per gene
# pdf(paste(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width=7,height = 5)
# df <- copy(peaks_w_orthgenes)[,c('EnsemblID','DA','peakID')]%>%na.omit()%>%unique()
# df <- df[,numb_peaks_per_gene:=.N,by=.(DA,EnsemblID)]
# ggplot(df,aes(x=numb_peaks_per_gene,fill=DA))+
# geom_bar(position=position_dodge())+
# scale_fill_manual(values=da_palette)+
# facet_wrap(DA~.,ncol=2)+
# theme_classic()+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black")
# )
# dev.off()


# ## further keep peaks with genes in same tad 
# ## get orthologous tads
# orth_tads <- read_tads('human_chimp_orth_tads.txt')[,c(1:3)]
# colnames(orth_tads) = range_keys
# setkeyv(orth_tads,range_keys)

# peaks_w_ensembl_in_tad<- foverlaps(copy(peaks_w_genes_ensembl),orth_tads,type='within')%>%na.omit()%>%dplyr::select(-c(contains('per_gene')))
# peaks_w_ensembl_in_tad <- peaks_w_ensembl_in_tad[
#   ,same_tad:=ifelse(distTSS<0 & start- (i.start+distTSS)>0,'y',
#   ifelse(distTSS>0 & end - (i.end+distTSS)>0,'y','n')
#   )
#   ][
#     same_tad=='y'
#     ][
#       ,c(range_keys[-1]):=NULL
# ]%>%unique()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

# ## prop peaks & genes retained
# ## count numb and prop of peaks and genes retained 
# numb_peaks_genes <- copy(peaks_w_genes)[,c('DA','peakID')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))
# numb_peaks_ensembl_tad <- copy(peaks_w_ensembl_in_tad)[,c('DA','peakID')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))

# purrr::map2(numb_peaks_ensembl_tad,numb_peaks_genes,`/`)

# numb_genes <- copy(peaks_w_genes)[,c('DA','gene')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))
# numb_genes_ensembl_tad <- copy(peaks_w_ensembl_in_tad)[,c('DA','gene')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))

# purrr::map2(numb_genes_ensembl_tad,numb_genes,`/`)

##--------------------------
## now focus on differences 
## between DE/non-DE genes
##--------------------------
# ## check how many peaks each gene has associated and stratify that by de/non-de  
# peaks_w_orthgenes <- peaks_w_orthgenes[,numb_peaks_per_gene_expr:=.N,by=.(EnsemblID,DE)]

# pdf(paste0(outplot_dir,'numb_peaks_per_gene_expr.pdf',sep=''),width = 7, height = 5)
# df <- copy(peaks_w_orthgenes)[,c('DE','EnsemblID','numb_peaks_per_gene_expr')]%>%unique()
# ggplot(peaks_w_orthgenes,aes(x=DE,y=numb_peaks_per_gene_expr,fill=DE))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# xlab(' ')+ylab('Number of peaks associated per gene')+
# scale_fill_manual(values=de_palette)+
# stat_compare_means(
#   method = "wilcox.test",
#   ref.group = 'non_de',
#   size=5
#   )+
#   theme_classic()+
#   theme(legend.position='none')
# dev.off()

# ## compare proportion da peaks for de/non-de genes
# prop_da_peaks <- copy(peaks_w_orthgenes)[
#   ,binned_dist:=ifelse(abs(distTSS)<=10000,10000,"> 10000")][
#     ,c('binned_dist','logFC','de_logFC','EnsemblID','peakID','DA','DE')
# ]%>%setorderv('binned_dist',1)%>%split(by='binned_dist')%>%lapply(function(x){
#   x<- x[
#   ,numb_peak_type_per_gene:=.N,by=.(EnsemblID,DA)
#   ][
#     ,totnumb_peaks_per_gene:=.N,by=.(EnsemblID)
#     ][
#       ,prop_peaktype_per_gene:=numb_peak_type_per_gene/totnumb_peaks_per_gene
#     ][
#       DA=='da'
#       ][
#       ,c('EnsemblID','numb_peak_type_per_gene','totnumb_peaks_per_gene','prop_peaktype_per_gene','DE')
#       ]%>%unique()
#   return(x)
# })
# prop_da_peaks <- Map(mutate,prop_da_peaks,binned_dist=names(prop_da_peaks))%>%rbindlist()
# prop_da_peaks$binned_dist <- factor(prop_da_peaks$binned_dist,levels=c('10000',"> 10000"))

# prop_da_peaks_stat=copy(prop_da_peaks)%>%split(by='binned_dist')%>%
# lapply(
#     function(x)x=compare_means(
#     prop_peaktype_per_gene~DE,
#     x,
#     method='wilcox.test', 
#     ref.group='non_de'
#     )
# )
# prop_da_peaks_stat = Map(mutate,prop_da_peaks_stat,binned_dist=names(prop_da_peaks_stat))%>%rbindlist()
# prop_da_peaks_stat = prop_da_peaks_stat[,c('group1','binned_dist','p.adj','p.signif')]%>%setnames(old=1,new='DE')

# pdf(paste0(outplot_dir,'prop_da_peaks_per_gene.pdf',sep=''),width = 7, height = 5)
# dodge <- position_dodge(width = 1)
# ggplot(prop_da_peaks,aes(x=binned_dist,y=prop_peaktype_per_gene,group = interaction(binned_dist,DE)))+
# geom_violin(trim=T,scale = "width",position = dodge,aes(fill = DE))+
# geom_boxplot(position = dodge,width=0.15,notch=T,fill = "white")+
# geom_text(data = prop_da_peaks_stat, aes(x = binned_dist, y = 1.1, label = p.signif),size=5)+
# xlab('Distance from TSS')+ylab('Proportion of DA peaks associated per gene')+
# scale_fill_manual(values=de_palette)+
# theme_classic()+
# theme(legend.position='bottom')
# dev.off()

##---------------------------
## 1 peak per gene
##---------------------------
# one_peak_per_gene <- copy(peaks_w_orthgenes)[,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]

# ## enrichment de genes with da peaks
# significant_da_genes <- copy(one_peak_per_gene)[DA=='da']$EnsemblID%>%unique()
# significant_de_genes <- copy(one_peak_per_gene)[DE=='de']$EnsemblID%>%unique()

# da_de_overlap <-  list(significant_da_genes,significant_de_genes)

# venn.diagram(
#     x = da_de_overlap,
#     category.names = c("DA","DE"),
#     filename = paste(outplot_dir,'da_de_overlap.png',sep=''),
#     output = TRUE ,
#     imagetype="png" ,
#     height = 700 , 
#     width = 700 , 
#     resolution = 400,
#     lwd = 1,
#    col=c(da_palette[1],de_palette[1]),
#     fill = c(alpha(da_palette[1],0.3),alpha(de_palette[1],0.3)),
#     cex = 0.5,
#     fontfamily = "sans",
#     cat.cex = 0.3,
#     cat.default.pos = "outer",
#     cat.pos = c(-20, 40),
#     cat.dist = c(0.055, 0.055),
#     cat.fontfamily = "sans",
#     cat.col = c(da_palette[1],de_palette[1])
# )

# ## test significance overlap
# da_de_genes <- Reduce(intersect,da_de_overlap)%>%unique()
# nonda_de_genes <- copy(one_peak_per_gene)[DA!='da'][DE=='de'][,EnsemblID]%>%unique()

# da_nonde_genes <- copy(one_peak_per_gene)[DA=='da'][DE!='de'][,EnsemblID]%>%unique()
# nonda_nonde_genes <- copy(one_peak_per_gene)[DA!='da'][DE!='de'][,EnsemblID]%>%unique()

# matrix <- data.frame(
#   as.numeric(length(da_de_genes)),
#   as.numeric(length(nonda_de_genes)), 
#   as.numeric(length(da_nonde_genes)),
#   as.numeric(length(nonda_nonde_genes))
# )%>%unlist()%>%matrix(nrow=2,byrow=T)

# enrichment_dade_overlap <-fisher.test(matrix) 
# p-value = 0.0001723
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.202687 1.845721
# sample estimates:
# odds ratio 
#   1.490688 

##------------------------------------------------------
## look at avg expression genes by peaks (DA vs non-DA)
##------------------------------------------------------
# all_genes_expr <- copy(de_genes)[,c('AveExpr')][,DA:='all_genes']

# comparisons = list(
#   c('non_da','da'),
#   c('da','all_genes'),
#   c('non_da','all_genes')
# )

# plot_de_genes <- function(df,xaxis,yaxis,ylab){
#   df_plot <- copy(df)[
#     ,xcol:=xaxis
#     ][
#       ,ycol:=yaxis
#       ]

#   p <- ggplot(df_plot,aes(x=xcol,y=ycol,fill=xcol))+
#   geom_violin(trim=T,scale = "width")+
#   geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
#   xlab(' ')+ylab(ylab)+
#   stat_compare_means(
#   method = "wilcox.test",
#   comparisons = comparisons,
#   size=5
#   )+
#     theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
#     )
#   return(p)
# }

# pdf(paste0(outplot_dir,'avg_expr_genes_by_peak_significance.pdf',sep=''),width = 10, height = 8)
# df <- rbind(one_peak_per_gene[,c('AveExpr','DA')],all_genes_expr)
# plot_de_genes(df,df$DA,df$AveExpr,'AveExpr')
# dev.off()


## get loeuf score for these genes
# gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
# colnames(gnomad_loeuf)[1] = 'EnsemblID'

# ## add loeuf score to genes
# one_peak_per_gene_loeuf <- copy(one_peak_per_gene)[gnomad_loeuf,on='EnsemblID',nomatch=0]

# all_genes_loeuf <- copy(gnomad_loeuf)[,DA:='all_genes']
# peak_genes_loeuf <- copy(one_peak_per_gene_loeuf)%>%dplyr::select(c(all_of(colnames(all_genes_loeuf))))

# pdf(paste0(outplot_dir,'loeuf_score_genes_by_peak.pdf',sep=''),width = 10, height = 8)
# df <- rbind(all_genes_loeuf,peak_genes_loeuf)
# plot_de_genes(df,df$DA,df$oe_lof_upper,'oe_lof_upper')
# dev.off()


# logfc_da_de_genes <- copy(one_peak_per_gene)[,da_logFC:=logFC]

# corr_dade=copy(logfc_da_de_genes)[DA=='da'][DE=='de']
# cor.test(corr_dade$de_logFC,corr_dade$logFC)
# # data:  corr_dade$de_logFC and corr_dade$logFC
# # t = 5.7815, df = 159, p-value = 3.808e-08
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.2801710 0.5368505
# # sample estimates:
# #       cor 
# # 0.4167843 
# corr_nonda_nonde=copy(logfc_da_de_genes)[DA!='da'][DE!='de']
# cor.test(corr_nonda_nonde$de_logFC,corr_nonda_nonde$logFC)
# # 	Pearson's product-moment correlation
# # data:  corr_nonda_nonde$de_logFC and corr_nonda_nonde$logFC
# # t = 5.1355, df = 7049, p-value = 2.89e-07
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  0.03776468 0.08427457
# # sample estimates:
# #        cor 
# # 0.06105277 

# count_genes_per_quadrant = function(x){
#   maxlogfc = 10
#   df = copy(x)[DA=='da'][DE=='de']%>%unique()
#   df = df[
#     ,quadrant:=ifelse((da_logFC>0 & de_logFC<0),'1',
#       ifelse((da_logFC>0 & de_logFC>0),'2',
#         ifelse((da_logFC<0 & de_logFC<0),'3','4'
#         )))
#         ][
#           ,maxDAfc:=maxlogfc
#           ][
#             ,minDAfc:=-maxlogfc
#             ][
#               ,maxDEfc:=maxlogfc
#               ][
#                 ,minDEfc:=-maxlogfc
#                 ][
#                   ,c("EnsemblID",'quadrant','maxDAfc','minDAfc','maxDEfc','minDEfc')
#                   ]%>%unique()
#   df = df[
#     ,numbgenes:=.N,by=.(quadrant)
#     ][
#       ,EnsemblID:=NULL
#       ]%>%unique()%>%setorderv('quadrant',1)

#   df = df[
#     ,y:=ifelse((quadrant=='1'|quadrant=='2'),maxDAfc,minDAfc)
#     ][
#       ,x:=ifelse((quadrant=='1'|quadrant=='3'),minDEfc,maxDEfc)
#       ][,dade_overlap:='overlap']
#   return(df)
# }

# numb_genes <- count_genes_per_quadrant(logfc_da_de_genes)

# logfc_da_de_genes <-logfc_da_de_genes[
#   ,dade_overlap:=as.factor(ifelse((DA=='da'& DE=='de'),'overlap','no_overlap'))
#   ][
#     ,c('dade_overlap','de_logFC','da_logFC')
# ]

# pdf(paste(outplot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 7,height = 7)
# ggplot(logfc_da_de_genes, aes(x=de_logFC, y=da_logFC,col=dade_overlap)) +
# geom_point(alpha=0.2)+
# scale_color_manual(values=c('darkgray','red'))+
# geom_text(data = numb_genes, aes(label = numb_genes$numbgenes,y = numb_genes$y, x = numb_genes$x))+
# geom_vline(xintercept=0,linetype='dashed')+
# geom_hline(yintercept=0,linetype='dashed')+
# geom_smooth(se=F,show.legend=T)+
# ylab('DA logFC')+xlab('DE logFC')+
# xlim(-10,10)+ylim(-10,10)+
# theme_classic()+
# theme(
#   legend.position ='none'
# )
# dev.off()

# ##---------------------------------------------------
# ## look at patterns accessibility around DE genes
# ##---------------------------------------------------
# library(edgeR)
# library(csaw)
# library(ComplexHeatmap)

# peak_file <- paste('output/final_peak_set/',genome,'_all_orthologous_peaks.txt',sep='')

# ## list bam
# standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
# param <- readParam(pe = "both",restrict=standard_chr,max.frag=1000)

# get_bams <- function(species){
#     bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
#     return(bams)
# }

# bams <- get_bams(genome)

# ## read consensus peak
# peaks <- fread(peak_file,sep='\t',header=T)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

# numb_all_peaks <- length(peaks)

# ##---------------------
# ## QC & Normalisation
# ##---------------------
# ## Use RPKM to account for differences in peak sizes
# reads_in_peaks <- regionCounts(bams, peaks,param=param)%>%asDGEList()
# peak_width <- copy(peaks)%>%as.data.table()%>%dplyr::pull('width')
# reads_in_peaks$counts <- rpkm(reads_in_peaks$counts,peak_width)

# rownames(reads_in_peaks) <- as.data.table(copy(peaks))$peakID
# colnames(reads_in_peaks) <- samples_names
# reads_in_peaks$samples <- mutate(reads_in_peaks$samples,'group'=ifelse(rownames(reads_in_peaks$samples) %like% 'C','chimp','human'))

# ## filter low reads 
# keep <- filterByExpr(reads_in_peaks, group=reads_in_peaks$group,min.total.count = 3,min.count=3)
# reads_in_peaks_filtered <- copy(reads_in_peaks)[keep,, keep.lib.sizes=FALSE]

# ## quantile normalise read counts
# reads_in_peaks_filtered_tmm <- calcNormFactors(copy(reads_in_peaks_filtered), method = "TMM")
# quant_norm_counts <- copy(reads_in_peaks_filtered_tmm)
# quant_norm_counts$counts <- normalizeQuantiles(reads_in_peaks_filtered_tmm$counts, ties=F)

# ## get row DA matrix corresponding to DA peaks associated with DE genes
# da_peak_de_genes <- copy(one_peak_per_gene)[DA=='da' & DE=='de']

# subset_peak_count_matrix <- copy(quant_norm_counts)
# subset_peak_count_matrix$counts <- subset(subset_peak_count_matrix$counts,rownames(subset_peak_count_matrix$counts) %in% da_peak_de_genes$peakID)

# peak_counts_df <- data.table(subset_peak_count_matrix$counts)[
#     ,peakID:=rownames(subset_peak_count_matrix$counts)
# ]

## add chrom state info to peaks 
chrom_state_dir <- '../data/iPSC_chrom_states_hg38'

ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

da_de_peaks <-copy(one_peak_per_gene)[EnsemblID %in% da_de_genes]

da_peak_de_genes_chromstate <- foverlaps(copy(da_de_peaks),ipsc_chromstate,type='any')[
    ,overlap:=ifelse(i.start<start,i.end-start,end-i.start)
    ][
        ,.SD[which.max(overlap)], by=.(peakID)
        ][
            ,c(range_keys[-1],'overlap'):=NULL
]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

# peak_counts_df <- peak_counts_df[da_peak_de_genes_chromstate,on='peakID',nomatch=0]%>%dplyr::select(c(1:12,'peak_chromstate'))

# mat_rownames = peak_counts_df$peak_chromstate
# peakcounts_mat <- peak_counts_df[,peak_chromstate:=NULL]%>%as.matrix()
# rownames(peakcounts_mat)=mat_rownames

# ## heatmap
# pdf(paste(outplot_dir,"heatmap_quantnormd_counts_peaks_w_de_genes.pdf",sep=''),width=10,height=7)
# x=round(log2(t(peakcounts_mat)+1),1)
# y=table(x)
# x[x > 5.4] <- 6
# names(chrom_state_colors) = chrom_states

# Heatmap(
#     x,
#     show_column_dend = F,
#     show_column_names = F,
#     column_order = order(as.factor(readr::parse_number(gsub("^.*\\.", "",colnames(x))))),
#     # column_split = as.factor(readr::parse_number(gsub("^.*\\.", "", colnames(x)))),
#     top_annotation = HeatmapAnnotation(
#     chrom_state = anno_simple(
#         gsub("^.*\\.", "", colnames(x)),
#         border=T,
#         height = unit(1,'cm'),
#         col=chrom_state_colors),
#         show_annotation_name = F
#     ),
#     name='log2 rpkm',
#     col=viridis::viridis(7)
# )
# dev.off()


## write file with all info to make visualization script much neater/quicker
write.table(da_peak_de_genes_chromstate,paste(tmp_files_dir,'da_peaks_w_de_genes_gviz.txt',sep=''),sep='\t',col.names=T,quote=F,row.names = F)