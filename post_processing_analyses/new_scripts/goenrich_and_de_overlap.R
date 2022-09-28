## script used to perform GO enrichment of common and species specific DA peaks using GREAT
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ChIPseeker)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'GO_enrich')
outfile_dir <- create_dir(outdir,'files/GO_enrich')

de_palette <- c('#DDA15E','#606C38')
names(de_palette) = c('de','non_de')

## get DA peaks
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

## use chipseeker to annotate peaks within 1Mb of distance to their putative target gene(s)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tss_region = 10000
peakAnno <- annotatePeak(
  makeGRangesFromDataFrame(da_results,keep.extra.columns=T), 
  tssRegion=c(-tss_region, tss_region),
  TxDb = txdb,
  level = "gene",
  overlap='all'
)%>%as.data.table()

## get genes IDs
symbols <- mapIds(org.Hs.eg.db, keys = peakAnno$geneId, keytype = "ENTREZID", column="SYMBOL")%>%na.omit()
symbols_df <- data.table(symbol=symbols,entrezid=names(symbols))%>%unique()
ensemblid <-mapIds(org.Hs.eg.db, keys = symbols_df$symbol, keytype = "SYMBOL", column="ENSEMBL")%>%na.omit()
ensemblid_df <- data.table(symbol=names(ensemblid),EnsemblID=ensemblid)

gene_ids <- copy(ensemblid_df)[symbols_df,on='symbol',nomatch=0]%>%setnames(old='entrezid',new='geneId')
peakAnno <- copy(peakAnno)[gene_ids,on='geneId',nomatch=0][,c('width','strand','geneLength','geneStrand','geneId','annotation'):=NULL]

## look @ the overlap between genes targeted by DA and non-DA peaks 
library(VennDiagram)
 
genes <- copy(peakAnno)%>%split(by='DA')%>%lapply(function(x)x=x[,EnsemblID]%>%unique()%>%na.omit())
venn.diagram(
    x = genes,
    category.names = c('da','non_da'),
    filename = paste(outplot_dir,'target_genes_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col=da_palette,
    fill = c(alpha(da_palette[[1]],0.3), alpha(da_palette[[2]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = da_palette
)

## some numbs and props 
numb_peaks_associated <- copy(peakAnno)[,c('peakID','DA','EnsemblID')]%>%unique()%>%na.omit()
numb_peaks_associated <-numb_peaks_associated[,EnsemblID:=NULL][,n:=.N,by=.(DA)][,c('DA','n')]%>%unique()
tot_numb_peaks <- copy(da_results)[,c('peakID','DA')][,tot_n:=.N,by=.(DA)][,c('DA','tot_n')]%>%unique()
prop <- numb_peaks_associated[tot_numb_peaks,on='DA',nomatch=0][,prop:=n/tot_n]
#       DA     n  tot_n    prop
# 1: non_da 93334 98914 0.9435874
# 2:     da 49869 52352 0.9525711

## Distribution number peaks per gene
pdf(paste(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width=7,height = 7)
df <- copy(peakAnno)[,c('EnsemblID','DA','peakID')]%>%na.omit()%>%unique()
df <- df[,numb_peaks_per_gene:=.N,by=.(DA,EnsemblID)][!numb_peaks_per_gene>100]
ggplot(df,aes(x=DA,y=numb_peaks_per_gene,fill=DA))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
scale_fill_manual(values=da_palette)+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'da',
  size=5
  )+
theme_classic()+
theme(
  legend.position ='none'
)
dev.off()

## get median number of peaks per gene
median(df[DA=='da']$numb_peaks_per_gene)
# [1] 5
median(df[DA!='da']$numb_peaks_per_gene)
# [1] 13

## assign peaks to genes
## group peaks into da/nonda only and da+nonda
peaks_w_genes <- copy(peakAnno)[,distTSS:=start-geneStart][,c(..range_keys,'DA','da_species','peakID','FDR','logFC','EnsemblID','symbol','distanceToTSS')]
gene_regulation <- copy(peaks_w_genes)[,c('DA','EnsemblID')]%>%unique()
gene_regulation <- gene_regulation[,numbpeak_types:=.N,by=.(EnsemblID)][,regulation:=ifelse(numbpeak_types>1,'nonda_da',DA)][,c('EnsemblID','regulation')]%>%unique()
peaks_w_genes <- peaks_w_genes[gene_regulation,on='EnsemblID',nomatch=0]

##---------------------------------------------
## GO enrichment using clusterProfiler
##---------------------------------------------
library(clusterProfiler)

genes_by_regulation <- copy(peaks_w_genes)[,c('regulation','EnsemblID')]%>%unique()%>%split(by='regulation')%>%lapply(function(x)x$EnsemblID)
compGOs <- compareCluster(
  geneCluster = genes_by_regulation,
  ont = "BP",
  fun = "enrichGO",
  keyType = "ENSEMBL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "fdr",
  OrgDb='org.Hs.eg.db'
)
pdf(paste(outplot_dir,'go_enrichment.pdf',sep=''),width=10,height = 15)
dotplot(compGOs, showCategory = 10, title = "GO Enrichment Analysis")
dev.off()

# ##---------------------------------------------
# ## loeuf score to check whether genes are
# ## less loss-of-function tolerant 
# ##---------------------------------------------
# gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
# colnames(gnomad_loeuf)[1] = 'EnsemblID'

# loeuf_score_mygenes <- copy(peaks_w_genes)[gnomad_loeuf,on='EnsemblID',nomatch=0][,type:=regulation][,c('type','EnsemblID','oe_lof_upper')]
# loeuf_score_all_genes <- copy(gnomad_loeuf)[,type:='all_genes'][,c('type','EnsemblID','oe_lof_upper')]

# all_loeuf_scores <- rbind(loeuf_score_mygenes,loeuf_score_all_genes)%>%na.omit()

# comparisons = list(
#   c('da','all_genes'),
#   c('non_da','all_genes'),
#   c('nonda_da','all_genes'),
#   c('non_da','da'),
#   c('nonda_da','da'),
#   c('nonda_da','non_da')
# )

# pdf(paste0(outplot_dir,'loeuf_score_target_genes.pdf',sep=''),width = 10, height = 8)
# loeuf_palette = c('#edf2f4','red',da_palette)
# names(loeuf_palette) = c('all_genes','nonda_da',names(da_palette))

# ggplot(all_loeuf_scores,aes(x=type,y=oe_lof_upper,fill=type))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# ylab('LOEUF score')+xlab(' ')+
# scale_fill_manual(values=loeuf_palette)+
# stat_compare_means(
#   method = "wilcox.test",
#   comparisons = comparisons,
#   size=5
#   )+
#   theme_classic()+
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
# )
# dev.off()


# # closest_peak <- copy(peaks_w_genes[regulation=='nonda_da'])[,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]


# # peaks_w_target_genes <- rbind(copy(closest_peak),copy(peaks_w_genes[regulation!='nonda_da']))[,regulation:=NULL]

##---------------
## export files
##---------------
## target genes
write.table(peaks_w_genes, file = paste(outfile_dir,'target_genes.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## GO tables
go_table <- copy(compGOs@compareClusterResult)%>%as.data.table()
go_table <- go_table[,geneID:=NULL]
write.xlsx(go_table,paste(outfile_dir,'GO_enrich_terms_da_non_da.xlsx',sep=''),append=T,overwrite=T)

##---------------------------------
## look at overlap with DE
##---------------------------------
## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
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

## some number and proportions
length(unique(peaks_w_genes_expr$EnsemblID))
# [1] 10667
length(unique(peaks_w_genes_expr$EnsemblID))/length(unique(peaks_w_genes$EnsemblID))
# [1] 0.5607717

length(unique(peaks_w_genes_expr[DE=='de']$EnsemblID))
# [1] 3466
length(unique(peaks_w_genes_expr[DE=='de']$EnsemblID))/length(unique(peaks_w_genes_expr$EnsemblID))
# [1] 0.3249273

## get number of peaks associated for DE/non-DE genes
pdf(paste(outplot_dir,'numb_peaks_per_de_gene.pdf',sep=''),width=7,height = 7)
df <- copy(peaks_w_genes_expr)[,c('EnsemblID','DE','peakID')]%>%na.omit()%>%unique()
df <- df[,numb_peaks_per_gene:=.N,by=.(DE,EnsemblID)][!numb_peaks_per_gene>100]
ggplot(df,aes(x=DE,y=numb_peaks_per_gene,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
scale_fill_manual(values=de_palette)+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'de',
  size=5
  )+
theme_classic()+
theme(
  legend.position ='none'
)
dev.off()

## for each regulation type keep the closest peak to each gene
closest_peaks_w_genes_expr <- copy(peaks_w_genes_expr)%>%split(by='regulation')%>%lapply(
  function(x) x<-x[ ,.SD[which.max(abs(logFC))], by=.(EnsemblID)]
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
  as.numeric(length(da_nonde_genes)),
  as.numeric(length(nonda_de_genes)), 
  as.numeric(length(nonda_nonde_genes))
)%>%unlist()%>%matrix(nrow=2,byrow=T)

enrichment_dade_overlap <-fisher.test(matrix,alternative='g')


## correlation atac/rna-seq logFC
# logfc_da_de_genes <- copy(closest_peaks_w_genes_expr)[,da_logFC:=logFC]
logfc_da_de_genes <- copy(closest_peaks_w_genes_expr)[,da_logFC:=logFC][DA=='da'][DE=='de']

cor.test(logfc_da_de_genes$de_logFC,logfc_da_de_genes$da_logFC)
# 	Pearson's product-moment correlation

# data:  logfc_da_de_genes$de_logFC and logfc_da_de_genes$da_logFC
# t = 1.2772, df = 3072, p-value = 0.2016
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.01232589  0.05834277
# sample estimates:
#        cor 
# 0.02303722 

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
ggplot(logfc_da_de_genes, aes(x=de_logFC, y=da_logFC)) +
geom_point(alpha=0.2)+
# scale_color_manual(values=c('darkgray','red'))+
geom_text(data = numb_genes, aes(label = numb_genes$numbgenes,y = numb_genes$y, x = numb_genes$x))+
geom_vline(xintercept=0,linetype='dashed')+
geom_hline(yintercept=0,linetype='dashed')+
ylab('DA logFC')+xlab('DE logFC')+
xlim(-10,10)+ylim(-10,10)+
theme_classic()+
theme(
  legend.position ='none'
)
dev.off()




