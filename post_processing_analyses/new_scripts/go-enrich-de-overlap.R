## script used to perform GO enrichment of common and species specific DA peaks using ChIPseeker
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ChIPseeker)
library(VennDiagram)
library(clusterProfiler)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'goEnrich')
outfile_dir <- create_dir(outdir,'files/goEnrich')

de_palette <- c('#DDA15E','#606C38')
names(de_palette) = c('de','non_de')

## get DA peaks
da_results <- read_da_results('da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

## split by da_species
# da_results_split = split(copy(da_results),by='da_species')

## use chipseeker to annotate peaks within 1Mb of distance to their putative target gene(s)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tss_region = 10000

peakAnno <- annotatePeak(
  makeGRangesFromDataFrame(copy(da_results),keep.extra.columns=T), 
  tssRegion=c(-tss_region, tss_region),
  TxDb = txdb,
  sameStrand=T,
  level = "gene",
  overlap='all'
)%>%as.data.table()

## get genes IDs
symbols <- mapIds(org.Hs.eg.db, keys = peakAnno$geneId, keytype = "ENTREZID", column="SYMBOL")%>%na.omit()
symbols_df <- data.table(symbol=symbols,entrezid=names(symbols))%>%unique()
ensemblid <- mapIds(org.Hs.eg.db, keys = symbols_df$symbol, keytype = "SYMBOL", column="ENSEMBL")%>%na.omit()
ensemblid_df <- data.table(symbol=names(ensemblid),EnsemblID=ensemblid)

gene_ids <- copy(ensemblid_df)[symbols_df,on='symbol',nomatch=0]%>%setnames(old='entrezid',new='geneId')
peakAnno <- copy(peakAnno)[gene_ids,on='geneId',nomatch=0][,c('width','strand','geneLength','geneStrand','geneId','annotation'):=NULL]

## proportion of peaks that have been associated with genes
nrow(peakAnno)
# [1] 76661
nrow(da_results)
# 1] 81065
nrow(peakAnno)/nrow(da_results)*100
# [1] 94.56732

## look @ the overlap between genes targeted by DA and non-DA peaks  
genes <- lapply(split(copy(peakAnno),by='da_species'),function(x)x=x[,EnsemblID]%>%unique()%>%na.omit())

venn.diagram(
    x = genes,
    category.names = names(genes),
    filename = paste(outplot_dir,'venn-overlap-target-genes.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col = da_species_palette, 
    fill = c(alpha(da_species_palette[[1]],0.3),alpha(da_species_palette[[3]],0.3),alpha(da_species_palette[[2]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = 'black'
)

## assign peaks to genes based on strongest peak (logfc)
peakAnno_strongest_peak <- copy(peakAnno)[ ,.SD[which.max(abs(logFC))], by=.(EnsemblID)][,bindist:=plyr::round_any(distanceToTSS, 1000)]

# regulation <- rbindlist(copy(peakAnno))[,c('DA','EnsemblID','da_species')]%>%unique()
# regulation = regulation[,npeaktype:=.N,by=.(EnsemblID)][,regulation:=ifelse(npeaktype==2,'da',DA)][,remove:=ifelse( (DA=='non_da'& regulation=='da'),'yes','no')][remove=='no'][,remove:=NULL]

# peakAnno_closest_peak <- rbindlist(copy(peakAnno))[regulation, on=c('EnsemblID','DA','da_species'),nomatch=0][,npeaktype:=NULL][,.SD[which.min(abs(distanceToTSS))], by=.(EnsemblID)]

## distribution peak-tss distances
pdf(paste0(outplot_dir,'qc-peak-tss-distance.pdf',sep=''),width = 10, height = 8)
df = copy(peakAnno_strongest_peak)[,c('EnsemblID','peakID','da_species','distanceToTSS')][abs(bindist)<=10000]%>%unique()

ggplot(df,aes(x=bindist,fill=da_species))+
geom_density(alpha=0.5)+
ylab('Density')+xlab(' ')+
scale_fill_manual(values=da_species_palette)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

nrow(peakAnno_strongest_peak[da_species=='chimp'])
nrow(peakAnno_strongest_peak[da_species=='human'])
nrow(peakAnno_strongest_peak[da_species=='non_da'])

##---------------------------------------------
## GO using clusterProfiler
##---------------------------------------------
genes_by_species <- copy(peakAnno_strongest_peak)[,c('da_species','EnsemblID')]%>%unique()%>%split(by='da_species')%>%lapply(function(x)x$EnsemblID)
compGOs <- compareCluster(
  geneCluster = genes_by_species,
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

##---------------------------------------------
## loeuf score to check whether genes are
## less loss-of-function tolerant 
##---------------------------------------------
gnomad_loeucompGOsf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
colnames(gnomad_loeuf)[1] = 'EnsemblID'

loeuf_score_mygenes <- copy(peakAnno_strongest_peak)[gnomad_loeuf,on='EnsemblID',nomatch=0][,type:=da_species][,c('type','EnsemblID','oe_lof_upper')]
loeuf_score_all_genes <- copy(gnomad_loeuf)[,type:='all_genes'][,c('type','EnsemblID','oe_lof_upper')]

all_loeuf_scores <- rbind(loeuf_score_mygenes,loeuf_score_all_genes)%>%na.omit()

comparisons = list(
  c('chimp','all_genes'),
  c('human','all_genes'),
  c('non_da','all_genes'),
  c('chimp','non_da'),
  c('human','non_da')
)

pdf(paste0(outplot_dir,'loeuf-score-target-genes.pdf',sep=''),width = 10, height = 8)
loeuf_palette = c('#edf2f4',da_species_palette)
names(loeuf_palette) = c('all_genes',names(species_palette))

ggplot(all_loeuf_scores,aes(x=type,y=oe_lof_upper,fill=type))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
ylab('LOEUF score')+xlab(' ')+
scale_fill_manual(values=loeuf_palette)+
stat_compare_means(
  method = "wilcox.test",
  comparisons = comparisons,
  size=5
  )+
  theme_classic()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

##---------------
## export files
##---------------
## target genes
write.table(peakAnno_strongest_peak, file = paste(outfile_dir,'peak-target-genes.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## GO tables
go_table <- copy(compGOs@compareClusterResult)%>%as.data.table()
go_table <- go_table[,geneID:=NULL]
write.xlsx(go_table,paste(outfile_dir,'go-enrich-terms-da-vs-nonda.xlsx',sep=''),append=T,overwrite=T)

##---------------------------------
## look at overlap with DE
##---------------------------------
## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
genExpr <- fread(
  "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,de_logFC:=-logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','t','B','genes','logFC'):=NULL
]

peakAnno_genExpr <- copy(peakAnno_strongest_peak)[genExpr,on='EnsemblID',nomatch=0]

## test enrichment DA peaks near DE genes
subset_peakAnno_genExpr = copy(peakAnno_genExpr)[de_logFC>=1 & da_logFC>=1 ]
da_nonde_genes = length(unique(copy(subset_peakAnno_genExpr)[DA=='da'& DE=='non_de']$EnsemblID))
da_de_genes = length(unique(copy(subset_peakAnno_genExpr)[DA=='da'& DE=='de']$EnsemblID))
nonda_de_genes = length(unique(copy(subset_peakAnno_genExpr)[DA=='non_da'& DE=='de']$EnsemblID))
nonda_nonde_genes = length(unique(copy(subset_peakAnno_genExpr)[DA=='non_da'& DE=='non_de']$EnsemblID))

matrix <- data.frame(
  as.numeric(da_de_genes),
  as.numeric(nonda_de_genes),
  as.numeric(da_nonde_genes), 
  as.numeric(nonda_nonde_genes)
)%>%unlist()%>%matrix(nrow=2,byrow=T)

fisher.test(matrix)
# data:  matrix
# p-value = 0.6813
# alternative hypothesis: true odds ratio is greater than 1
# 95 percent confidence interval:
#  0.8962623       Inf
# sample estimates:
# odds ratio 
#   0.977224 

logfc_DaDe <- copy(peakAnno_genExpr)[,da_logFC:=logFC][,divergent:=as.factor(ifelse(abs(de_logFC)>=1 & abs(da_logFC)>=1,'y','n'))]

count_genes = function(x){
  maxlogfc = 10

  df = copy(x)[divergent=='y'][
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

numb_genes <- count_genes(logfc_DaDe)
numb_genes

pdf(paste(outplot_dir,'logfc-da-de.pdf',sep=''),width= 7,height = 7)
ggplot(logfc_DaDe, aes(x=de_logFC, y=da_logFC)) +
geom_point(alpha=0.2,aes(color=divergent))+
scale_color_manual(values=c('darkgray','red'))+
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
