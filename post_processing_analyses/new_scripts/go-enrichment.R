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
library(biomaRt)
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

## GET PEAKS
da_results <- read_da_results('da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

random_peaks <- fread('./output/files/random_regions/random_genomic_regions.txt',sep='\t',header=T)
random_peaks <- random_peaks[,peakID:=paste('peak_',1:nrow(random_peaks),sep='')]

## LOOK FOR PEAK TARGET GENES USING CHIPSEEKER 
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db) 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tss_region = 10000

annotate <- function(peak){
  peakAnno <- annotatePeak(
    makeGRangesFromDataFrame(copy(peak),keep.extra.columns=T), 
    tssRegion=c(-tss_region, tss_region),
    TxDb = txdb,
    sameStrand=T,
    level = "gene",
    overlap='all'
  )%>%as.data.table()

  symbols <- mapIds(org.Hs.eg.db, keys = peakAnno$geneId, keytype = "ENTREZID", column="SYMBOL")%>%na.omit()
  symbols_df <- data.table(symbol=symbols,entrezid=names(symbols))%>%unique()
  ensemblid <- mapIds(org.Hs.eg.db, keys = symbols_df$symbol, keytype = "SYMBOL", column="ENSEMBL")%>%na.omit()
  ensemblid_df <- data.table(symbol=names(ensemblid),EnsemblID=ensemblid)

  gene_ids <- copy(ensemblid_df)[symbols_df,on='symbol',nomatch=0]%>%setnames(old='entrezid',new='geneId')
  peakAnno <- copy(peakAnno)[gene_ids,on='geneId',nomatch=0][,c('width','strand','geneLength','geneStrand','geneId','annotation'):=NULL]
  return(peakAnno)
}

allPeaks <- list(da_results,random_peaks[,DA:='random'])
names(allPeaks) = c('mypeaks','randomPeaks')

allPeaks <- lapply(allPeaks,function(x)annotate(x))

## proportion of peaks that have been associated with genes
x = length(unique(allPeaks[[1]]$peakID))
y = length(unique(da_results$peakID))
x/y*100
# 94.75016

## look @ the overlap between genes targeted by DA and non-DA peaks  
genes <- copy(allPeaks[[1]])%>%setorderv('DA')%>%split(by='DA')%>%lapply(function(x)x=x[,EnsemblID]%>%unique()%>%na.omit())

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
    col = da_palette, 
    fill = c(alpha(da_palette[[1]],0.3),alpha(da_palette[[2]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = 'black'
)

## LOOK @ ORTHOLOGY INFORMATION FOR TARGET GENES
## get orthologous genes for each species and retain only genes with orthologous information
orthologous_attributes = paste(
        'ptroglodytes',c(
            'homolog_ensembl_gene','homolog_associated_gene_name',
            'homolog_orthology_type','homolog_perc_id'
            ),
            sep='_'
)

mart <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')

humanChimpOrthologs <-  getBM(
        c('external_gene_name','ensembl_gene_id',orthologous_attributes),
        mart = mart,uniqueRows=T
)%>%as.data.table()
humanChimpOrthologs <- humanChimpOrthologs%>%na.omit()

genes <- lapply(genes,function(x)dt = data.table(ensembl_gene_id=x))
genes <- Map(mutate,genes,da=names(genes))%>%rbindlist()

humanChimpOrthologs <- humanChimpOrthologs[genes,on='ensembl_gene_id',nomatch=0]

## Plot orthology type and homology % for genes targeted by DA and nonDA peaks
orthologyType <- copy(humanChimpOrthologs)[
  ,numb_orthology_type :=.N,by=.(da,ptroglodytes_homolog_orthology_type)
  ][
    ,numb_da := .N,by=.(da) 
    ][
      ,prop := numb_orthology_type/numb_da
]

## plot the % of one2one orthologous genes between anchor and non-anchor sets of genes
pdf(paste(outplot_dir,'distrib-homology-percent.pdf',sep='/'),width=7,height=7)
dt <- copy(orthologyType)[,c('da','ptroglodytes_homolog_perc_id','external_gene_name')]%>%unique()
ggplot(dt,aes(x=da,y=ptroglodytes_homolog_perc_id,fill=da))+
geom_violin(scale = "width")+xlab(' ')+ylab('Homology percentage')+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
scale_fill_manual(values=da_palette)+
stat_compare_means(
  method = "wilcox.test",
  label.y = max(dt$ptroglodytes_homolog_perc_id+0.1),
  size=5)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

pdf(paste(outplot_dir,'distr-homology-type.pdf',sep='/'),width=7,height=7)
dt <- copy(orthologyType)[,c('da','ptroglodytes_homolog_orthology_type','prop')]%>%unique()
ggplot(dt,aes(x=da, y=prop,fill=ptroglodytes_homolog_orthology_type))+
geom_bar(stat = 'identity',position='stack')+xlab(' ')+ylab('Proportion of genes')+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

## RUN GO ENRICHMENT ANALYSIS USING CLUSTERPROFILER
genesSplit <- copy(genes)%>%split(by='da')%>%lapply(function(x)x$ensembl_gene_id)
compGOs <- compareCluster(
  geneCluster = genesSplit,
  ont = "BP",
  fun = "enrichGO",
  keyType = "ENSEMBL",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "fdr",
  OrgDb='org.Hs.eg.db'
)

pdf(paste(outplot_dir,'go_enrichment.pdf',sep=''),width=10,height = 15)
dotplot(compGOs, showCategory = 15, title = "GO Enrichment Analysis")
dev.off()

## assign peaks to genes based on strongest peak (logfc)
# peakAnno_strongest_peak <- copy(peakAnno)[ ,.SD[which.max(abs(logFC))], by=.(EnsemblID)][,bindist:=plyr::round_any(distanceToTSS, 1000)]

# regulation <- rbindlist(copy(peakAnno))[,c('DA','EnsemblID','da_species')]%>%unique()
# regulation = regulation[,npeaktype:=.N,by=.(EnsemblID)][,regulation:=ifelse(npeaktype==2,'da',DA)][,remove:=ifelse( (DA=='non_da'& regulation=='da'),'yes','no')][remove=='no'][,remove:=NULL]

# peakAnno_closest_peak <- rbindlist(copy(peakAnno))[regulation, on=c('EnsemblID','DA','da_species'),nomatch=0][,npeaktype:=NULL][,.SD[which.min(abs(distanceToTSS))], by=.(EnsemblID)]

## distribution peak-tss distances
pdf(paste0(outplot_dir,'qc-peak-tss-distance.pdf',sep=''),width = 10, height = 8)
df = copy(allPeaks[[1]])[,c('EnsemblID','peakID','DA','distanceToTSS')][abs(distanceToTSS)<=1000000]%>%unique()
ggplot(df,aes(x=distanceToTSS,fill=DA))+
geom_density(alpha=0.5)+
ylab('Density')+xlab(' ')+
scale_fill_manual(values=da_palette)+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

# nrow(peakAnno_strongest_peak[da_species=='chimp'])
# nrow(peakAnno_strongest_peak[da_species=='human'])
# nrow(peakAnno_strongest_peak[da_species=='non_da'])


##---------------------------------------------
## loeuf score to check whether genes are
## less loss-of-function tolerant 
##---------------------------------------------
gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('chromosome','start_position','end_position','gene_id','oe_lof_upper'))%>%setnames(old=c(1:4),new=c(range_keys,'EnsemblID'))

loeuf_score_mygenes <- allPeaks[[1]][gnomad_loeuf,on='EnsemblID',nomatch=0][,type:=DA][,c('type','EnsemblID','oe_lof_upper')]
loeuf_score_all_genes <- copy(gnomad_loeuf)[,type:='all_genes'][,c('type','EnsemblID','oe_lof_upper')]

all_loeuf_scores <- rbind(loeuf_score_mygenes[[1]],loeuf_score_all_genes)%>%na.omit()

comparisons = list(
  c('da','all_genes'),
  c('da','non_da'),
  c('non_da','all_genes')
)

pdf(paste0(outplot_dir,'loeuf-score-target-genes.pdf',sep=''),width = 10, height = 8)
loeuf_palette = c('#edf2f4',da_palette)
names(loeuf_palette) = c('all_genes',names(da_palette))

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
write.table(peakAnno, file = paste(outfile_dir,'peak-target-genes.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## GO tables
go_table <- copy(compGOs@compareClusterResult)%>%as.data.table()
go_table <- go_table[,geneID:=NULL]
write.xlsx(go_table,paste(outfile_dir,'go-enrich-terms-da-vs-nonda.xlsx',sep=''),append=T,overwrite=T)
