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

# table_dir <- './post_processing_analyses/output/tables/GO_enrich/'
# target_genes_dir <- './post_processing_analyses/output/GO/target_genes/'

## get DA peaks
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','peak_species','logFC')]

## use chipseeker to annotate peaks within 1Mb of distance to their putative target gene(s)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

tss_region = 1000000

peakAnno <- annotatePeak(
  makeGRangesFromDataFrame(da_results,keep.extra.columns=T), 
  tssRegion=c(-tss_region, tss_region),
  TxDb=txdb, 
  annoDb="org.Hs.eg.db"
)%>%as.data.table()%>%setnames(old='ENSEMBL',new='EnsemblID')

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
#        DA     n tot_n      prop
# 1: non_da 52851 55731 0.9483232
# 2:     da 16592 17701 0.9373482

## number peaks per gene
pdf(paste(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width=10,height = 7)
df <- copy(peakAnno)[,c('EnsemblID','DA','peakID')]%>%na.omit()%>%unique()
df <- df[,numb_peaks_per_gene:=.N,by=.(DA,EnsemblID)]
ggplot(df,aes(x=numb_peaks_per_gene,fill=DA))+
geom_bar(position=position_dodge())+
scale_fill_manual(values=da_palette)+
facet_wrap(DA~.,ncol=2)+
theme_classic()+
theme(
  legend.key = element_rect(fill = "white", colour = "black")
)
dev.off()

# mean(df[DA=='da']$numb_peaks_per_gene)
# [1] 5.25
# mean(df[DA!='da']$numb_peaks_per_gene)
# [1] 8.15

## assign closest peak to gene
# one_peak_per_gene <- copy(peakAnno)[,distTSS:=start-geneStart][,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]

peak_regulation <- copy(peakAnno)[,c('EnsemblID','DA')]%>%unique()
peak_regulation <- peak_regulation[,reg:=.N,by=.(EnsemblID)][,peaktype:=ifelse(reg==1,DA,'da')]

one_peak_per_gene <- copy(peakAnno)[,distTSS:=start-geneStart][,.SD[which.min(abs(distTSS))], by=.(EnsemblID)]
# [peak_regulation,on=c('EnsemblID','DA'),nomatch=0][,DA:=NULL][,DA:=peaktype][,peaktype:=NULL]

pdf(paste(outplot_dir,'dist_peak_tss.pdf',sep=''),width=10,height = 7)
df <- copy(one_peak_per_gene)[,c('distTSS','DA','peakID')]%>%na.omit()%>%unique()
df <- df[,log10dist:=log10(abs(distTSS)+0.1)]
ggplot(df,aes(x=log10dist,fill=DA))+
geom_density(alpha=0.5)+
scale_fill_manual(values=da_palette)+
theme_classic()+
theme(
  legend.key = element_rect(fill = "white", colour = "black")
)
dev.off()

# ## rank genes by enrichment of DA peaks
# df <- copy(peakAnno)[,c('ENSEMBL','DA','peakID')]%>%na.omit()%>%unique()
# df <- df[,numb_peaks_per_gene:=.N,by=.(DA,ENSEMBL)][,peakID:=NULL]
# df_da <- copy(df[DA=='da'])
# df_nonda <- copy(df[DA!='da'])

# observed_value = as.data.table(table(copy(df_da)$ENSEMBL))%>%setnames(old=c(1:2),new=c('ensembl','score'))%>%split(by='ensembl')

# list_df_permuted_scores = list()
# for(i in 1:10000){
#       permuted_scores=table(sample(unique(df$ENSEMBL), nrow(df), replace = T))
#       list_df_permuted_scores[[i]] = as.data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('ensembl','permuted_score'))
# }
# list_df_permuted_scores = rbindlist(list_df_permuted_scores)[ensembl%in%names(observed_value)]%>%split(by='ensembl')

# enrichment = purrr::map2(list_df_permuted_scores,observed_value,function(p,o){
#   zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
#   pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#   stat =  data.table(zscore = zscore,pval=pvalue)
#   return(stat)
#   }
# )
# enrichment <- Map(mutate,enrichment,ensembl=names(enrichment))%>%rbindlist()%>%setorderv('zscore',1)
# # enrichment <- enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]

##---------------------------------------------
## GO enrichment using clusterProfiler
##---------------------------------------------
library(clusterProfiler)

gene_regulation <- copy(one_peak_per_gene)[,c('DA','EnsemblID')]%>%unique()%>%split(by='DA')%>%lapply(function(x)x$EnsemblID)

lapply(gene_regulation,function(x)length(x))
# $da
# [1] 3900
# $non_da
# [1] 14435
length(setdiff(gene_regulation[[1]],gene_regulation[[2]]))
# [1] 14435

compGOs <- compareCluster(
  geneCluster = gene_regulation,
  ont = "BP",
  fun = "enrichGO",
  keyType = "ENSEMBL",
  pvalueCutoff  = 0.05,
  # universe = c(gene_regulation[[1]],gene_regulation[[2]]),
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
gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
colnames(gnomad_loeuf)[1] = 'EnsemblID'

loeuf_score_mygenes <- copy(one_peak_per_gene)[gnomad_loeuf,on='EnsemblID',nomatch=0][,type:=DA][,c('type','EnsemblID','oe_lof_upper')]
loeuf_score_all_genes <- copy(gnomad_loeuf)[,type:='all_genes'][,c('type','EnsemblID','oe_lof_upper')]

all_loeuf_scores <- rbind(loeuf_score_mygenes,loeuf_score_all_genes)%>%na.omit()

comparisons = list(
  c('da','all_genes'),
  c('non_da','all_genes'),
  c('non_da','da')
)

pdf(paste0(outplot_dir,'loeuf_score_target_genes.pdf',sep=''),width = 10, height = 8)
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
target_genes <- copy(one_peak_per_gene)[,c(..range_keys,'peakID','DA','EnsemblID','SYMBOL','geneId','geneStart')]%>%unique()
target_genes <- target_genes[,distTSS:=start-geneStart][,geneStart:=NULL]

write.table(target_genes, file = paste(outfile_dir,'target_genes.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)

## GO tables
go_table <- copy(compGOs@compareClusterResult)%>%as.data.table()
go_table <- go_table[,geneID:=NULL]
write.xlsx(go_table,paste(outfile_dir,'GO_enrich_terms_da_non_da.xlsx',sep=''),append=T,overwrite=T)


# ## GREAT enrichments
# ## generate random background of 100000 regions located within 1Mb of distaance from 
# ## hg38 genes
# # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# # txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# # cres <- promoters(txdb, upstream=10000, downstream=10000)%>%as.data.table()
# # cres <- cres[!seqnames%like% "M|Un|_"]

# # # random_cres <- sample_n(cres,10000)%>%as.data.table()
# # random_cres <- copy(cres)
# # random_cres <- random_cres[,peakID:=paste('bgpeak_',1:nrow(random_cres),sep='')][,c(..range_keys,'peakID')]

# test_peaks <- split(copy(da_results),by='DA')
# background <- copy(da_results)

# get_enrichment = function(gr){
#     enrich = submitGreatJob(
#     gr=gr,
#     bg=background,
#     species = genome,
#     rule= "twoClosest", 
#     adv_twoDistance = 1000,  
#     includeCuratedRegDoms = T,
#     max_tries=1,
#     request_interval=10
#     )
    
#     go_enrich = getEnrichmentTables(enrich,ontology='GO Biological Process')%>%rbindlist()
#     go_enrich = go_enrich[Hyper_Adjp_BH <= 0.05]
    
#     # get list of putatve target genes 
#     targets = plotRegionGeneAssociationGraphs(enrich,type=1,request_interval = 10)%>%as.data.table() %>% na.omit() %>% unique()
#     great_enrich = list(go_enrich,targets)
#     return(great_enrich)
# }

# go_enrichment <- lapply(test_peaks,function(x)get_enrichment(x))

# go_tables <- list(go_enrichment[[1]][[1]],go_enrichment[[2]][[1]])
# names(go_tables) = names(go_enrichment)

# ## plot GO results 
# go_bp <- Map(mutate,go_tables,DA = as.factor(names(go_tables)))
# go_bp <- lapply(go_bp,function(x)
#  x=x[
#    1:15,c(1,2,13,14)
#    ][
#      ,log10_adj_p := -log10(Hyper_Adjp_BH)
#      ]
# )%>%rbindlist()%>%na.omit()%>%setorderv('log10_adj_p',1)%>%setorderv('DA',1)

# go_enrichment_plot = function(df){
#   p <- ggplot(df, aes(x=factor(df$name,levels=df$name),y=log10_adj_p,fill=DA)) +
#   geom_bar(stat = 'identity',position = 'dodge',col='black')+
#   scale_fill_manual(values=da_palette)+
#   xlab(" ") +ylab("\n -Log10 (P adj.) \n ") +
#   theme(
#     legend.position='bottom',
#     legend.background=element_rect(),
#     axis.text.x=element_text(angle=0, hjust=1.10),
#     axis.text.y=element_text(angle=0, vjust=0.8),
#     axis.title=element_text(),
#     axis.line = element_line(color = "black",size = 0, linetype = "solid"),
#     panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_blank(),
#     title=element_text(),
#     strip.background = element_blank()
#     ) +
#   coord_flip()
#   return(p)
  
# }

# pdf(paste(outplot_dir,'go_barplot.pdf',sep=''),width=10,height = 15)
# go_enrichment_plot(go_bp)
# dev.off()

# ## check ontology index of significant GO terms 
# library(ontologyIndex)
# data(go)

# significant_gos <- copy(go_tables)%>%lapply(function(x)x[,'ID']%>%unique()%>%split(by='ID'))

# ancestor_terms <- copy(significant_gos)%>%lapply(
#   function(x) x<-lapply(x,function(y)
#   y<-length(get_ancestors(go,as.character(y))))
# )

# ancestor_terms <- lapply(ancestor_terms,function(x){
#   goids <- names(x)
#   table <- lapply(x,function(y)y=data.table(ancestor_terms=y))%>%rbindlist()
#   table <- table[,goids:=goids]
#   return(table)
# })

# ancestor_terms <- Map(mutate,ancestor_terms,DA=names(ancestor_terms))%>%rbindlist()

# pdf(paste0(outplot_dir,'ancestor_terms.pdf',sep=''),width = 10, height = 8)
# ggplot(ancestor_terms,aes(x=DA,y=ancestor_terms,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# scale_fill_manual(values=da_palette)+
# geom_hline(yintercept=median(ancestor_terms[DA=='da']$ancestor_terms),linetype='dashed')+
# stat_compare_means(
#   method = "wilcox.test",
#   ref.group = 'da',
#   size=5
#   )+xlab(' ')+ylab('Number of ancestor terms per significant GO \n associated with target genes')+
# theme_bw()+
# theme(
#   legend.position = 'none',
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
# )
# dev.off()

# ##-------------------------
# ## check the target genes
# ##-------------------------
# target_genes <- list(go_enrichment[[1]][[2]],go_enrichment[[2]][[2]])
# names(target_genes) = names(go_enrichment)

# ## get number and % peaks with target genes
# numb_peaks_w_gene <- copy(target_genes)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
# numb_peaks <- copy(test_peaks)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
# purrr::map2(numb_peaks_w_gene,numb_peaks,`/`)
# numb_genes <- copy(target_genes)%>%lapply(function(x)x=x[,'gene']%>%unique()%>%nrow())

# target_genes <- Map(mutate,target_genes,DA=names(target_genes))%>%rbindlist()

# ## plot venn diagramm of shared genes between peaks
# library(VennDiagram)
 
# genes <- copy(target_genes)%>%split(by='DA')%>%lapply(function(x)x=x[,gene]%>%unique())
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

# gene_regulation <- copy(target_genes)[,c('gene','DA')]%>%unique()
# gene_regulation <- gene_regulation[
#   ,peaktype_assoc:=.N,by=.(gene)
#   ][
#     ,regulation:=ifelse(peaktype_assoc>1,'nonda_da',DA)
#     ][
#       ,peaktype_assoc:=NULL
# ]

# target_genes <- target_genes[
#   gene_regulation,on=c('DA','gene'),nomatch=0
#   ][
#     ,log10_abs_dist:= log10(abs(distTSS)+1)
#     ][
#       ,c('width','strand'):=NULL
# ]

# ## add ensembl ID info
# # library(biomaRt)
# # ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl") 

# # genes_assoc_w_signf_go <- lapply(
# #   go_tables,function(x){ 
# #     go_id_genes <-copy(x)[,'ID']%>%split(by='ID')%>%
# #     lapply(
# #       function(y)y<-getBM(
# #         attributes = c('hgnc_symbol'), 
# #         filters = 'go', 
# #         values = x$ID, 
# #         mart = ensembl
# #         )
# #         )%>%rbindlist()
# #   }
# # )

# ## get loeuf score genes targeted by da, non-da and da + non-da peaks 
# gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene','gene_id','oe_lof_upper'))
# colnames(gnomad_loeuf)[2] = 'EnsemblID'

# target_genes_w_loeuf <- target_genes[gnomad_loeuf,on='gene',nomatch=0]%>%unique()%>%na.omit()
  
# ## divide genes by peak 
# target_genes_w_loeuf <- target_genes_w_loeuf[
#   ,c('gene','oe_lof_upper','regulation')
# ]%>%unique()


# all_genes_loeuf <- copy(gnomad_loeuf)[,EnsemblID:=NULL][,regulation:='all_genes']%>%na.omit()

# comparisons = list(
#   c('da','all_genes'),
#   c('non_da','all_genes'),
#   c('nonda_da','all_genes')
# )

# pdf(paste0(outplot_dir,'loeuf_score_target_genes.pdf',sep=''),width = 10, height = 8)
# df <- rbind(all_genes_loeuf,target_genes_w_loeuf%>%dplyr::select(c(all_of(colnames(all_genes_loeuf)))))
# ggplot(df,aes(x=regulation,y=oe_lof_upper,fill=regulation))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# geom_hline(yintercept=median(all_genes_loeuf$oe_lof_upper),linetype='dashed',size=0.5)+
# ylab('LOEUF score')+xlab('Type of regulation')+
# stat_compare_means(
#   method = "wilcox.test",
#   comparisons = comparisons,
#   size=5
#   )+
#   theme_bw()+
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
# )
# dev.off()

# ##  plot distance between peak and target genes
# target_genes <- target_genes[,log10_abs_dist:= log10(abs(distTSS)+1)]

# # pdf(paste(outplot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
# # df <- copy(target_genes)[,log10_dist:=ifelse(distTSS<0,-log10_abs_dist,log10_abs_dist)]
# # ggplot(df, aes(x=log10_dist,fill=DA)) +
# #     geom_density(alpha=0.5)+scale_fill_manual(values=da_palette)
# # dev.off()

# pdf(paste(outplot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
# ggplot(target_genes,aes(x=DA,y=log10_abs_dist,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# geom_hline(yintercept=median(target_genes[DA=='da']$log10_abs_dist),linetype='dashed',size=0.5)+
# xlab('')+ylab('Log10 bp dist from TSS')+
# scale_fill_manual(values=da_palette)+
# stat_compare_means(
#   method = "wilcox.test",
#   ref.group = "da",
#   label.y = (max(target_genes$log10_abs_dist)+0.5),
#   size=5
#   )+
# theme_bw()+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()

# # ## look at the distribution across genomic elements of DA vs non-DA peaks associated with genes 
# # library(annotatr)

# # annots <- c(
# #   'hg38_basicgenes',
# #   'hg38_genes_intergenic',
# #   'hg38_genes_intronexonboundaries'
# # )
# # annotations <- build_annotations(genome = 'hg38', annotations = annots)

# # annotated_peaks <- annotate_regions(
# #     regions = makeGRangesFromDataFrame(target_genes,keep.extra.columns=T),
# #     annotations = annotations,
# #     ignore.strand = TRUE,
# #     quiet = FALSE
# # )%>%as.data.table()

# # annotated_peaks <- annotated_peaks[,c(..range_keys,"annot.id","DA","regulation")][,annotation_id:=gsub('\\:.*','',annot.id)]%>%unique()
# # test <- copy(annotated_peaks)[
# #   ,numbpeaks:=.N,by=.(DA)
# #   ][
# #     ,numbpeaks_annot:=.N,by=.(DA,annotation_id)
# #     ][
# #       ,proppeaks_annot:=numbpeaks_annot/numbpeaks
# #       ][
# #         ,c('DA','annotation_id','proppeaks_annot')
# # ]%>%unique()

# # pdf(paste(outplot_dir,'prop_peak_target_genes_annotation.pdf',sep=''),width=10,height = 7)
# # ggplot(test,aes(x=reorder(annotation_id, -proppeaks_annot),,y=proppeaks_annot,fill=DA))+
# # geom_bar(stat='identity',position=position_dodge())+
# # xlab('')+ylab('Proportion peaks')+
# # scale_fill_manual(values=da_palette)+
# # theme_bw()+
# # theme(
# #   legend.key = element_rect(fill = "white", colour = "black"),
# #   axis.line = element_blank()
# # )
# # dev.off()

# ##---------------
# ## export files
# ##---------------
# ## target genes
# write.table(target_genes, file = paste(outfile_dir,'target_genes.txt',sep=''),col.names = T, row.names = F, sep = "\t", quote = F)
# ## GO tables
# write.xlsx(go_tables,paste(outfile_dir,'GO_enrich_terms_da_non_da.xlsx',sep=''),append=T,overwrite=T)