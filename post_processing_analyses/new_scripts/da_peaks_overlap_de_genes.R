## script used to look for overlap between DA peaks and DE genes
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(VennDiagram)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'rna_seq')
target_genes_dir <- './output/files/GO_enrich'
tmp_files_dir <- './output/temp_files/'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'logFC','DA','peakID','peak_species','FDR'))
setkeyv(da_results,range_keys)
da_results <- da_results%>%setorderv('DA',1)%>%split(by='DA')

## read target genes (need to add logFC info to these)
target_genes <- list.files(target_genes_dir,full.names= T,recursive=F,pattern='target')%>%lapply(
  function(x)fread(x,sep='\t',header=T,select=c('peakID','distTSS','gene','go_signif','DA')))%>%rbindlist()%>%split(by='DA')
  
peaks_w_genes <- purrr::map2(da_results,target_genes,function(x,y)x[y,on=c('peakID','DA'),nomatch=0])%>%rbindlist()
##-----------------
## read DE genes 
##-----------------
## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
de_genes <- fread(
  "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom.sex.interact_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,de_logFC:=-logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','t','B','genes','logFC'):=NULL
]

##-----------------------------------------
## keep genes with ensembl ids from v 86
##-----------------------------------------
# ensembl_86_ids <- fread("../rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

library(biomaRt)

ensembl <- useMart('ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl") 

ensembl_ids <- getBM(
        attributes = c('hgnc_symbol','ensembl_gene_id'), 
        filters = 'hgnc_symbol', 
        values = unique(peaks_w_genes$gene), 
        mart = ensembl
)%>%as.data.table()
colnames(ensembl_ids) = c('gene','EnsemblID')

peaks_w_genes_ensembl <- copy(peaks_w_genes)[ensembl_ids,on='gene',nomatch=0][de_genes,on='EnsemblID',nomatch=0]%>%unique()%>%setorderv('DA')

## count numb and prop of peaks and genes retained 
numb_peaks_genes <- copy(peaks_w_genes)[,c('DA','peakID')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))
numb_peaks_ensembl <- copy(peaks_w_genes_ensembl)[,c('DA','peakID')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))

purrr::map2(numb_peaks_ensembl,numb_peaks_genes,`/`)

numb_genes <- copy(peaks_w_genes)[,c('DA','gene')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))
numb_genes_ensembl <- copy(peaks_w_genes_ensembl)[,c('DA','gene')]%>%unique()%>%split(by='DA')%>%lapply(function(x)nrow(x))

purrr::map2(numb_genes_ensembl,numb_genes,`/`)

## check how many peaks each gene has associated and stratify that by de/non-de  
peaks_w_genes_ensembl <- peaks_w_genes_ensembl[,numb_peaks_per_gene:=.N,by=.(gene,DE)]

pdf(paste0(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width = 10, height = 8)
ggplot(peaks_w_genes_ensembl,aes(x=DE,y=numb_peaks_per_gene,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
geom_hline(yintercept=median(peaks_w_genes_ensembl[DE=='de']$numb_peaks_per_gene),linetype='dashed')+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

## compare proportion da peaks per de/non de gene
peaks_w_genes_ensembl <- peaks_w_genes_ensembl[
  ,numb_peak_type_per_gene:=.N,by=.(gene,DE,DA)
  ][
    ,prop_peaktype_per_gene:=numb_peak_type_per_gene/numb_peaks_per_gene
]

pdf(paste0(outplot_dir,'prop_da_peaks_per_gene.pdf',sep=''),width = 10, height = 8)
df <- copy(peaks_w_genes_ensembl)[DA=='da'][,c('prop_peaktype_per_gene','gene','DA','DE')]%>%unique()
ggplot(df,aes(x=DE,y=prop_peaktype_per_gene,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
geom_hline(yintercept=median(df[DE=='de']$prop_peaktype_per_gene),linetype='dashed')+
ylab('proportion DA peaks per gene')+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

## plot avg logFC between peaks associated same gene
pdf(paste0(outplot_dir,'avg_logfc_peaks_per_gene.pdf',sep=''),width = 10, height = 8)
df <- copy(peaks_w_genes_ensembl)[
  ,.SD[which.max(abs(logFC))], by=.(gene)
]
ggplot(df,aes(x=DE,y=logFC,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
geom_hline(yintercept=median(df[DE=='de']$logFC),linetype='dashed')+
ylab('max ATAC logFC peaks per gene')+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

## keep peaks with genes in same tad
orth_tads <- read_tads('human_chimp_orth_tads.txt')[,c(1:3)]
colnames(orth_tads) = range_keys
setkeyv(orth_tads,range_keys)

peaks_w_ensembl_in_tad<- foverlaps(copy(peaks_w_genes_ensembl),orth_tads,type='within')%>%na.omit()%>%dplyr::select(-c(contains('per_gene')))
peaks_w_ensembl_in_tad <- peaks_w_ensembl_in_tad[
  ,same_tad:=ifelse(distTSS<0 & start- (i.start+distTSS)>0,'y',
  ifelse(distTSS>0 & end - (i.end+distTSS)>0,'y','n')
  )
  ][
    same_tad=='y'
    ][
      ,c(range_keys[-1]):=NULL
]%>%unique()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

## prop peaks & genes retained
100*(length(unique(peaks_w_ensembl_in_tad$peakID))/length(unique(peaks_w_genes_ensembl$peakID)))
100*(length(unique(peaks_w_ensembl_in_tad$EnsemblID))/length(unique(peaks_w_genes_ensembl$EnsemblID)))

one_peak_per_gene <- copy(peaks_w_ensembl_in_tad)
one_peak_per_gene<- one_peak_per_gene%>%split(by='gene')%>%lapply(
  function(x){
    x<-x[
      ,da_length:=length(levels(as.factor(x$DA)))
      ][
        !(da_length==2 & DA=='non_da')
      ]
  }
)%>%rbindlist()

one_peak_per_gene <- one_peak_per_gene[
        ,.SD[which.min(abs(distTSS))], by=.(gene)
        ][
        ,.SD[which.min(abs(distTSS))], by=.(peakID)
]%>%dplyr::select(-c('same_tad','da_length'))

##-----------------------------------------
## enrichment de genes with da peaks
##-----------------------------------------
significant_da_genes <- copy(one_peak_per_gene)[DA=='da']$EnsemblID%>%unique()
significant_de_genes <- copy(one_peak_per_gene)[DE=='de']$EnsemblID%>%unique()

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
   col=c('#fde725ff',"#440154ff"),
    fill = c(alpha('#fde725ff',0.3),alpha("#440154ff",0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = c('#fde725ff',"#440154ff")
)

## test significance overlap
da_de_genes <- Reduce(intersect,da_de_overlap)%>%unique()
nonda_de_genes <- copy(one_peak_per_gene)[DA!='da'][DE=='de'][,EnsemblID]%>%unique()

da_nonde_genes <- copy(one_peak_per_gene)[DA=='da'][DE!='de'][,EnsemblID]%>%unique()
nonda_nonde_genes <- copy(one_peak_per_gene)[DA!='da'][DE!='de'][,EnsemblID]%>%unique()

matrix <- data.frame(
  as.numeric(length(da_de_genes)),
  as.numeric(length(nonda_de_genes)), 
  as.numeric(length(da_nonde_genes)),
  as.numeric(length(nonda_nonde_genes))
)%>%unlist()%>%matrix(nrow=2,byrow=T)

enrichment_dade_overlap <-fisher.test(matrix) # 0.01792 

## correlation heatmap between DA logFC and DE logFC 
library(ComplexHeatmap)
corr_logFCs <- copy(one_peak_per_gene)[DA=='da' & DE=='de'][,c('de_logFC','logFC','gene','distTSS')]%>%unique()
corr_logFCs <- corr_logFCs[
    ,.SD[which.min(abs(distTSS))], by=.(gene)
]
pdf(paste0(outplot_dir,'heatmap_logfold_changes.pdf',sep=''),width = 10, height = 8)
Heatmap(cor(corr_logFCs[,-c('gene','distTSS')]),name = "Pearson corr",col=viridis::viridis(100)) 
dev.off()

## assign to each gene the peak logFC, for genes with > 1 peak associated calculate the median logFC
# one_peak_per_gene <- copy(peaks_w_genes_ensembl)%>%
# dplyr::select(-c(contains('per_gene')))%>%
# split(by='EnsemblID')%>%lapply(function(x)x<-x[,da_logFC:=median(logFC)][,c('logFC','distTSS'):=NULL])%>%rbindlist()

# one_peak_per_gene <- copy(peaks_w_genes_ensembl)[
#   ,.SD[which.max(abs(logFC))], by=.(EnsemblID)
# ]%>%dplyr::select(-c(contains('per_gene')))


# one_peak_per_gene <- copy(peaks_w_ensembl_in_tad)
# one_peak_per_gene<- one_peak_per_gene[
#   ,da_length:=length(levels(as.factor(one_peak_per_gene$DA))),by=.(gene)
#   ][
#       !(da_length==2 & DA=='non_da') ## if a gene has da and non da peaks associated, keep the das
#       ][
#         ,.SD[which.min(abs(distTSS))], by=.(peakID)
# ]%>%dplyr::select(-c('same_tad','da_length'))

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


# logfc_da_de_genes <- copy(one_peak_per_gene)[,da_logFC:=logFC][DA=='da'][DE=='de']

# count_genes_per_quadrant = function(x){
#   df = copy(x)[
#     ,quadrant:=ifelse((da_logFC>0 & de_logFC<0),'1',
#       ifelse((da_logFC>0 & de_logFC>0),'2',
#         ifelse((da_logFC<0 & de_logFC<0),'3','4'
#         )))
#         ][
#           ,maxDAfc:=max(da_logFC)
#           ][
#             ,minDAfc:=min(da_logFC)
#             ][
#               ,maxDEfc:=max(de_logFC)
#               ][
#                 ,minDEfc:=min(de_logFC)
#                 ][
#                   ,c("EnsemblID",'quadrant','maxDAfc','minDAfc','maxDEfc','minDEfc')
#                   ]%>%unique()
#   df = df[
#     ,numbgenes:=.N,by=.(quadrant)
#     ][
#       ,EnsemblID:=NULL
#       ]%>%unique()%>%setorderv('quadrant',1)

#   df = df[
#     ,y:=ifelse((quadrant=='1'|quadrant=='2'),maxDAfc+0.5,minDAfc-0.5)
#     ][
#       ,x:=ifelse((quadrant=='1'|quadrant=='3'),minDEfc+0.5,maxDEfc-0.5)
#       ]
#   return(df)
# }

# logfc_plot = function(x){
#   df=copy(x)
#   numb_genes = count_genes_per_quadrant(x)
#   p = ggplot(df, aes(x=de_logFC, y=da_logFC)) +
#     geom_point()+
#     geom_text(data = numb_genes, aes(label = numb_genes$numbgenes,y = numb_genes$y, x = numb_genes$x))+
#     geom_vline(xintercept=0,linetype='dashed')+
#     geom_hline(yintercept=0,linetype='dashed')+
#     ylab('DA logFC')+xlab('DE logFC')+
#     theme()
  
#   return(p)
# }

# pdf(paste(outplot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 7,height = 7)
# logfc_plot(logfc_da_de_genes)
# dev.off()

##---------------------------------------------------
## look at patterns accessibility around DE genes
##---------------------------------------------------
library(edgeR)
library(csaw)

orthologous_peak_file <- paste('output/final_peak_set/',genome,'_all_orthologous_peaks.txt',sep='')

## list bam
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(pe = "both",restrict=standard_chr,max.frag=1000)

get_bams <- function(species){
    bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}

bams <- get_bams(genome)

## read consensus peak
orthologous_peaks <- fread(orthologous_peak_file,sep='\t',header=T)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

## count reads in peaks (rpkm)
reads_in_peaks <- regionCounts(bams, orthologous_peaks,param=param)

## remove peaks with low counts
abundances <- aveLogCPM(asDGEList(reads_in_peaks))
keep <- abundances > 0

reads_in_peaks_filtered <- reads_in_peaks[keep,]

reads_in_peaks_dge <- asDGEList(reads_in_peaks_filtered)
rownames(reads_in_peaks_dge) <- rowData(reads_in_peaks_filtered)$peakID
colnames(reads_in_peaks_dge) <- samples_names
reads_in_peaks_dge$samples <- mutate(reads_in_peaks_dge$samples,'group'=ifelse(rownames(reads_in_peaks_dge$samples) %like% 'C','chimp','human'))

peak_width <- copy(orthologous_peaks)%>%as.data.table()
peak_width <- peak_width[peakID %in% rowData(reads_in_peaks_filtered)$peakID]%>%dplyr::pull('width')
reads_in_peaks_dge$counts = rpkm(reads_in_peaks_dge$counts,peak_width) 

## quantile norm with voom
quant_norm <- function(filtered_counts){
    qn_counts = copy(filtered_counts)
    qn_counts$counts = normalizeQuantiles(qn_counts$counts, ties=TRUE)
    return(qn_counts)
}

quant_norm_counts = quant_norm(reads_in_peaks_dge)

## keep row DA matrix corresponding to DA peaks associated with DE genes
da_peak_de_genes <- copy(one_peak_per_gene)[DA=='da' & DE=='de']

subset_peak_count_matrix <- copy(quant_norm_counts)
subset_peak_count_matrix$counts <- subset(subset_peak_count_matrix$counts,rownames(subset_peak_count_matrix$counts) %in% da_peak_de_genes$peakID)

peak_counts_df <- data.table(subset_peak_count_matrix$counts)[
    ,peakID:=rownames(subset_peak_count_matrix$counts)
]

## add chrom state info to peaks 
chrom_state_dir <- '../data/iPSC_chrom_states_hg38'

ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

da_peak_de_genes_chromstate <- foverlaps(copy(da_peak_de_genes),ipsc_chromstate,type='any')[
  ,c(range_keys[-1]):=NULL
]%>%na.omit()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))

da_peak_de_genes_chromstate <- da_peak_de_genes_chromstate[
    ,pleiotropy:=.N,by=.(peakID,chrom_state)
    ][
      ,.SD[which.max(abs(pleiotropy))], by=.(peakID)
      ][
          ,peak_chromstate:=paste(peakID,chrom_state,sep='.')
          ][
            ,c('chrom_state','cell_type'):=NULL
]%>%unique()

peak_counts_df <- peak_counts_df[da_peak_de_genes_chromstate,on='peakID',nomatch=0]%>%dplyr::select(c(1:12,'peak_chromstate'))

mat_rownames = peak_counts_df$peak_chromstate
peakcounts_mat <- peak_counts_df[,peak_chromstate:=NULL]%>%as.matrix()
rownames(peakcounts_mat)=mat_rownames

## heatmap
pdf(paste(outplot_dir,"heatmap_quantnormd_counts_peaks_w_de_genes.pdf",sep=''),width=10,height=7)
x=round(log2(t(peakcounts_mat)+1),1)
y=table(x)
x[x > 5.4] <- 6
names(chrom_state_colors) = chrom_states

Heatmap(
    x,
    show_column_dend = F,
    show_column_names = F,
    column_order = order(as.factor(readr::parse_number(gsub("^.*\\.", "",colnames(x))))),
    # column_split = as.factor(readr::parse_number(gsub("^.*\\.", "", colnames(x)))),
    top_annotation = HeatmapAnnotation(
    chrom_state = anno_simple(
        gsub("^.*\\.", "", colnames(x)),
        border=T,
        height = unit(1,'cm'),
        col=chrom_state_colors),
        show_annotation_name = F
    ),
    name='log2 rpkm',
    col=viridis::viridis(7)
)
dev.off()


## write file with all info to make visualization script much neater/quicker
write.table(da_peak_de_genes_chromstate,paste(tmp_files_dir,'da_peaks_w_de_genes_gviz.txt',sep=''),sep='\t',col.names=T,quote=F,row.names = F)