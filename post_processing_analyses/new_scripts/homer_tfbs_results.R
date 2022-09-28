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
library(GenomicRanges)
library(openxlsx)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

da_peaks_dir <- './output/DA/peaks/'
homer_dir <- './output/homer/homer_output/'
outplot_dir <- create_dir(plot_dir,'tfbs')
outfile_dir <- create_dir(outdir,'tables/tfbs')

##-------------
## read files
##-------------
## homer results
homer_files <- list.files(homer_dir,recursive=T,full.names=T,pattern='knownResults.txt')

read_homer_results <- function(file){
  results <- fread(file,sep='\t',header=T)[,gene_symbol:=gsub("\\/.*","",`Motif Name`)]
  motif_family <- stringr::str_extract(string = results$gene_symbol,pattern = "(?<=\\().*(?=\\))")
  results <- results[
    ,motif_family:=motif_family
    ][
      ,gene_symbol:=gsub("\\(.*","",gene_symbol)
      ][
        ,prop_motif_in_target:=as.numeric(gsub("\\%.*","",`% of Target Sequences with Motif`))
        ][
          ,c('Motif Name','Consensus','% of Target Sequences with Motif'):=NULL
          ]%>%setnames(old=c(1:3),new=c('pval','logpval','qval'))
  return(results)
}

homer_results <- read_homer_results(homer_files)

# ## all homer motifs
# all_homer_motifs<- universalmotif::read_homer('/data/projects/punim0595/dvespasiani/homer/downloaded_homer.motifs')
# all_homer_motifs_dt <- copy(all_homer_motifs)%>%lapply(function(y)data.table(motif_name=y@name))%>%rbindlist()%>%unique()

##----------------
## filter motifs
##----------------
## remove motifs with fdr adj.pval > 0.01 
## remove motifs present in less than 5% of target sequences
## then modify motif names to retrieve their ensembl IDs
filtered_motifs <- copy(homer_results)[
        ,padj:=p.adjust(pval,'fdr')
        ][
          ,log10pval:=-log10(ifelse(padj==0, 1e-320,padj)) ## done this to control R from converting <<< 0 values to 0
          ][
            log10pval > -log10(1e-15) #homer suggested threshold
              ][
                prop_motif_in_target>=10 & motif_family !='?'
                ][
                  !gene_symbol %like% 'Unknown|OCT4-SOX2-TCF-NANOG' ## *** ##
                  ][
                    ,gene_symbol:=ifelse(gene_symbol %like% 'BORIS','CTCFL',gene_symbol)
                    ][
                      ,gene_symbol:=toupper(gene_symbol)
                      ][
                        ,gene_symbol:= gsub("\\-.*","",gene_symbol)
                        ][
                          ,gene_symbol:= gsub("\\:.*","",gene_symbol)
                          ][
                            ,gene_symbol:= gsub("\\+.*","",gene_symbol)
                            ][
                              ,gene_symbol:= gsub("\\..*","",gene_symbol)
]

## export filtered motifs in excel file for supplementary
write.xlsx(filtered_motifs[,c('gene_symbol','motif_family','pval','padj')],paste(outfile_dir,'homer_filtered_motifs.xlsx',sep=''),append=F,overwrite=T)

###  ** those motifs are already invidually enriched plus I individually analyse them (see below) ** ###

## numb and prop
numb_initial_motifs <- length(unique(homer_results$gene_symbol))
numb_filtered_motifs <-length(unique(filtered_motifs$gene_symbol))
round(numb_filtered_motifs/numb_initial_motifs *100,2)
#  35.18
length(unique(filtered_motifs$motif_family))
# [1] 53

##-----------------------
## read expression data 
##-----------------------
library(biomaRt)
rnaseq_dir = '../rna_seq/'
counts_dir = 'output/PostAlignment/'

exon_lengths <- fread(paste(rnaseq_dir,'orthoexon_files/metaOrthoExon_hg38_panTro5_ensembl86_0.97_ortho_FALSE_pc.txt',sep=''),header=T)

read_expression <- function(genome,file,colname){
  featurecounts_results <- fread(paste(rnaseq_dir,genome,counts_dir,file,sep=''), header=T)
  gene_names <- featurecounts_results$GeneID
  cpm_counts <- featurecounts_results[,GeneID:=NULL]%>%as.matrix()%>%edgeR::cpm()%>%as.data.table()
  cpm_counts <- cpm_counts[,mean_expr:=round(rowMeans(.SD),2)][,EnsemblID:=gene_names]%>%setnames(old='mean_expr',new=colname)
  return(cpm_counts)
}

chimp_cpm_counts <- read_expression('./pantro5/',"ipsc_chimp_counts.txt",'mean_expr_chimp')
human_cpm_counts <- read_expression('./hg38/',"ipsc_human_counts.txt",'mean_expr_human')

ipsc_counts <- merge(chimp_cpm_counts,human_cpm_counts,by="EnsemblID", all=T)

## filter ipsc_counts by the number of genes retained at the end of DE testing
de_info <- fread(
  "../rna_seq/de_output/irenes_files/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,de_logFC:=logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','AveExpr','t','B','genes','logFC'):=NULL
]

ipsc_counts<- merge(ipsc_counts,de_info,on='EnsemblID')

## get gene symbols
# ensembl.mart <- useMart(
#     # host='https://dec2016.archive.ensembl.org',
#     biomart='ENSEMBL_MART_ENSEMBL', 
#     dataset='hsapiens_gene_ensembl'
# )

# ensembl_id_gene_symbol <- getBM(
#     attributes = c('ensembl_gene_id','hgnc_symbol'), 
#     filters = 'ensembl_gene_id', 
#     values = ipsc_counts$EnsemblID,
#     mart = ensembl.mart
# )%>%as.data.table()
# names(ensembl_id_gene_symbol) = c('EnsemblID','gene_symbol')

library(org.Hs.eg.db) # remember to install it if you don't have it already
symbols <- mapIds(org.Hs.eg.db, keys = ipsc_counts$EnsemblID, keytype = "ENSEMBL", column="SYMBOL")
symbols <- data.table(EnsemblID=names(symbols),gene_symbol=symbols)

## add expression info
filtered_motifs_expr <- copy(filtered_motifs)[,c('log10pval','gene_symbol','motif_family','prop_motif_in_target')]%>%full_join(symbols,by='gene_symbol')%>%as.data.table()  
filtered_motifs_expr <- filtered_motifs_expr[,keep:=ifelse(is.na(log10pval),'no','yes')][keep=='yes'][,keep:=NULL]

mean_expr <- copy(ipsc_counts)%>%dplyr::select(c('EnsemblID','adj.P.Val',contains('mean')))

filtered_motifs_expr <- merge(filtered_motifs_expr,mean_expr,by='EnsemblID',all.x=T)
filtered_motifs_expr <- filtered_motifs_expr[,tfsymbol_fam:=paste(gene_symbol,motif_family,sep=':')]

## correlation TF enrichment TF expression
tf_w_expr_only <- copy(filtered_motifs_expr)%>%na.omit()

cor.test(tf_w_expr_only$log10pval,tf_w_expr_only$mean_expr_human)
# 	Pearson's product-moment correlation

# data:  tf_w_expr_only$log10pval and tf_w_expr_only$mean_expr_human
# t = 1.3528, df = 75, p-value = 0.1802
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.07213612  0.36567449
# sample estimates:
#       cor 
# 0.1543363 

cor.test(tf_w_expr_only$log10pval,tf_w_expr_only$mean_expr_chimp)
# 	Pearson's product-moment correlation

# data:  tf_w_expr_only$log10pval and tf_w_expr_only$mean_expr_chimp
# t = 1.1593, df = 75, p-value = 0.25
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.09409951  0.34636179
# sample estimates:
#       cor 
# 0.1326756 

##--------------------
## plot results 
##--------------------
# pvals <- copy(genwide_motif_enrichment)[,c('DE','tf_symbol')][,p_sign:=ifelse(DE=='de','*',NA)]%>%dplyr::pull('p_sign')
# names(pvals) = genwide_motif_enrichment$tf_symbol

motifs_expr <- copy(filtered_motifs_expr)[log10pval>-log10(1e-200)]

matzscore_rownames <- copy(motifs_expr$tfsymbol_fam)
mat_zscore <- copy(motifs_expr)[,log10pval]%>%as.matrix()
rownames(mat_zscore) = matzscore_rownames

# colfun_zscore = circlize::colorRamp2(c(min(mat_zscore), 0, max(mat_zscore)), c("#0582CA", "white", "#d62828"))

## cluster colors
qual_col_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_palette$maxcolors, rownames(qual_col_palette)))

cluster_colors = sample(col_vector,length(unique(gsub(".*:","",rownames(mat_zscore)))))
names(cluster_colors) = unique(gsub(".*:","",rownames(mat_zscore)))

## de pval annotation
de_pval <- copy(motifs_expr$adj.P.Val)
de_pval[is.na(de_pval)] <- 1
de_pval <- -log10(de_pval)
names(de_pval) <- gsub(".*:","",rownames(mat_zscore))

heatmap_zscore <- Heatmap(
  mat_zscore,
  name='log10 FDR-adjusted pvalue',
  col=magma(10),
  heatmap_legend_param = list(direction = "horizontal")
)

mat_expr <- copy(motifs_expr)[,c('mean_expr_chimp','mean_expr_human')]%>%as.matrix()
mat_expr <- log2(mat_expr+2)
mat_expr[is.na(mat_expr)]<-0
rownames(mat_expr) = matzscore_rownames

colfun_expr = circlize::colorRamp2(c(seq(min(mat_expr),max(mat_expr),2.5)), c('gray82',"#FFFF3F",'#BFD200',"#2B9348"))

col_depval <- circlize::colorRamp2(c(seq(min(de_pval),max(de_pval),2.5)), c('gray82',c("#00b4d8","#0077b6","#023e8a","#03045e")))

# c("#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4","#90e0ef","#ade8f4","#caf0f8")


# depval_colors <- viridis(length(de_pval))
# names(depval_colors) = names(de_pval)

heatmap_expr <- Heatmap(
  mat_expr,
  name='avg gene expression \n (log2 CPM)',
  col=colfun_expr,
  heatmap_legend_param = list(direction = "horizontal"),
  right_annotation = HeatmapAnnotation(
    which='row',
    de_pvalue = anno_simple(de_pval,width = unit(1,'cm'),col=col_depval),
    cluster = anno_simple(gsub(".*:","",rownames(mat_zscore)),width = unit(1,'cm'),col=cluster_colors),
    show_legend=T 
    )
)

lgd_pvalue = Legend(
  title = "log10 DE adj.p", 
  col_fun = col_depval, at = c(seq(min(de_pval),max(de_pval),2.5)), 
  direction = "horizontal"
)

pdf(paste(outplot_dir,'heatmap_expressed_enriched_motifs.pdf',sep=''),width=7,height = 10)
draw(
  heatmap_zscore+heatmap_expr,
  annotation_legend_list = list(lgd_pvalue),
  heatmap_legend_side = "bottom", 
  annotation_legend_side = "bottom"
)
dev.off()

##------------------------
## pluripotency motifs
##------------------------
pluripotency_motif_files <- list.files('output/homer/homer_output',full.names=T,recursive=T,pattern=c('discovered.txt'))
plurip_filnames = gsub("\\/.*","",stringr::str_remove(pluripotency_motif_files,'output/homer/homer_output/'))

pluripotency_motifs<- lapply(pluripotency_motif_files,function(x){
  results <- fread(x,sep='\t',header=T,drop=c("Motif Name",'Strand'))%>%setnames(old=1,new='peakID')
  results <- results[,motif_length:=nchar(Sequence)][,rank_score:=round(MotifScore/max(MotifScore),2)]
  return(results)
  }
)
names(pluripotency_motifs) = plurip_filnames
pluripotency_motifs <- Map(mutate,pluripotency_motifs,tf=names(pluripotency_motifs))%>%rbindlist()

## add genomic region
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','da_species','logFC')]

pluripotency_motifs <- pluripotency_motifs[
  da_results,on='peakID',nomatch=0
  ][
    ,peak_center:=start+(end-start)
    ][
      ,motif_start:=peak_center+Offset
      ][
        ,motif_end:=motif_start+motif_length
]

## plot distribution number of motifs by rank scores
pdf(paste(outplot_dir,'distribution_ranked_motif_scores.pdf',sep=''),width=10,height=7)
df <- copy(pluripotency_motifs)[,c(..range_keys,'rank_score')]%>%unique()
ggplot(df,aes(x=rank_score))+geom_bar()+
ylab('number of motifs')+
scale_x_continuous(
  name ='Ranked motif score',
  limit =  range(df$rank_score),
  breaks =  seq(0,1,0.1))+
theme_classic()+
theme()
dev.off()

## proportion of original peaks with 
numb_peaks_w_plurip_tf <- copy(top_pluripotency_motifs$peakID)%>%unique()%>%length()
numb_peaks <- copy(da_results$peakID)%>%length()
numb_peaks_w_plurip_tf/numb_peaks*100
#[1] 79.69008

tfs <- copy(pluripotency_motifs)%>%split(by='tf')%>%lapply(function(x)x=x[,peakID]%>%unique())
plurip_tf_palette <- c('#c1121f','#723d46','#e29578','#606c38','lightgray')
names(plurip_tf_palette) = c('nanog','oct4','sox2','ctcf','other')

pdf(paste(outplot_dir,'pluripotency_tf_upset.pdf',sep=''),width=10,height=7)
UpSetR::upset(UpSetR::fromList(tfs),nsets = 6,order.by = "freq")
dev.off()

## look whether they are enriched within bivalent states
chrom_state_dir <- '../data/iPSC_chrom_states_hg38'
ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

## overlap ranges
annotate_peaks <- function(peaks){
  annotated <- foverlaps(copy(peaks),ipsc_chromstate,type='any')%>%na.omit()
  annotated <- annotated[
    ,overlap:=ifelse(i.start<start,i.end-start,end-i.start),by=.(cell_type)
    ][
        ,.SD[which.max(overlap)], by=.(peakID,cell_type)
        ][
            ,c(range_keys[-1],'overlap'):=NULL
            ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  return(annotated)
}

mypeaks_annotation <- annotate_peaks(da_results)
mypeaks_annotation <- mypeaks_annotation[,pluripotency:=ifelse(peakID %in% pluripotency_motifs$peakID,'yes','no')]

count_peaks_chromstate <- function(peaks){
  counts <- copy(peaks)
  counts<-counts[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state,cell_type)
      ][
        ,numb_peaks:=.N,by=.(cell_type)
        ][
            ,c('cell_type','chrom_state','numb_peaks_chromstate','numb_peaks')
            ]%>%unique()
  return(counts)
}

mypeaks_plurip_tf_counts <- count_peaks_chromstate(mypeaks_annotation[pluripotency=='yes'])%>%split(by='cell_type')
mypeaks_other_tf_counts <- count_peaks_chromstate(mypeaks_annotation[pluripotency!='yes'])%>%split(by='cell_type')

plurip_tf_chromstate_or <- purrr::map2(mypeaks_plurip_tf_counts,mypeaks_other_tf_counts,function(x,y){
  or<-calculate_or(x,y,'chrom_state')%>%adjust_pvalues()
  or <- or[,6:=NULL]
  return(or)
})
plurip_tf_chromstate_or <- Map(mutate,plurip_tf_chromstate_or,cell_type=names(plurip_tf_chromstate_or))%>%rbindlist()

## plot the results
plot_or <- function(peaks){
  p <- ggplot(peaks,aes(x=factor(elements,levels=chrom_states),y=log(odds_ratio),fill=elements))+
    geom_violin(trim=T,scale = "width")+
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=0.2)+
    scale_fill_manual(values = chrom_state_colors)+
    geom_hline(yintercept=0,linetype='dashed')+
    xlab('chromatin state') + ylab('log OR')+
    theme_classic()+
    theme(
      legend.position = "bottom",
      axis.text.x =element_blank(),
      axis.ticks.x =element_blank()
    )
  return(p)
}

pdf(paste(outplot_dir,'or_peaks_w_plurip_tf_vs_peaks_without_chromstate.pdf',sep=''),width=7,height = 7)
plot_or(plurip_tf_chromstate_or)
dev.off()

## check whether these peaks are more conserved compared to all the others
library(phastCons7way.UCSC.hg38)

da_results <- da_results[,pluripotency:=ifelse(peakID %in% pluripotency_motifs$peakID,'yes','no')]

## expand peaks (i.e. get all bp within peaks from start to end) 
get_phastCons_score <-function(df){
    expanded_peaks <- copy(df)[,c(..range_keys,'peakID')]%>%unique()
    expanded_peaks <- expanded_peaks[,list(start = seq(start, end)), by = peakID] 
    expanded_peaks <- expanded_peaks[
        df[,c('seqnames','peakID')],on='peakID',allow.cartesian=T
        ][
            ,end:=start
    ]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

    phastCons <- phastCons7way.UCSC.hg38

    scores <- gscores(phastCons, expanded_peaks)%>%as.data.table()
    rm(expanded_peaks)
    scores <- scores[
            ,avg_phastcons:=mean(default),by=.(peakID)
            ][
                ,c('peakID','avg_phastcons')
    ]%>%unique()%>%na.omit()
    final_df <- copy(df)[scores,on='peakID',nomatch=0]
    return(final_df)
}

conservation_scores = get_phastCons_score(da_results)
conservation_peaks_w_tf <- copy(conservation_scores)[pluripotency=='yes'][pluripotency_motifs[,c('tf','peakID')],on='peakID',nomatch=0]
conservation_peaks_wout_tf <- copy(conservation_scores)[pluripotency=='no'][,tf:='other']

## plot conservation score 
comparisons = list(
  c("nanog",'other'),
  c('oct4','other'),
  c('ctcf','other'),
  c('sox2','other')
)
pdf(paste0(outplot_dir,'avg_phastcons_peaks_w_pluriptf.pdf',sep=''),width = 7, height = 7)
df <- rbind(conservation_peaks_w_tf,conservation_peaks_wout_tf)
tf_order <- c('nanog', 'oct4','ctcf','sox2','other')
ggplot(df,aes(x=factor(tf,level=tf_order),y=avg_phastcons,fill=tf))+
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
  scale_fill_manual(values=plurip_tf_palette)+
  xlab(' ')+ylab('avg phastCons')+
  stat_compare_means(
  method = "wilcox.test",
  comparisons = comparisons,
  size=5
  )+
  theme_classic()+
  theme(
    legend.position = "bottom",
    axis.text.x =element_blank(),
    axis.ticks.x =element_blank()
  )
dev.off()

## look if DA peaks are enriched for these TFs (by permutation)
permute_devtfs <- function(x){
  numb_permutations = 1000
  tf_da_observed_value = as.data.table(table(copy(unique(x[,c('peakID','tf','DA')][DA=='da']))$tf))%>%setnames(old=c(1:2),new=c('tf','score'))%>%split(by='tf')

  tf_da_list_permutations = list()
  for(i in 1:numb_permutations){
        new_da_col <- sample(x$DA,nrow(x),replace=F)
        new_df <- copy(x)[,DA:=new_da_col]
        permuted_scores=table(unique(new_df[,c('peakID','tf','DA')][DA=='da'])$tf)
        tf_da_list_permutations[[i]] = as.data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('tf','permuted_score'))
  }
  tf_da_permutation <- rbindlist(tf_da_list_permutations)%>%split(by='tf')

  enrichment = purrr::map2(tf_da_permutation,tf_da_observed_value,function(p,o){
    zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
    pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
    stat =  data.table(zscore = zscore,pval=pvalue)
    return(stat)
    }
  )
  enrichment <- Map(mutate,enrichment,tf=names(enrichment))%>%rbindlist()%>%setorderv('zscore',1)
  enrichment <- enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]
  
  return(enrichment)
}
# numb_permutations = 10000
# tf_da_observed_value = as.data.table(table(copy(unique(top_pluripotency_motifs[,c('peakID','tf','DA')][DA=='da']))$tf))%>%setnames(old=c(1:2),new=c('tf','score'))%>%split(by='tf')

# tf_da_list_permutations = list()
# for(i in 1:numb_permutations){
#       new_da_col <- sample(pluripotency_motifs$DA,nrow(pluripotency_motifs),replace=F)
#       new_df <- copy(pluripotency_motifs)[,DA:=new_da_col]
#       permuted_scores=table(unique(new_df[,c('peakID','tf','DA')][DA=='da'])$tf)
#       tf_da_list_permutations[[i]] = as.data.table(permuted_scores)[,permutation:=i]%>%setnames(old=c(1:2),new=c('tf','permuted_score'))
# }
# tf_da_permutation <- rbindlist(tf_da_list_permutations)%>%split(by='tf')

# tf_da_enrichment = purrr::map2(tf_da_permutation,tf_da_observed_value,function(p,o){
#   zscore = (o$score-mean(p$permuted_score))/sd(p$permuted_score)
#   pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
#   stat =  data.table(zscore = zscore,pval=pvalue)
#   return(stat)
#   }
# )
# tf_da_enrichment <- Map(mutate,tf_da_enrichment,tf=names(tf_da_enrichment))%>%rbindlist()%>%setorderv('zscore',1)
# tf_da_enrichment <- tf_da_enrichment[,abslog_zscore:=log(abs(zscore))][,log_zscore:=ifelse(zscore<0,-abslog_zscore,abslog_zscore)]

alltfs_daenrichment <- permute_devtfs(pluripotency_motifs) 
toptfs_daenrichment <- permute_devtfs(top_pluripotency_motifs) 

plot_enrichments <- function(enrich_res,column,palette,xlab){
  df <- copy(enrich_res)[,column_to_reorder:=column]
  p <- ggplot(df,aes(x=reorder(column_to_reorder,-log_zscore),y=log_zscore,fill=column_to_reorder))+
    geom_bar(stat='identity',position =  position_dodge())+
    scale_fill_manual(values = palette)+
    geom_hline(yintercept=0,linetype='dashed')+
    xlab(xlab) + ylab('log zscore')+
    theme_classic()+
    theme(
      legend.position = "bottom",
      axis.text.x =element_blank(),
      axis.ticks.x =element_blank()
    )
  return(p)
}


pdf(paste(outplot_dir,'enrichm_pluripotency_tf_davsnonda.pdf',sep=''),width=8,height = 5)
plot_enrichments(enrich_res=alltfs_daenrichment,column=alltfs_daenrichment$tf,palette=plurip_tf_palette[-5],xlab='TF')
dev.off()


pdf(paste(outplot_dir,'enrichm_top_pluripotency_tf_davsnonda.pdf',sep=''),width=8,height = 5)
plot_enrichments(enrich_res=toptfs_daenrichment,column=toptfs_daenrichment$tf,palette=plurip_tf_palette[-5],xlab='TF')
dev.off()

# ## for each of these TFs make volcano plot of logFC
# volcano_plot <-function(df){
#     plot <- ggplot(df) + 
#     geom_point(aes(x = logFC,y =-log10(FDR), col = tf),alpha=0.2)+
#     # scale_color_manual(values = da_palette) + 
#     geom_hline(yintercept=-log10(0.01), linetype='dashed', color='black', size=0.5)+
#     geom_vline(xintercept=0, linetype='dashed', color='black', size=0.5)+
#     scale_colour_manual(values=plurip_tf_palette)+
#     theme_classic()+
#     xlim(-8,+8)+
#     theme(
#         legend.position = "bottom",
#         axis.ticks.x =element_blank()
#         )
#     return(plot)
# }

# ctcf_motifs <- copy(pluripotency_motifs)[tf=='ctcf']
# nanog_motifs <- copy(pluripotency_motifs)[tf=='nanog']
# oct4_motifs <- copy(pluripotency_motifs)[tf=='oct4']
# sox2_motifs <- copy(pluripotency_motifs)[tf=='sox2']

# pdf(paste(outplot_dir,'ctcf_volcano.pdf',sep=''),width=7,height = 7)
# volcano_plot(ctcf_motifs)
# dev.off()

# pdf(paste(outplot_dir,'nanog_volcano.pdf',sep=''),width=7,height = 7)
# volcano_plot(nanog_motifs)
# dev.off()

# pdf(paste(outplot_dir,'oct4_volcano.pdf',sep=''),width=7,height = 7)
# volcano_plot(oct4_motifs)
# dev.off()

# pdf(paste(outplot_dir,'sox2_volcano.pdf',sep=''),width=7,height = 7)
# volcano_plot(sox2_motifs)
# dev.off()
# # ##----------------
# ## filter motifs
# ##----------------
# ## FDR adjust pvalues
# ## remove motifs with fdr adj.pval > 0.01 
# ## remove motifs present in less than 5% of target sequences

# filtered_motifs <- copy(homer_known_results)%>%lapply(
#   function(x){
#     x<-x[
#       ,log_pval:=ifelse(`Log P-value`< -729,-745,`Log P-value`) ## done this to control R from converting <<< 0 values to 0
#       ][
#         ,log10_fdr_adjP:=-log10(p.adjust(exp(log_pval),'fdr'))
#         ][
#           log10_fdr_adjP > 2
#           ][
#             prop_motif_in_target>5
#             ][
#               !motif_name %like% 'Unknown'
#               ][
#                 ,motif_name:=ifelse(motif_name %like% 'BORIS','CTCFL',motif_name)
#                 ]
#   }
# )

# number_filtered_motifs = copy(filtered_motifs)%>%lapply(function(x)as.data.table(x)[,'motif_name']%>%unique()%>%nrow())
# ##--------------
# ## plot results
# ##--------------
# ## important: GFY is a general factor Y, which refers to a motif for which HOMER doesnt know which if the factor that binds it (see http://homer.ucsd.edu/homer/motif/motifDatabase.html)
# ## however, these motifs are real and some TFs have been described binding them (see http://genesdev.cshlp.org/content/24/14/1479.long)
# ## what I'll do is to report the family for this as 'unknown' and in the text I will mention this and cite the above paper.
# motifs_table = copy(filtered_known_motif)%>%lapply(function(x)x=x[,c('motif_name','motif_family','species','log10_fdr_adjP')])%>%rbindlist()
# motifs_table = motifs_table[,motif_family:=ifelse(motif_family=='?','Unknown',motif_family)][,motif_name_family:=paste(motif_name,motif_family,sep=':')][,c('motif_name','motif_family'):=NULL]

# ## create matrix for heatmap
# matrix = copy(motifs_table)[
#   ,.SD[which.max(log10_fdr_adjP)], by=.(species,motif_name_family)
# ]
# matrix = dcast(matrix, motif_name_family~species,value.var='log10_fdr_adjP')
# matrix[is.na(matrix)]=0
# matrix_rownames = matrix$motif_name_family
# matrix = matrix[,motif_name_family:=NULL]%>%as.matrix()
# rownames(matrix) = matrix_rownames

# filtered_matrix = copy(matrix)
# filtered_matrix=filtered_matrix[rowSums(filtered_matrix)>100, ]
# filtered_matrix[!rownames(filtered_matrix) %like% 'Zic2', ]

# ## get random colors for motif families 
# n_tf_families = length(unique(gsub(".*:","",rownames(filtered_matrix))))

# qual_col_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_palette$maxcolors, rownames(qual_col_palette)))

# tf_family_colors = sample(col_vector,n_tf_families)
# names(tf_family_colors) = unique(gsub(".*:","",rownames(filtered_matrix)))


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
#                       col= tf_family_colors
#                       ),
#                       show_annotation_name = F)
#                       )
# }

# pdf(paste(outplot_dir,'homer_known_heatmap.pdf',sep=''),width= 10,height = 20)
# homer_heatmap(filtered_matrix)
# dev.off()


# ##----------------------------
# ## create motif logo object
# ##----------------------------
# ## for this read again HOMER results but using marge 
# ## this allows to retrieve the PWMs for plotting
# library(marge)
# library(motifStack)

# chimp_known_pwms = read_known_results(paste(homer_dir,'chimp/',sep=''))
# common_known_pwms = read_known_results(paste(homer_dir,'common/',sep=''))
# human_known_pwms = read_known_results(paste(homer_dir,'human/',sep=''))

# my_motifs = list(chimp_known_pwms,common_known_pwms,human_known_pwms)%>%lapply(function(y)y=dplyr::select(y,c('motif_name','motif_pwm'))%>%as.data.table())
# my_motifs = lapply(my_motifs,function(x)x=x[,motif_name:=ifelse(motif_name %like% 'BORIS','CTCFL',motif_name)][motif_name %in% gsub('\\:.*','',rownames(filtered_matrix))])%>%rbindlist()%>%tibble()%>%unique()
# names(my_motifs$motif_pwm)=my_motifs$motif_name

# motif_names=as.list(my_motifs$motif_name)
# motifs_pwm=copy(my_motifs$motif_pwm)

# get_motif_logo=function(pwm,names){
#   motifs=copy(pwm)%>%as.matrix()%>%t()
#   colnames(motifs)=paste('V',1:dim(motifs)[2],sep='')

#   motif_object= new("pcm", mat=motifs, name=names)
#   return(motif_object)
# }

# my_motif_object= purrr::map2(motifs_pwm,motif_names,function(x,y)get_motif_logo(x,y))

# pdf(paste(outplot_dir,'logo_homer_top_known_motifs.pdf',sep=''),width= 15,height = 5)
# lapply(my_motif_object,function(x)plot(x))
# dev.off()

# ## now make a barplot with proportion motif families
# motif_families=data.table(motif_name_family=rownames(filtered_matrix))
# motif_families = motif_families[
#   ,c("motif_name", "motif_family") := tstrsplit(motif_name_family, ":", fixed=TRUE)
#   ][
#     ,'motif_name_family':=NULL
#     ][
#       ,numb_motifs_by_family:=.N,by=.(motif_family)
#       ][
#         ,numb_motifs:=.N
#         ][
#           ,prop_family:=numb_motifs_by_family/numb_motifs
#           ][
#             ,mock_col:='tf_family'
#           ]
# prop_motif_families = copy(motif_families)[,c('prop_family','motif_family','mock_col')]%>%unique()

# pdf(paste(outplot_dir,'prop_families_homer_top_known_motifs.pdf',sep=''),width= 5,height = 10)
# ggplot(prop_motif_families,aes(x=mock_col,y=prop_family,fill=motif_family))+
# geom_bar(position="fill", stat="identity")+
# scale_fill_manual(values = tf_family_colors)+ theme_void()
# dev.off()

