## rerun da analysis 
## script to perform DA analysis using csaw
## here i use RPKM 
library(data.table)
library(magrittr)
library(dplyr)
library(edgeR)
library(GenomicRanges)
library(csaw)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(viridis)
library(circlize)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

peak_file <- paste('output/final_peak_set/',genome,'_all_orthologous_peaks.txt',sep='')

target_genes_dir <- './output/files/GO_enrich'
outplot_dir <- create_dir(plot_dir,'rna_seq')

## list bam
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(pe = "both",restrict=standard_chr,minq=20, dedup=TRUE)

get_bams <- function(species){
    bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}

bams <- get_bams(genome)

## read consensus peak
peaks <- fread(peak_file,sep='\t',header=T)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

## count reads in peaks (rpkm)
reads_in_peaks <- regionCounts(bams, peaks,param=param)

## remove peaks with low counts
abundances <- aveLogCPM(asDGEList(reads_in_peaks))
keep <- abundances > 0

reads_in_peaks_filtered <- reads_in_peaks[keep,]

reads_in_peaks_dge <- asDGEList(reads_in_peaks_filtered)
rownames(reads_in_peaks_dge) <- rowData(reads_in_peaks_filtered)$peakID
colnames(reads_in_peaks_dge) <- samples_names
reads_in_peaks_dge$samples <- mutate(reads_in_peaks_dge$samples,'group'=ifelse(rownames(reads_in_peaks_dge$samples) %like% 'C','chimp','human'))

peak_width <- copy(peaks)%>%as.data.table()
peak_width <- peak_width[peakID %in% rowData(reads_in_peaks_filtered)$peakID]%>%dplyr::pull('width')
reads_in_peaks_dge$counts = rpkm(reads_in_peaks_dge$counts,peak_width) 

##------------------
## Normalisation
##------------------
## quantile norm with voom
quant_norm <- function(filtered_counts){
    qn_counts = copy(filtered_counts)
    qn_counts$counts = normalizeQuantiles(qn_counts$counts, ties=TRUE)
    return(qn_counts)
}

quant_norm_counts = quant_norm(reads_in_peaks_dge)

## read target genes and keep only the DE ones 
target_genes <- list.files(target_genes_dir,full.names= T,recursive=F,pattern='target')%>%lapply(
  function(x)fread(x,sep='\t',header=T))%>%rbindlist()

## add chrom state info to peaks 
chrom_state_dir <- '../data/iPSC_chrom_states_hg38'

ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
setkeyv(ipsc_chromstate,range_keys)

target_genes_annotation <- foverlaps(target_genes,ipsc_chromstate,type='any')[
  ,c(range_keys[-1]):=NULL
]%>%na.omit()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  
target_genes_annotation <- target_genes_annotation[
    ,pleiotropy:=.N,by=.(peakID,chrom_state)
    ][
      ,.SD[which.max(abs(pleiotropy))], by=.(peakID)
      ][
          ,peak_chromstate:=paste(peakID,chrom_state,sep='.')
]

## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
de_genes <- fread(
  "../rna_seq/de_output/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,logFC:=-logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','t','B','genes'):=NULL
]

ensembl_86_ids <- fread("../rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

target_genes_ensembl <- copy(target_genes_annotation)[ensembl_86_ids,on='gene',nomatch=0][de_genes,on='EnsemblID',nomatch=0]%>%unique()
target_de_genes <- copy(target_genes_ensembl)[DE=='de']

## keep row DA matrix corresponding to peaks with target genes
test <- copy(quant_norm_counts)
test$counts <- subset(test$counts,rownames(test$counts) %in% target_de_genes$peakID)

test_df <- data.table(test$counts)[
    ,peakID:=rownames(test$counts)
    ][
        target_de_genes,on='peakID',nomatch=0
]%>%dplyr::select(c(1:12,'peak_chromstate'))
mat_rownames = test_df$peak_chromstate
test_mat <- test_df[,peak_chromstate:=NULL]%>%as.matrix()
rownames(test_mat)=mat_rownames

## heatmap
pdf(paste(outplot_dir,"heatmap_quantnormd_counts_peaks_w_de_genes.pdf",sep=''),width=10,height=7)
x=round(log2(t(test_mat)+1),1)
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
    col=viridis(7)
)
dev.off()
