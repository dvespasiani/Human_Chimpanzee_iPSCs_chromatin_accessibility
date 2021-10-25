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
library(mixOmics)
library(ComplexHeatmap)
library(viridis)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

peak_file <- paste('output/final_peak_set/',genome,'_all_orthologous_peaks.txt',sep='')

outfile_dir <- create_dir(da_dir,genome)
outplot_dir <- create_dir(plot_dir,paste('DA/',genome,sep=''))

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
# %>%dplyr::pull('width')
peak_width <- peak_width[peakID %in% rowData(reads_in_peaks_filtered)$peakID]%>%dplyr::pull('width')
reads_in_peaks_dge$counts = rpkm(reads_in_peaks_dge$counts,peak_width) 

##----------
## QCs
##----------
pdf(paste(outplot_dir,genome,'_',"peak_logrpkm_counts.pdf",sep=''),width=7,height=7)
plotDensities(log(reads_in_peaks_dge$counts+1), group=rownames(reads_in_peaks_dge$samples),col=samples_palette) 
dev.off()

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

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(outplot_dir,genome,'_',"peak_logrpkm_counts_quantnormd.pdf",sep=''),width=7,height=7)
plotDensities(log(quant_norm_counts$counts+1), group=rownames(quant_norm_counts$samples),col=samples_palette) 
dev.off()

##--------------
## Run PCA 
##--------------
run_pca <-function(counts){
    pca <- pca(t(counts$counts),ncomp=5,scale=F) 
    group = c(rep('chimp',6),rep('human',6))
    return(pca)
}
pca <- run_pca(quant_norm_counts)

pdf(paste(outplot_dir,genome,'_',"peak_pca.pdf",sep=''),width = 7, height = 7)
plotIndiv(pca,comp = c(1,2),legend=T,ind.names = T)
dev.off()

## heatmap
plot_heatmap <- function(counts){
    Heatmap(
        cor(counts$counts,method='pearson'), 
        name = "Pearson corr", col=viridis(100)
        ) 
}

pdf(paste(outplot_dir,genome,'_',"peak_heatmap.pdf",sep=''),width=7,height=7)
plot_heatmap(quant_norm_counts)
dev.off()

##------------
## DA testing
##------------
design_matrix <- model.matrix(~0+group, data=quant_norm_counts$samples)  
colnames(design_matrix) <- c("chimp", "human")

y <- estimateDisp(quant_norm_counts, design_matrix) 
fit <- glmQLFit(y, design_matrix, robust=TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp, levels=design_matrix))

## get df with final results_df
peaks_df <- copy(peaks)%>%as.data.table()

final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'peakID'= rownames(final_results))%>%as.data.table()

final_results <-  final_results[
        ,FDR:=p.adjust(final_results$PValue,method = 'fdr')
        ][
            ,DA:=as.factor(ifelse(FDR<=0.05,'da','non_da'))
            ][
                peaks_df,on='peakID',nomatch=0
            ][
                ,c('width','strand'):=NULL
                ][
                    ,peak_species := ifelse(DA == 'non_da' , 'common',ifelse(DA=='da'& logFC <0,'chimp','human'))
] 

write.table(final_results,paste(outfile_dir,'da_results.txt',sep=''),sep='\t',col.names=T,row.names = F,quote=F)
##--------------------
## Generate MA plot
##--------------------
ma_plot = function(df){
    plot <- ggplot(df,aes(x = logCPM, y = logFC, col = DA)) + 
        geom_point() + scale_color_manual(values = da_palette) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess",se = F) + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)
        return(plot)
}

pdf(paste(outplot_dir,genome,'_',"peak_ma_plot.pdf",sep=''),width = 7, height = 7)
ma_plot(final_results)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x = logFC,y =-log10(FDR), col = peak_species))+
    scale_color_manual(values = species_palette) + 
    geom_hline(yintercept=1.3, linetype='dashed', color='black', size=1)
    
    return(plot)
}

pdf(paste(outplot_dir,genome,'_',"peak_volcano_plot.pdf",sep=''),width = 10, height = 7)
volcano_plot(final_results)
dev.off()
