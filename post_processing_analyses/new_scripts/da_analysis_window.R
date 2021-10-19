## script to perform DA analysis using csaw
library(data.table)
library(magrittr)
library(dplyr)
library(edgeR)
library(GenomicRanges)
library(csaw)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(liftOver)
library(rtracklayer)
library(mixOmics)
library(ComplexHeatmap)
library(viridis)

options(width=150)

## directories
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

bamDir  = "/output/Post_alignment/Files"
result_dir = './post_processing_analyses/output'
peakDir = '/output/PeakCalling/Files/'
plot_dir = paste0(result_dir,'/plots/DA',sep='')
peak_outdir = paste0(result_dir,'/DA/peaks/',sep='')

genome = 'hg38'

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

## list bam
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(pe = "both",restrict=standard_chr,minq=20, dedup=TRUE)

get_bams=function(species){
    bams=list.files(paste0(species,bamDir), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}

bams <- get_bams(genome)

## get distribution fragment sizes 
# get_distr_fragments <- function(list_bams){
#     fragment_sizes <- lapply(list_bams,function(x)getPESizes(x,param=param)$sizes)
#     names(fragment_sizes) <- sample_names

#     plots <- purrr::map2(fragment_sizes,names(fragment_sizes),function(f,n){
#         ggplot()+aes(f)+
#         geom_density()+
#         xlab(n)+ylab('frequency')
#         }
#     )
#     figure <- ggarrange(plotlist=plots,ncol=2,nrow=6)

#     return(figure)

# }

# hg38_fragment_sizes <- get_distr_fragments(hg38_bams)
# pantro5_fragment_sizes <- get_distr_fragments(pantro5_bams)

# pdf(paste(plot_dir,'hg38_fragment_size_distribution.pdf',sep=''),width=7,height=7)
# hg38_fragment_sizes
# dev.off()

# pdf(paste(plot_dir,'pantro5_fragment_size_distribution.pdf',sep=''),width=7,height=7)
# pantro5_fragment_sizes
# dev.off()


## count reads in sliding window
window_counts <- windowCounts(bams, width=1000, param=param)
windows <- rowRanges(window_counts)%>%as.data.table()
windows <- windows[,windowID:=paste('window_',1:nrow(windows),sep='')]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

## filter out low counts (i.e. uninteresting windows) by local enrichment
neighbor_regions <- suppressWarnings(resize(windows, width=2000, fix="center")) # change width parameter as desired
wider_window_counts <- regionCounts(bams, regions=neighbor_regions, param=param) # count reads in neighborhoods

filter_stat <- filterWindowsLocal(window_counts, wider_window_counts) 

pdf(export_file(plot_dir,genome,"local_enrichment_counts.pdf"),width=7,height=7)
hist(filter_stat$filter, xlab="Log-fold change from local background", 
    breaks=100, main="", col="grey80", xlim=c(0, 5))
abline(v=log2(4), col="red", lwd=2)
dev.off()

## filter counts
## threshold of 3-fold increase in enrichment over 2kb neighborhood abundance
window_counts_dge <- asDGEList(window_counts)
rownames(window_counts_dge) <- windows@elementMetadata$windowID
colnames(window_counts_dge) <- samples_names
window_counts_dge$counts <- window_counts_dge$counts[filter_stat$filter > log2(4),] ## filtering step

density_plot <- function(counts){
    plotDensities(log(counts$counts+1), group=rownames(counts$samples),col=samples_palette) 
}

pdf(export_file(plot_dir,genome,"logcounts.pdf"),width=7,height=7)
density_plot(window_counts_dge)
dev.off()

## quantile normalization with voom
quant_norm <- function(filtered_counts){
    qn_counts = copy(filtered_counts)
    qn_counts$counts = normalizeQuantiles(qn_counts$counts, ties=TRUE)
    qn_counts$samples <- mutate(qn_counts$samples,'group'=ifelse(rownames(qn_counts$samples) %like% 'C','chimp','human'))
    return(qn_counts)
}

quant_norm_counts = quant_norm(window_counts_dge)

pdf(export_file(plot_dir,genome,"logcounts_quantnormd.pdf"),width=7,height=7)
density_plot(quant_norm_counts)
dev.off()

##-------------------
## Run PCA 
##-------------------
run_pca <-function(quant_norm_counts){
    pca <- pca(t(quant_norm_counts$counts),ncomp=5,scale=F) 
    group = c(rep('chimp',6),rep('human',6))
    return(pca)
}
pca <- run_pca(quant_norm_counts)

pdf(export_file(plot_dir,genome,"pca.pdf"),width = 7, height = 7)
plotIndiv(pca,comp = c(1,2),legend=T,ind.names = T)
dev.off()

## heatmap
plot_heatmap <- function(quant_norm_counts){
    Heatmap(
        cor(quant_norm_counts$counts,method='pearson'), 
        name = "Pearson corr", col=viridis(100)
        ) 
}

pdf(export_file(plot_dir,genome,"heatmap.pdf"),width=7,height=7)
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
windows_df <- copy(windows)%>%as.data.table()

final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'windowID'= rownames(final_results))%>%as.data.table()

final_results <-  final_results[
        ,FDR:=p.adjust(final_results$PValue,method = 'fdr')
        ][
            ,DA:=as.factor(ifelse(FDR<=0.05,'significant','non_significant'))
            ][
                windows_df,on='windowID',nomatch=0
            ][
        ,c('width','strand'):=NULL
] 


outfile = export_file(paste(result_dir,'/DA',sep=''),genome,'da_results.txt')
write.table(final_results,outfile,sep='\t',col.names=T,row.names = F,quote=F)

##--------------------
## Generate MA plot
##--------------------
ma_plot = function(df){
    plot <- ggplot(df,aes(x = logCPM, y = logFC, col = DA)) + 
        geom_point() + scale_color_manual(values = c("black", "red")) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess",se = FALSE) + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)
        return(plot)
}

pdf(export_file(plot_dir,genome, "ma_plot.pdf"),width = 7, height = 7)
ma_plot(final_results)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x = logFC,y =-log10(FDR), col = DA))+
    scale_color_manual(values = c("black", "red")) + 
    geom_hline(yintercept=1.3, linetype='dashed', color='black', size=1)
    # +facet_wrap(species~.,ncol=3)
    
    return(plot)
}

pdf(export_file(plot_dir,genome,"volcano_plot.pdf"),width = 10, height = 7)
volcano_plot(final_results)
dev.off()
