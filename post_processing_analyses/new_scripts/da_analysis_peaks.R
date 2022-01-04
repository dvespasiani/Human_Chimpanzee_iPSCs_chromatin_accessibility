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
# peak_file <- paste('../',genome,'/output/PeakCalling/Files/merged_sample_macs2_default_peaks.narrowPeak',sep='')

outfile_dir <- create_dir(da_dir,genome)
outplot_dir <- create_dir(plot_dir,paste('DA/',genome,sep=''))

## list bam
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(pe = "both",restrict=standard_chr,max.frag=1000)

get_bams <- function(species){
    bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}

bams <- get_bams(genome)

## read consensus peak
peaks <- fread(peak_file,sep='\t',header=T)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

# peaks <- fread(peak_file,sep='\t',select=c(1:3),header=F,col.names=c(range_keys))%>%makeGRangesFromDataFrame()%>%reduce()%>%as.data.table()
# peaks <- peaks[,peakID:=paste('peak_',1:nrow(peaks),sep='')][,c('width','strand'):=NULL]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

numb_all_peaks <- length(peaks)

## count reads in peaks (rpkm)
reads_in_peaks <- regionCounts(bams, peaks,param=param)

## remove peaks with low counts
abundances <- aveLogCPM(asDGEList(reads_in_peaks))
keep <- abundances > 0

numb_keep_peaks <- length(keep[keep==TRUE])
prop_keep_peaks <- round(numb_keep_peaks/numb_all_peaks *100,2)
numb_keep_peaks
prop_keep_peaks

reads_in_peaks_filtered <- normFactors(reads_in_peaks[keep,])

reads_in_peaks_dge <- asDGEList(reads_in_peaks_filtered)
rownames(reads_in_peaks_dge) <- rowData(reads_in_peaks_filtered)$peakID
colnames(reads_in_peaks_dge) <- samples_names
reads_in_peaks_dge$samples <- mutate(reads_in_peaks_dge$samples,'group'=ifelse(rownames(reads_in_peaks_dge$samples) %like% 'C','chimp','human'))

peak_width <- copy(peaks)%>%as.data.table()
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

## add sex covariate 
quant_norm_counts$sex <- c('F','M','F','F','M','M','M','F','M','F','M','F') # order follows the samples_names vector

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(outplot_dir,genome,'_',"peak_logrpkm_counts_quantnormd.pdf",sep=''),width=7,height=7)
plotDensities(log(quant_norm_counts$counts+1), group=rownames(quant_norm_counts$samples),col=samples_palette) 
dev.off()

##--------------
## Run PCA 
##--------------
run_pca <-function(counts){
    pca <- pca(t(counts$counts),ncomp=10,scale=F) 
    group = c(rep('chimp',6),rep('human',6))
    return(pca)
}
pca <- run_pca(quant_norm_counts)

pdf(paste(outplot_dir,genome,'_',"peak_pca.pdf",sep=''),width = 7, height = 7)
plotIndiv(pca,comp = c(1,2),legend=T,ind.names = T)
dev.off()

## check whether PCs are associated with things such as lib.sizes, FRiP, sex, other than species 
## read frip scores 
frip_files <- list.files(paste('../',genome,'/output/PeakCalling/qc/',sep=''),recursive = F,full.names=T,pattern="^H.*_default.frip.txt$|C.*_default.frip.txt$")
frip_scores <- lapply(frip_files,function(x)fread(x,header=F,select='V7',col.names='frip'))%>%rbindlist()
frip_scores <- frip_scores[,samples:=samples_names]

get_pc_associations <- function(pca_results){
     pcs <- copy(pca_results$variates$X)%>%as.data.frame()
     all_pcs_associations <- data.frame()
     pca_ncomp <- pca_results$ncomp
 
     get_neglog10pval <- function(x){
         neglog10pval <- -log10(anova(lm(pcs[,i] ~x))$Pr[1]) 
         return(neglog10pval)
     }
 
     for (i in 1:pca_ncomp){
         species_assoc <- get_neglog10pval(quant_norm_counts$samples$group)
         sex_assoc <- get_neglog10pval(quant_norm_counts$sex)
         libsize_assoc <- get_neglog10pval(quant_norm_counts$samples$lib.size) 
         frip_assoc <- get_neglog10pval(frip_scores$frip) 
         single_pc_assoc <- c(paste('pc',i,sep=''), species_assoc, sex_assoc, libsize_assoc, frip_assoc)
         all_pcs_associations <- rbind(all_pcs_associations, single_pc_assoc)
     }
 
     names(all_pcs_associations) <- c("PC", "species", "sex", "libsize", "frip")
    return(all_pcs_associations)
}

pcs_associations <- get_pc_associations(pca)

## plot PCs associations
top_pc_associations <- copy(pcs_associations)%>%as.data.table()
top_pc_associations <- top_pc_associations[
    ,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 2:5
    ][
        , .(Max = do.call(max, .SD)), .SDcols = 2:5, .(PC,maximum_column)
        ][
            ,Max:=round(as.numeric(Max),2)
            ]
top_pc_associations <- top_pc_associations[,PC:=factor(PC,levels=top_pc_associations$PC)]

pdf(paste(outplot_dir,genome,'_',"top_pc_association.pdf",sep=''),width=7,height=7)
ggplot(top_pc_associations,aes(x=PC,y=Max,fill=maximum_column))+
geom_histogram(stat='identity', position=position_dodge(width=0.5))+
xlab('PC')+
ylab('-log10 pval of most associated factor')+
geom_text(aes(label=maximum_column), position='dodge', vjust=-0.25,size=5)+
geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.5)+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position='bottom')
dev.off()

##------------------------------
## heatmap correlation counts
##------------------------------
plot_heatmap <- function(counts){
    Heatmap(
        cor(counts$counts,method='pearson'), 
        name = "Pearson corr", col=viridis(10)
        ) 
}

pdf(paste(outplot_dir,genome,'_',"peak_heatmap.pdf",sep=''),width=7,height=7)
plot_heatmap(quant_norm_counts)
dev.off()

##------------
## DA testing
##------------
## find best fitting model
model_species_only <- model.matrix(~0+quant_norm_counts$samples$group)  
model_species_frip <- model.matrix(~0+quant_norm_counts$samples$group+frip_scores$frip)  
model_species_sex <- model.matrix(~0+quant_norm_counts$samples$group+quant_norm_counts$sex)  
model_species_sex_frip <- model.matrix(~0+quant_norm_counts$samples$group+quant_norm_counts$sex+frip_scores$frip)  

designlist <- list(
  species_only = model_species_only,
  species_frip = model_species_frip,
  species_sex = model_species_sex,
  species_sex_frip = model_species_sex_frip
#   Full=cbind(Int=1,A=A,B=B,AB=A*B)
)

compare_models <- selectModel(quant_norm_counts$counts,designlist)

pdf(paste(outplot_dir,genome,'_',"compare_model_design.pdf",sep=''),width = 7, height = 7)
aic_score =  as.data.table(table(compare_models$pref))
colnames(aic_score) = c('model_design','aic_score')
ggplot(aic_score,aes(x=reorder(model_design,-aic_score),y=aic_score)) +
geom_histogram(stat='identity')+
xlab('model design')+ylab('AIC score')+
theme_bw()
dev.off()

## get best model with lowest aic score
colnames(model_species_sex_frip) <- c("chimp", "human",'sex_M','frip')

y <- estimateDisp(quant_norm_counts, model_species_sex_frip) 
fit <- glmQLFit(y, model_species_sex_frip, robust=TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp,levels=model_species_sex_frip))

## get df with final results_df
peaks_df <- copy(peaks)%>%as.data.table()

final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'peakID'= rownames(final_results))%>%as.data.table()

final_results <-  final_results[
        ,FDR:=p.adjust(final_results$PValue,method = 'fdr')
        ][
            ,DA:=as.factor(ifelse(FDR<=0.05, 'da','non_da'))
            ][
                peaks_df,on='peakID',nomatch=0
            ][
                ,c('width','strand'):=NULL
                ][
                    ,peak_species := ifelse(DA == 'non_da' ,'non_da',ifelse(DA=='da'& logFC <0,'da_chimp','da_human'))
] 

## QCs
## check distribution Pvalues
pdf(paste(outplot_dir,genome,'_',"da_test_distribution_raw_pvals.pdf",sep=''),width = 7, height = 7)
ggplot(final_results,aes(x=round(PValue,4))) +geom_histogram(binwidth=0.01)+
xlab('raw pvalue')+
theme_bw()
dev.off()

## check distribution peak lengths between da and nonda peaks
pdf(paste(outplot_dir,genome,'_',"da_test_peak_sizes.pdf",sep=''),width = 7, height = 7)
df <- copy(final_results)[,width:=end-start]
ggplot(df,aes(x=width,fill=peak_species,col=peak_species)) +
geom_density(alpha=0.5) +
geom_vline(xintercept=150, linetype='dashed', color='black', size=0.5)+
geom_vline(xintercept=300, linetype='dashed', color='black', size=0.5)+
theme_bw()
dev.off()

write.table(final_results,paste(outfile_dir,'da_results.txt',sep=''),sep='\t',col.names=T,row.names = F,quote=F)
##--------------------
## Generate MA plot
##--------------------
ma_plot = function(df){
    plot <- ggplot(df,aes(x = logCPM, y = logFC, col = peak_species)) + 
        geom_point() + 
        # scale_color_manual(values = da_palette) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess",se = F) + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)+
        theme_bw()
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
    # scale_color_manual(values = da_palette) + 
    geom_hline(yintercept=1.3, linetype='dashed', color='black', size=1)+
    geom_vline(xintercept=0, linetype='dashed', color='black', size=1)+
    theme_bw()
    
    return(plot)
}

pdf(paste(outplot_dir,genome,'_',"peak_volcano_plot.pdf",sep=''),width = 7, height = 7)
volcano_plot(final_results)
dev.off()
