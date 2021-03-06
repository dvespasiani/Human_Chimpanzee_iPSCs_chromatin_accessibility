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

outplot_dir <- create_dir(plot_dir,'DA')

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

numb_all_peaks <- length(peaks)

##---------------------
## QC & Normalisation
##---------------------
## Use RPKM to account for differences in peak sizes
reads_in_peaks <- regionCounts(bams, peaks,param=param)%>%asDGEList()
peak_width <- copy(peaks)%>%as.data.table()%>%dplyr::pull('width')
reads_in_peaks$counts <- rpkm(reads_in_peaks$counts,peak_width)

rownames(reads_in_peaks) <- as.data.table(copy(peaks))$peakID
colnames(reads_in_peaks) <- samples_names
reads_in_peaks$samples <- mutate(reads_in_peaks$samples,'group'=ifelse(rownames(reads_in_peaks$samples) %like% 'C','chimp','human'))

pdf(paste(outplot_dir,"density_reads_in_peaks_rpkm.pdf",sep=''),width=7,height=7)
plotDensities(log(reads_in_peaks$counts+1), group=rownames(reads_in_peaks$samples),col=samples_palette) 
dev.off()

## filter low reads 
keep <- filterByExpr(reads_in_peaks, group=reads_in_peaks$group,min.total.count = 3,min.count=3)
reads_in_peaks_filtered <- copy(reads_in_peaks)[keep,, keep.lib.sizes=FALSE]
dim(reads_in_peaks_filtered)

## quantile normalise read counts
reads_in_peaks_filtered_tmm <- calcNormFactors(copy(reads_in_peaks_filtered), method = "TMM")
quant_norm_counts <- copy(reads_in_peaks_filtered_tmm)
quant_norm_counts$counts <- normalizeQuantiles(reads_in_peaks_filtered_tmm$counts, ties=F)

## add sex covariate 
quant_norm_counts$sex <- c('F','M','F','F','M','M','M','F','M','F','M','F') # order follows the samples_names vector

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(outplot_dir,"density_filtered_reads_in_peaks_rpkm_quantnorm.pdf",sep=''),width=7,height=7)
plotDensities(log(quant_norm_counts$counts+1), group=rownames(quant_norm_counts$samples),col=samples_palette) 
dev.off()

##------------------------------
## Clustering: PCA + Heatmap 
##------------------------------
pca_rpkm_quantnorm_counts <- pca(t(quant_norm_counts$counts),ncomp=5,scale=F) 

pca_plot <- function(pca_results,components,cols){
    prop_exp_var = pca_results$prop_expl_var$X[c(components)]

    label_plot <-function(comp){
        label_comp <-  paste(paste(names(comp),paste('(',round(comp*100,2),'%',')',sep=''),sep= ' '))
        return(label_comp)
    }

    df <- copy(pca_results$variates$X)%>%as.data.frame()
    rownames_df <- rownames(df)
    df <- df%>%dplyr::select(c(all_of(c(components))))%>%mutate("names"=rownames_df)%>%as.data.table()
    df <- df[,group:=ifelse(names %like% 'C','chimp','human')]%>%setnames(old=c(1:2),new=c('x','y'))

    plot <-ggplot(df,aes(x=x,y=y,col=group))+
    geom_point(size=3)+
    scale_color_manual(values = cols) +
    xlab(label_plot(prop_exp_var[1]))+
    ylab(label_plot(prop_exp_var[2]))+
    theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
  return(plot)
}

pdf(paste(outplot_dir,"pca_rpkm_quantnorm_counts.pdf",sep=''),width = 7, height = 7)
pca_plot(pca_rpkm_quantnorm_counts,c(1:2),c('#14213d','#fca311'))
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

pcs_associations <- get_pc_associations(pca_rpkm_quantnorm_counts)

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

pdf(paste(outplot_dir,"pca_top_pc_association.pdf",sep=''),width=7,height=7)
ggplot(top_pc_associations,aes(x=PC,y=Max,fill=maximum_column))+
geom_histogram(stat='identity', position=position_dodge(width=0.5))+
xlab('PC')+
ylab('-log10 pval of most associated factor')+
geom_text(aes(label=maximum_column), position='dodge', vjust=-0.25,size=5)+
geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.5)+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position='bottom')
dev.off()

## Heatmap 
pdf(paste(outplot_dir,"heatmap_quantnorm_rpkm.pdf",sep=''),width=7,height=7)
Heatmap(cor(quant_norm_counts$counts,method='pearson'), name = "Pearson corr", col=viridis(10)) 
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

pdf(paste(outplot_dir,"compare_model_design.pdf",sep=''),width = 7, height = 7)
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
pdf(paste(outplot_dir,"da_test_distribution_raw_pvals.pdf",sep=''),width = 7, height = 7)
ggplot(final_results,aes(x=round(PValue,4))) +geom_histogram(binwidth=0.01)+
xlab('raw pvalue')+
theme_bw()
dev.off()

## check distribution peak lengths between da and nonda peaks
pdf(paste(outplot_dir,"da_test_peak_sizes.pdf",sep=''),width = 7, height = 7)
df <- copy(final_results)[,width:=end-start]
ggplot(df,aes(x=width,fill=peak_species,col=peak_species)) +
geom_density(alpha=0.5) +
geom_vline(xintercept=150, linetype='dashed', color='black', size=0.5)+
geom_vline(xintercept=300, linetype='dashed', color='black', size=0.5)+
theme_bw()
dev.off()

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

pdf(paste(outplot_dir,"da_test_ma_plot.pdf",sep=''),width = 7, height = 7)
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

pdf(paste(outplot_dir,"da_test_volcano_plot.pdf",sep=''),width = 7, height = 7)
volcano_plot(final_results)
dev.off()

## export DA resutls
write.table(final_results,paste(da_dir,'da_results.txt',sep=''),sep='\t',col.names=T,row.names = F,quote=F)
