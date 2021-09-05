## script to perform DA analysis using csaw
## here i use RPKM 

library(data.table)
library(magrittr)
library(dplyr)
library(edgeR)
library(GenomicRanges)
library(csaw)
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(liftOver)
library(rtracklayer)
library(mixOmics)

options(width=150)

## directories
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

bamDir  = "/output/Post_alignment/Files"
result_dir = './post_processing_analyses/output'
peakDir = '/output/PeakCalling/Files/'
plot_dir = paste0(result_dir,'/plots/DA/',sep='')
peak_outdir = paste0(result_dir,'/DA/peaks/',sep='')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

## read peak
read_peaks=function(species){
    peaks=fread(
        paste(species,peakDir,'merged_macs2_default_peaks_filtered_sorted.narrowPeak.gz',sep=''),
        sep='\t',header=F,select=c(1:4),col.names=c(range_keys,'peakID'))%>%setorderv(range_keys,1)%>%
        makeGRangesFromDataFrame(keep.extra.columns=T)
    return(peaks)
}

pantro5_peaks=read_peaks('panTro5')
hg38_peaks=read_peaks('hg38')

## list bam
get_bams=function(species){
    bams=list.files(paste0(species,bamDir), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}
pantro5_bams=get_bams('panTro5')
hg38_bams=get_bams('hg38')

standard.chr <- paste0("chr", c(1:23, "X", "Y",'2A','2B')) # only use standard chromosomes
param <- readParam(restrict= standard.chr )

## get overlaps
get_overlaps <- function(bams,peaks){
    overlaps =  regionCounts(bams,peaks,param=param)
    
    samples = sub("\\_.*", "", sub('.*\\/', '', bams))
    samples = gsub("*.bam", "", samples)
    
    overlaps_DGE <- asDGEList(overlaps)

    peak_width = copy(peaks)%>%as.data.table()%>%dplyr::pull('width')
    overlaps_DGE$counts = rpkm(overlaps_DGE$counts,peak_width) ## this passed manual check

    rownames(overlaps_DGE) = peaks$peakID
    colnames(overlaps_DGE) = samples
    overlaps_DGE$samples=mutate(overlaps_DGE$samples,'group'=ifelse(rownames(overlaps_DGE$samples) %like% 'C','chimp','human'))


    return(overlaps_DGE)
} 

hg38_counts <- get_overlaps(hg38_bams, hg38_peaks)
pantro5_counts <- get_overlaps(pantro5_bams, pantro5_peaks)

##----------
## Some QCs
##----------
density_plot <- function(counts){
    plotDensities(log(counts$counts+1), group=rownames(counts$samples),col=brewer.pal(12, 'Set3')) 
}
## plot density of log-RPKM for each sample and decide cutoff for min number of reads in peaks
pdf(paste(plot_dir,"hg38_logRPKM_read_densities.pdf",sep=''),width=7,height=7)
density_plot(hg38_counts)
dev.off()

pdf(paste(plot_dir,"panTro5_logRPKM_read_densities.pdf",sep=''),width=7,height=7)
density_plot(pantro5_counts)
dev.off()

## filter counts
## by retaining only peaks with log-RPKM > -xxx
filter_counts <- function(counts){
    keep_exprs <- filterByExpr(counts$counts,group=counts$samples$group)
    filtered_counts <- counts[keep_exprs,, keep.lib.sizes=FALSE]
    return(filtered_counts)
}

hg38_filtered_counts <- filter_counts(hg38_counts)
panTro5_filtered_counts <- filter_counts(pantro5_counts)

## plot density of log-RPKM for each sample
pdf(paste(plot_dir,"hg38_logRPKM_filtered_read_densities.pdf",sep=''),width=7,height=7)
density_plot(hg38_filtered_counts)
dev.off()

pdf(paste(plot_dir,"panTro5_logRPKM_filtered_read_densities.pdf",sep=''),width=7,height=7)
density_plot(panTro5_filtered_counts)
dev.off()

##------------------
## Normalisation
##------------------
# ## because there is a bit too much varibility across samples we are using quantile normalisation (voom) rather than TMM
# ## look at it yourself
# tmm_normd_filtered_counts <- calcNormFactors(filtered_counts,method='TMM')

# ## plot the density of log-RPKM for the tmm norm filtered counts for each sample
# pdf(paste(plot_dir,"hg38_logRPKM_tmm_normalised_filtered.pdf",sep=''),width=7,height=7)
# plotDensities(log(tmm_normd_filtered_counts$counts+1), group=rownames(tmm_normd_filtered_counts$samples),col=brewer.pal(12, 'Set3')) 
# dev.off()

## quantile norm with voom
quant_norm <- function(filtered_counts){
    quant_norm_counts = copy(filtered_counts)
    quant_norm_counts$counts =normalizeQuantiles(quant_norm_counts$counts, ties=TRUE)
    return(quant_norm_counts)
}

hg38_quant_norm_counts = quant_norm(hg38_filtered_counts)
pantro5_quant_norm_counts = quant_norm(panTro5_filtered_counts)


## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(plot_dir,"hg38_logRPKM_quant_normalised_filtered.pdf",sep=''),width=7,height=7)
density_plot(hg38_quant_norm_counts)
dev.off()

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(plot_dir,"panTro5_logRPKM_quant_normalised_filtered.pdf",sep=''),width=7,height=7)
density_plot(pantro5_quant_norm_counts)
dev.off()

##-------------------
## Run PCA 
##-------------------
run_pca <-function(quant_norm_counts){
    pca <- pca(t(quant_norm_counts$counts),ncomp=5,scale=T) 
    group = c(rep('chimpanzee',6),rep('human',6))
    return(pca)
}
hg38_pca <- run_pca(hg38_quant_norm_counts)
pantro5_pca <- run_pca(pantro5_quant_norm_counts)

## plot PCA results
pdf(paste0(plot_dir,'hg38_PCA.pdf'),width = 7, height = 7)
plotIndiv(hg38_pca,group=group,comp = c(1,2),legend=T,ind.names = T)
dev.off()

pdf(paste0(plot_dir,'panTro5_PCA.pdf'),width = 7, height = 7)
plotIndiv(pantro5_pca,group=group,comp = c(1,2),legend=T,ind.names = T)
dev.off()


## calculate correaltion across samples for reads in peaks
## separate species-specific from common peaks
## but also make a single comprehensive heatmap
library(ComplexHeatmap)
library(viridis)

plot_heatmap <- function(quant_norm_counts){
    Heatmap(
        cor(quant_norm_counts$counts,method='pearson'), 
        name = "Pearson corr", col=viridis(100)
        ) 
}
pdf(paste(plot_dir,"hg38_PearsonCor_read_in_all_peaks.pdf",sep=''),width=7,height=7)
plot_heatmap(hg38_quant_norm_counts)
dev.off()

pdf(paste(plot_dir,"panTro5_PearsonCor_read_in_all_peaks.pdf",sep=''),width=7,height=7)
plot_heatmap(pantro5_quant_norm_counts)
dev.off()

##------------------------------------------
## DIFFERENTIAL ACCESSIBILITY ANALYSIS
##------------------------------------------
da_testing <- function(quant_norm_counts){
  design_matrix <- model.matrix(~0+group, data=quant_norm_counts$samples)  
colnames(design_matrix) <- c("chimp", "human")

y <- estimateDisp(quant_norm_counts, design_matrix) 
fit <- glmQLFit(y, design_matrix, robust=TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp, levels=design_matrix))
fdr_results <- getBestTest(rownames(results$table),results$table)
fdr_results <- data.table(peakID=rownames(fdr_results),FDR=fdr_results$FDR,direction=fdr_results$direction)

final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'peakID'=rownames(final_results))%>%as.data.table()
final_results <-  final_results[
        fdr_results,on='peakID',nomatch=0
        ][
            ,significant:=ifelse(FDR<0.05,'significant','non_significant')
            ]  
return(final_results)
}

hg38_da <- da_testing(hg38_quant_norm_counts)
pantro5_da <- da_testing(pantro5_quant_norm_counts)

# ## setup design matrix and fit model to data
# design_matrix <- model.matrix(~0+group, data=filtered_counts_quant_norm$samples)  ## 0+ in the formula is an instruction not to include an intercept column and instead to include a column for each group
# colnames(design_matrix) <- c("chimp", "human")

# ## model the dispersion/variance in accessibility assuming most regions arent DA  
# ## and stabilize dispersion estimates with empirical bayes
# ## edgeR determines the common dispersion across features.
# ## It determines the common variance across features and uses this to model the mean-variance relationship.
# ## This is then used to calculate a dispersion estimate per feature which is the necessary to test for DE/DA

# ## The variance for each feature measures the degree of inter-library variation for that feature
# ## edgeR calculates the common dispersion across feature which, in turn, gives an idea of the overall variability across the genome for the dataset

# y <- estimateDisp(filtered_counts_quant_norm, design_matrix) 

# ## Fits a quasi-likelihood negative binomial generalized log-linear model to data
# fit <- glmQLFit(y, design_matrix, robust=TRUE)

# ## test for DA
# results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp, levels=design_matrix))

# ## correct for multiple hypothesis
# fdr_results = getBestTest(rownames(results$table),results$table)
# fdr_results = data.table(peakID=rownames(fdr_results),FDR=fdr_results$FDR,direction=fdr_results$direction)

# ## Concatenate all relevant statistical data and peak info (e.g. regions and species) into a single data.table
# final_results = data.frame(results$table)
# final_results = dplyr::mutate(final_results,'peakID'=rownames(final_results))%>%as.data.table()
# final_results = final_results[
#     # genomic_regions,on='peakID',nomatch=0
#     # ][
#         fdr_results,on='peakID',nomatch=0
#         ][
#             ,significant:=ifelse(FDR<0.05,'significant','non_significant')
# ]
# # [
# #                 ,species:=ifelse(peakID %in% rownames(common_counts_combined),'common',
# #                             ifelse(grepl("c",peakID),'chimp','human'))
# #             ][
# #                 , c("human_peakID", "chimp_peakID") := tstrsplit(peakID, ".", fixed=TRUE)
# #                 ][
# #                     ,human_peakID := ifelse(species=='common',sub(".*?_", "", human_peakID),peakID)
# #                     ][
# #                         ,chimp_peakID := ifelse(species=='common',sub(".*?_", "", chimp_peakID),peakID)
# # ]

##--------------------
## Generate MA plot
##--------------------
ma_plot = function(df){
    plot <- ggplot(df,
        aes(x = logCPM, y = logFC, col = factor(significant, levels=c("non_significant", "significant")))) + 
        geom_point() + scale_color_manual(values = c("black", "red")) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)
        return(plot)
}

pdf(paste0(plot_dir,'hg38_rpkm_quant_norm_normalisation_plot.pdf'),width = 7, height = 7)
ma_plot(hg38_da)
dev.off()

pdf(paste0(plot_dir,'panTro5_rpkm_quant_norm_normalisation_plot.pdf'),width = 7, height = 7)
ma_plot(pantro5_da)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x=logFC,y=-log10(FDR)))+
    geom_hline(yintercept=1.3, linetype='dashed', color='black', size=1)
    # +facet_wrap(species~.,ncol=3)
    
    return(plot)
}

pdf(paste0(plot_dir,'hg38_rpkm_quant_norm_volcano_humanVSchimp.pdf'),width = 10, height = 7)
volcano_plot(hg38_da)
dev.off()

pdf(paste0(plot_dir,'panTro5_rpkm_quant_norm_volcano_humanVSchimp.pdf'),width = 10, height = 7)
volcano_plot(pantro5_da)
dev.off()

# ##---------------------------------
# ## plot mean-var trend with voom
# ##----------------------------------
# v_no_normd <- voomLmFit(
#     filtered_counts, design_matrix, 
#     plot=F,save.plot = T,
#     normalize.method='none'
# )

# pdf(paste0(plot_dir,'voom_mean_var_trend_no_normd.pdf'),width = 7, height = 7)
# plot(v_no_normd$voom.xy$x, v_no_normd$voom.xy$y, xlab = v_no_normd$voom.xy$xlab, ylab = v_no_normd$voom.xy$ylab)
# title("voom: Mean-variance trend")
# lines(v_no_normd$voom.line, col = "red")
# dev.off()

# v_quantile_normd <- voomLmFit(
#     filtered_counts, design_matrix, 
#     plot=F,save.plot = T,
#     normalize.method='quantile'
# )

# pdf(paste0(plot_dir,'voom_mean_var_trend_quant_normd.pdf'),width = 7, height = 7)
# plot(v_quantile_normd$voom.xy$x, v_quantile_normd$voom.xy$y, xlab = v_quantile_normd$voom.xy$xlab, ylab = v_quantile_normd$voom.xy$ylab)
# title("voom: Mean-variance trend")
# lines(v_quantile_normd$voom.line, col = "red")
# dev.off()

##------------------
## write results
##------------------
## now create a single file with all the relevant info, i.e.
## 1) da results
## 2) peakID
## 3) genomic positions
## 4) species 
 
human_genomic_regions = copy(as.data.table(human_and_common_peaks))[,peakID:= gsub("\\..*","",peakID)][,c('width','strand'):=NULL]
chimp_genomic_regions = copy(as.data.table(chimp_and_common_peaks))[,peakID:= gsub("\\..*","",peakID)][,c('width','strand'):=NULL]

common_peaks_results = copy(final_results)[species=='common']%>%
    inner_join(human_genomic_regions,by=c('human_peakID'='peakID','species'))%>%
        inner_join(chimp_genomic_regions,by=c('chimp_peakID'='peakID','species'))%>%
            setnames(
                old = c(paste(range_keys,'x',sep='.'),paste(range_keys,'y',sep='.')), 
                new = c(paste('human',range_keys,sep='_'),paste('chimp',range_keys,sep='_')))%>%
                    dplyr::select(c(contains('human'),contains('chimp'),everything(),-'peakID'))

human_specific_results  = copy(final_results)[
    species=='human'
    ][
        ,peakID:= gsub("\\..*","",human_peakID)
        ][
            ,c('human_peakID','chimp_peakID') := NULL
        ]%>%inner_join(human_genomic_regions,by=c('peakID','species'))%>%
        dplyr::select(c(all_of(range_keys),everything()))

chimp_specific_results  = copy(final_results)[
    species=='chimp'
    ][
        ,peakID:= gsub("\\..*","",chimp_peakID)
        ][
            ,c('human_peakID','chimp_peakID') := NULL
        ]%>%inner_join(chimp_genomic_regions,by=c('peakID','species'))%>%
        dplyr::select(c(all_of(range_keys),everything()))

all_results = list(chimp_specific_results,common_peaks_results,human_specific_results)
names(all_results) = c('chimp','common','human')

mapply(write.table, all_results, file = paste0(peak_outdir,names(all_results),'_da_results',sep='.txt'),sep='\t',col.names=T,row.names=F,quote=F)

## read common regions and add the DA results to them
## this will make your life easier with dowstream analyses
common_regions = fread(paste(peakDir,'common_regions_hg38.bed',sep=''),sep='\t',header=T)

common_regions_da = copy(common_regions)[
    all_results[[2]],on=c('human_peakID','chimp_peakID'),nomatch=0
    ][
        ,peakID:=paste(paste('C_',chimp_peakID,sep=''),paste('H_',human_peakID,sep=''),sep='.')
]%>%dplyr::select(-c(contains('chimp'),contains('human')))

write.table(common_regions_da,paste0(peak_outdir,'common_regions.txt',sep=''),sep='\t',col.names=T,row.names=F,quote=F)

## Plot proportion of DA/non-DA peaks
prop_peaks = function(x){
    x=x[
        ,numb_peaktype :=.N,by=.(significant)
        ][
            ,tot_peaks :=.N
            ][
                ,prop_peaktype:=numb_peaktype/tot_peaks
                ]
    return(x)
}

prop_da_peaks = copy(all_results)%>%
lapply(
    function(x)x=x[
        ,numb_peaktype :=.N,by=.(significant)
        ][
            ,tot_peaks :=.N
            ][
                ,prop_peaktype:=numb_peaktype/tot_peaks
                ][
                    ,c('species','prop_peaktype','significant','numb_peaktype')
                    ]%>%unique()
)%>%rbindlist()

pdf(paste(plot_dir,'prop_da_nonda_peaks_species.pdf',sep=''),width=7,height=7)
ggplot(prop_da_peaks,aes(x=species,y=prop_peaktype,fill=significant))+
    geom_bar(position="stack", stat="identity")
dev.off()
