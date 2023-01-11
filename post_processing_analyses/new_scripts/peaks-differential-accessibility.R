## script for DA testing of peaks between humans and chimps

library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(edgeR)
library(GenomicRanges)
library(csaw)
library(mixOmics)
library(preprocessCore)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility')

source('./post_processing_analyses/scripts/utils.R')

peak_inputdir = paste(basedir,'post_processing_analyses/output/files/consensusPeaks/',sep='')
outplot_dir = create_dir(paste(plot_dir,sep=''),'DA')

##---------------------
## reusable functions
##---------------------
make_dt <- function(counts){
    dt <-copy(as.data.table(counts))
    long_dt <- melt(dt,measure.vars=names(dt))
    long_dt <- long_dt[,group:=ifelse(variable %like% 'C','chimp','human')]
    return(long_dt)
}

plot_boxes <- function(dt,ylab){
    p <- ggplot(dt,aes(x=variable,y=value,fill=variable))+
    geom_boxplot()+
    scale_fill_manual(values=samples_palette)+
    xlab(' ')+ylab(ylab)+
    theme_classic()+
    theme(
        legend.position = "bottom",
        axis.ticks.x =element_blank()
    )
  return(p)
}

plot_densities <-function(dt){
    p <- ggplot(dt,aes(x=value,col=variable))+
    geom_density()+
    scale_color_manual(values=samples_palette)+
    xlab('log2 CPM')+ylab('Density')+
    theme_classic()+
    theme(
        legend.position = "bottom",
        axis.ticks.x =element_blank()
    )
    return(p)
}


## Read the mapped consensus peaks
## for this and all downstream analyses use peaks in hg38 coords. 
## so get the panTro5 peaks mapped to hg38 coordinates and then combine the 2 data tables 
get_peaks = function(file){
    f = fread(paste(peak_inputdir,file,sep=''),header=T)[seqnames %in% standard_chr]
    setkeyv(f,range_keys)
    f = f%>%setorderv(c(range_keys,'peakID'),1)
    return(f)
}

human_cons_peaks = get_peaks('human_consensus_peaks_liftForward_liftBack.bed')# hg38 coords
human_map_peaks_panTro5 = get_peaks('human_consensus_peaks_liftForward.bed') # pantro5 coords
human_map_peaks_panTro5 <- human_map_peaks_panTro5[peakID %in% human_cons_peaks$peakID]
# human_map_peaks_hg38 = copy(human_cons_peaks)[peakID %in% human_map_peaks_panTro5$peakID]

chimp_cons_peaks = get_peaks('chimp_consensus_peaks_liftForward_liftBack.bed') # pantro5 coords
chimp_map_peaks_hg38 = get_peaks('chimp_consensus_peaks_liftForward.bed') # hg38 coords
chimp_map_peaks_hg38 <- chimp_map_peaks_hg38[peakID %in% chimp_cons_peaks$peakID]
# chimp_map_peaks_panTro5 = copy(chimp_cons_peaks)[peakID %in% chimp_map_peaks_hg38$peakID]

## get common peaks by looking at overlap betweeen genomic regions
## and create 2 common peaks dt (1 x species) with the common genomic regions in the relative species genomic coordinates
common_peaks_hg38 = foverlaps(chimp_map_peaks_hg38,human_cons_peaks,type='any')%>%na.omit()%>%unique()%>%setnames(old=c('peakID','i.peakID'),new=c('human_peakID','chimp_peakID'))

unique_common_peaks = copy(common_peaks_hg38)[,c('human_peakID','chimp_peakID','i.support','support')]%>%unique()
unique_common_peaks = unique_common_peaks[
    ,duphuman:=.N,by=.(chimp_peakID)
    ][
        ,dupchimp:=.N,by=.(human_peakID)
        ][
            duphuman == 1 & dupchimp==1
            ][
                i.support > 1 & support >1
                ][
                    ,peaktype:='common'
]

common_peaks_hg38 = copy(human_cons_peaks)[peakID %in% unique_common_peaks$human_peakID][,peaktype:='common']
common_peaks_panTro5 = copy(chimp_cons_peaks)[peakID %in% unique_common_peaks$chimp_peakID][,peaktype:='common']

## now get species-specific peaks in their relative species coordinates
human_sp_hg38 = copy(human_map_peaks_hg38)[! peakID %in% common_peaks_hg38$peakID][,peaktype:='human_sp']
chimp_sp_panTro5 = copy(chimp_map_peaks_panTro5)[! peakID %in% common_peaks_panTro5$peakID][,peaktype:='chimp_sp']

## repeat this step but with the mixed species-genomic coordinates files
human_sp_panTro5 = copy(human_map_peaks_panTro5)[ peakID %in% human_sp_hg38$peakID][,peaktype:='human_sp']
chimp_sp_hg38 = copy(chimp_map_peaks_hg38)[ peakID %in% chimp_peaks_panTro5$peakID][,peaktype:='chimp_sp']

## now make 2 final datatables with all common and sp peaks in the 2 genomic coordinates (you need this for DA)
make_final_peaks <- function(peaks){
    return <- copy(peaks)
    return <- return[,peakID:=paste('peak_',1:nrow(return),sep='')]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
    return(return)
}

final_hg38_peaks = make_final_peaks(rbind(common_peaks_hg38,human_sp_hg38,chimp_sp_hg38))
final_panTro5_peaks = make_final_peaks(rbind(common_peaks_panTro5,chimp_sp_panTro5,human_sp_panTro5))

# pdf(paste(outplot_dir,"qc-dist-peaksizes-pre-da.pdf",sep=''),width=7,height=5)
# options(scipen=999)
# df <- as.data.table(copy(final_peaks))
# ggplot(df,aes(x=width,fill=peaktype))+
# geom_density(alpha=0.5)+
# scale_fill_manual(values=c('#83c5be','#e29578','#006d77'))+
# xlab('Peak size')+ ylab('Density')+
# theme_classic()+
# theme()
# dev.off()

##---------------------------------------
## Differential Accessibility
##---------------------------------------
param <- readParam(pe = "both",restrict=standard_chr,max.frag=1000)

get_bamReads = function(organimsDir,pattern){
    bamReads = list.files(paste0(organimsDir,bamDir), 
        recursive = T,full.names = T,pattern=pattern)
    return(bamReads)
}

hg38_bams = get_bamReads('./hg38/',"^H.*_tn5_shifted_sorted.bam$")
pantro5_bams = get_bamReads('./panTro5/',"^C.*_tn5_shifted_sorted.bam$")

## count number reads overlapping orthologous consensus peaks 
## and convert to RPKM
raw_counts_human <- regionCounts(hg38_bams, final_hg38_peaks,param=param)%>%asDGEList()
colnames(raw_counts_human) = samples_names[7:12]

human_peaks_width <- copy(final_hg38_peaks@ranges@width)

rpkm_counts_human <- copy(raw_counts_human)
rpkm_counts_human$counts <- rpkm(raw_counts_human$counts,gene.length=human_peaks_width)

raw_counts_chimp <- regionCounts(pantro5_bams, final_panTro5_peaks,param=param)%>%asDGEList()
colnames(raw_counts_chimp) = samples_names[1:6]

chimp_peaks_width <-copy(final_panTro5_peaks@ranges@width)
rpkm_counts_chimp <- copy(raw_counts_chimp)
rpkm_counts_chimp$counts <- rpkm(raw_counts_chimp$counts,gene.length=chimp_peaks_width)


## merge matrices 
merged_raw_counts = cbind(raw_counts_chimp,raw_counts_human)
merged_raw_counts$samples <- mutate(merged_raw_counts$samples,'group'=ifelse(rownames(merged_raw_counts$samples) %like% 'C','chimp','human'))

mat_rownames = final_hg38_peaks@elementMetadata@listData$peakID
rownames(merged_raw_counts) = mat_rownames

merged_rpkm_counts = cbind(rpkm_counts_chimp,rpkm_counts_human)
merged_rpkm_counts$samples <- mutate(merged_rpkm_counts$samples,'group'=ifelse(rownames(merged_rpkm_counts$samples) %like% 'C','chimp','human'))
rownames(merged_rpkm_counts) =mat_rownames

##---------------------
## Filtering
##---------------------
prior_counts = 0.5

## plot log2 CPM before and after filtering
pdf(paste(outplot_dir,"qc-density-reads-in-peaks.pdf",sep=''),width=7,height=7)
plot_densities(make_dt(cpm(merged_raw_counts$counts,log=T,prior.counts = prior_counts)))
dev.off()

## filter low reads
dt_list <- list()
for (i in 1:12){
    false = table(rowSums(merged_raw_counts$counts==0)==i)[1]
    true = table(rowSums(merged_raw_counts$counts==0)==i)[2]
    dt_list[[i]] <- data.table(false=false,true=true,numb_samples=i)
}
dt_list = rbindlist(dt_list)
dt_list[is.na(dt_list)]=0

keep <- filterByExpr(merged_raw_counts)
raw_counts_filtered <- copy(merged_raw_counts)[keep,, keep.lib.sizes=FALSE]
rpkm_counts_filtered <- copy(merged_rpkm_counts)[keep,, keep.lib.sizes=FALSE]
dim(merged_raw_counts)
dim(raw_counts_filtered)

pdf(paste(outplot_dir,"qc-density-reads-in_peaks-filtered.pdf",sep=''),width=7,height=7)
plot_densities(make_dt(cpm(raw_counts_filtered$counts,log=T,prior.counts = prior_counts)))
dev.off()

##------------------
## Normalisation
##------------------
# cpm_counts_filtered = copy(raw_counts_filtered)
# cpm_counts_filtered$counts = cpm(cpm_counts_filtered$counts,log=T,prior.counts= prior_counts)

pdf(paste(outplot_dir,"qc-distribution-unnormalised-filtered-counts.pdf",sep=''),width=7,height=7)
plot_boxes(make_dt(log2(rpkm_counts_filtered$counts+prior_counts)),ylab='log2 RPKM')
dev.off()

# ## tmm normalisation
# tmm_norm_counts <- calcNormFactors(copy(raw_counts_filtered),method='TMM')

# pdf(paste(outplot_dir,"distribution_tmm_normalised_filtered_counts.pdf",sep=''),width=7,height=7)
# plot_boxes(make_dt(cpm(tmm_norm_counts$counts,log=T,prior.counts = prior_counts)),ylab='log2 CPM')
# dev.off()

# ## loess
# loess_norm_counts <- copy(raw_counts_filtered)
# loess_norm_counts$counts <- normalizeBetweenArrays(loess_norm_counts$counts, method = "cyclicloess")

# pdf(paste(outplot_dir,"distribution_loess_normalised_filtered_counts.pdf",sep=''),width=7,height=7)
# plot_boxes(make_dt(cpm(loess_norm_counts$counts,log=T,prior.counts = prior_counts)),ylab='log2 CPM')
# dev.off()

## get RPKM and quantile normalise these
quantnorm_rpkm_counts <- copy(rpkm_counts_filtered)
quantnorm_rpkm_counts$counts <- normalizeBetweenArrays(quantnorm_rpkm_counts$counts, method = "quantile")

pdf(paste(outplot_dir,"qc-distribution-quantile-normalised-filtered-counts.pdf",sep=''),width=7,height=7)
plot_boxes(make_dt(log2(quantnorm_rpkm_counts$counts+prior_counts)),ylab='log2 RPKM')
dev.off()

##------------------------------
## Clustering: PCA + Heatmap 
##------------------------------
pca <- pca(t(log2(quantnorm_rpkm_counts$counts+prior_counts)),ncomp=10,scale=T) 

pca_plot <- function(pca_results,components){
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
    scale_color_manual(values = da_species_palette[-1]) +
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

pdf(paste(outplot_dir,"qc-pca-screeplot.pdf",sep=''),width = 7, height = 7)
plot(pca)
dev.off()

pdf(paste(outplot_dir,"pca-quantnorm-counts.pdf",sep=''),width = 7, height = 7)
pca_plot(pca,c(1:2))
dev.off()

## add sex covariate 
# quant_norm_counts <- DGEList(quant_norm_counts$counts)
quantnorm_rpkm_counts$sex <- c('F','M','F','F','M','M','M','F','M','F','M','F') # order follows the samples_names vector
# merged_filterd_rpkm$samples <- mutate(merged_filterd_rpkm$samples,'group'=ifelse(rownames(quant_norm_counts$samples) %like% 'C','chimp','human'))

## check whether PCs are associated with things such as lib.sizes, FRiP, sex, other than species 
## ADD frip for each sample
# frip_files <- list.files(paste(genome,'/output/PeakCalling/qc',sep=''),recursive = F,full.names=T,pattern="^H.*_default.frip.txt$|C.*_default.frip.txt$")
human_frips <- list.files('./hg38/output/PeakCalling/qc',recursive = F,full.names=T,pattern="^H.*_default.frip.txt$")
chimp_frips <- list.files('./panTro5/output/PeakCalling/qc',recursive = F,full.names=T,pattern="^C.*_default.frip.txt$")

frip_files = c(human_frips,chimp_frips)
frip_scores <- lapply(frip_files,function(x)fread(x,header=F,select='V7',col.names='frip'))
names(frip_scores) = gsub("\\_.*","",gsub(".*/","",frip_files))
frip_scores = Map(mutate,frip_scores,samples=names(frip_scores))%>%rbindlist()

## add FRiP covariate
## NB: order frip scores matching colnames count matrix
quantnorm_rpkm_counts$frip <- c(frip_scores$frip[order(match(frip_scores$samples,colnames(quantnorm_rpkm_counts$counts)))])

get_pc_associations <- function(pca_results){
     pcs <- copy(pca_results$variates$X)%>%as.data.frame()
     all_pcs_associations <- data.frame()
     pca_ncomp <- pca_results$ncomp
 
     get_neglog10pval <- function(x){
         neglog10pval <- -log10(anova(lm(pcs[,i] ~x))$Pr[1]) 
         return(neglog10pval)
     }
 
     for (i in 1:pca_ncomp){
         species_assoc <- get_neglog10pval(quantnorm_rpkm_counts$samples$group)
         sex_assoc <- get_neglog10pval(quantnorm_rpkm_counts$sex)
         libsize_assoc <- get_neglog10pval(quantnorm_rpkm_counts$samples$lib.size) 
         frip_assoc <- get_neglog10pval(frip_scores$frip) 
         single_pc_assoc <- c(paste('pc',i,sep=''), species_assoc, sex_assoc, libsize_assoc, frip_assoc)
         all_pcs_associations <- rbind(all_pcs_associations, single_pc_assoc)
     }
 
     names(all_pcs_associations) <- c("PC", "species", "sex", "libsize", "frip")
    return(all_pcs_associations)
}

pcs_associations <- get_pc_associations(pca)

## plot PCs associations
top_pc_associations <- copy(pcs_associations)%>%as.data.table()%>%melt(id.vars='PC')
top_pc_associations <- top_pc_associations[,.SD[which.max(value)], by=.(PC)][,value:=round(as.numeric(value),2)]
top_pc_associations <- top_pc_associations[,PC:=factor(PC,levels=top_pc_associations$PC)]

pdf(paste(outplot_dir,"qc-top-pc-association.pdf",sep=''),width=7,height=7)
ggplot(top_pc_associations,aes(x=PC,y=value,fill=variable))+
geom_histogram(stat='identity', position=position_dodge(width=0.5))+
xlab('PC')+
ylab('-log10 pval of most associated factor')+
geom_text(aes(label=variable), position='dodge', vjust=-0.25,size=5)+
geom_hline(yintercept=-log10(0.05),linetype='dashed',size=0.5)+
theme_classic()+
  theme(
    legend.position = "bottom",
    axis.ticks.x =element_blank()
  )
dev.off()

## Heatmap 
pdf(paste(outplot_dir,"heatmap-quantnorm-pearson-corr.pdf",sep=''),width=7,height=7)
Heatmap(cor(log2(quantnorm_rpkm_counts$counts+prior_counts),method='pearson'), name = "Pearson corr", col=viridis(10)) 
dev.off()

# pdf(paste(outplot_dir,"heatmap_quantnorm_spearman_corr.pdf",sep=''),width=7,height=7)
# Heatmap(cor(log2(quantnorm_rpkm_counts$counts+prior_counts),method='spearman'), name = "Spearman corr", col=viridis(10)) 
# dev.off()

##------------
## DA testing
##------------
## find best fitting model
model_species_only <- model.matrix(~0+quantnorm_rpkm_counts$samples$group)  
model_species_frip <- model.matrix(~0+quantnorm_rpkm_counts$samples$group+frip_scores$frip)  
model_species_sex <- model.matrix(~0+quantnorm_rpkm_counts$samples$group+quantnorm_rpkm_counts$sex)  
model_species_sex_frip <- model.matrix(~0+quantnorm_rpkm_counts$samples$group+quantnorm_rpkm_counts$sex+frip_scores$frip)  

designlist <- list(
  species_only = model_species_only,
  species_frip = model_species_frip,
  species_sex = model_species_sex,
  species_sex_frip = model_species_sex_frip
#   Full=cbind(Int=1,A=A,B=B,AB=A*B)
)

compare_models <- selectModel(quantnorm_rpkm_counts$counts,designlist)

pdf(paste(outplot_dir,"qc-compare-model-design.pdf",sep=''),width = 7, height = 7)
aic_score =  as.data.table(table(compare_models$pref))
colnames(aic_score) = c('model_design','aic_score')
ggplot(aic_score,aes(x=reorder(model_design,-aic_score),y=aic_score)) +
geom_histogram(stat='identity')+
xlab('model design')+ylab('AIC score')+
theme_bw()
dev.off()

## get best model with lowest aic score
colnames(model_species_sex_frip) <- c("chimp", "human",'sex_M','frip')

y <- estimateDisp(quantnorm_rpkm_counts, model_species_sex_frip) 
fit <- glmQLFit(y, model_species_sex_frip, robust=TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp,levels=model_species_sex_frip))
# results <- glmTreat(fit, contrast=makeContrasts(human-chimp,levels=model_species_sex_frip),lfc = log2(2))

## get df with final results_df
## peaks to keep 
human_peaks_to_keep <- copy(as.data.table(final_hg38_peaks))[peakID %in% rownames(results$table),][,c('width','strand'):=NULL]
chimp_peaks_to_keep <- copy(as.data.table(final_panTro5_peaks))[peakID %in% rownames(results$table),][,c('width','strand'):=NULL]

## these are in hg38 coordinates
final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'peakID'= human_peaks_to_keep$peakID)%>%as.data.table()

pval_threshold = 0.01

final_results <-  final_results[
        ,FDR:=p.adjust(final_results$PValue,method = 'fdr')
        ][
            ,DA:=as.factor(ifelse(FDR<=pval_threshold, 'da','non_da'))
            ][
                    ,da_species := ifelse(DA == 'non_da' ,'non_da',ifelse(DA=='da'& logFC <0,'chimp','human'))
] 

nrow(final_results[DA=='da'])
# [1] 25774
nrow(final_results[DA=='da'])/nrow(final_results)
# [1] 0.1858135
nrow(final_results[DA=='da'][logFC>0])/nrow(final_results[DA=='da'])
# [1] 0.5716614

final_results <- final_results[human_peaks_to_keep,on='peakID',nomatch=0]

## QCs
## check distribution Pvalues
pdf(paste(outplot_dir,"qc-distribution-raw-pvals.pdf",sep=''),width = 7, height = 7)
ggplot(final_results,aes(x=round(PValue,4))) +geom_histogram(binwidth=0.01)+
xlab('raw pvalue')+
theme_classic()+
  theme(
    legend.position = "bottom",
    axis.ticks.x =element_blank()
  )
dev.off()

## check distribution peak lengths between da and nonda peaks
pdf(paste(outplot_dir,"qc-da-peak-sizes.pdf",sep=''),width = 7, height = 7)
df <- copy(final_results)[,width:=end-start]
ggplot(df,aes(x=width,fill=da_species)) +
geom_density(alpha=0.5) +
scale_fill_manual(values=da_species_palette)+
theme_classic()+
  theme(
    legend.position = "bottom",
    axis.ticks.x =element_blank()
  )
dev.off()

##--------------------
## Generate MA plot
##--------------------
ma_plot = function(df){
    plot <- ggplot(df,aes(x = logCPM, y = logFC, col = da_species)) + 
        geom_point(alpha=0.2) + 
        scale_color_manual(values = da_species_palette) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess",se = F) + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)+
        theme_classic()+
        theme(
            legend.position = "bottom",
            axis.ticks.x =element_blank()
            )
        return(plot)
}

pdf(paste(outplot_dir,"ma-plot.pdf",sep=''),width = 7, height = 7)
ma_plot(final_results)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x = logFC,y =-log10(FDR), col = da_species),alpha=0.2)+
    geom_hline(yintercept=-log10(pval_threshold), linetype='dashed', color='black', size=0.5)+
    geom_vline(xintercept=0, linetype='dashed', color='black', size=0.5)+
    scale_colour_manual(values=da_species_palette)+
    theme_classic()+
    xlim(-8,+8)+
    theme(
        legend.position = "bottom",
        axis.ticks.x =element_blank()
        )
    return(plot)
}

pdf(paste(outplot_dir,"volcano-plot.pdf",sep=''),width = 7, height = 7)
volcano_plot(final_results)
dev.off()

## export DA resutls
write.table(final_results,paste(outdir,'files/da_results.txt',sep=''),sep='\t',col.names=T,row.names = F,quote=F)

