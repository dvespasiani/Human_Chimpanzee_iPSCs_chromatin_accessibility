## input: the human and chimpanzee peaks called on the merged bams in the respective coordinates 
## output: the set of peaks that have orthologous regions in the other species
## qc: check how many individuals support each consensus peak
## use liftover back and forth

library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(liftOver)
library(rtracklayer)
library(UpSetR)
library(ComplexHeatmap)
library(viridis)
library(edgeR)
library(GenomicRanges)
library(csaw)
library(mixOmics)
library(preprocessCore)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility')

scripts_dir <- './post_processing_analyses/scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

# peakDir = "output/PeakCalling/Files/"
peak_inputdir <- './post_processing_analyses/output/files/orthologous_consensus_filtered_peaks/'
outplot_dir <- create_dir(paste('./post_processing_analyses/',plot_dir,sep=''),'DA')

human_chimp_col <- c('#219ebc','#fb8500')
names(human_chimp_col) = c('chimp','human')

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

##-------------------------------------
## Read filtered orth consensus peaks 
##-------------------------------------
read_peaks <- function(peak_file){
    peaks <- fread(paste(peak_inputdir,peak_file,sep=''),header=T)
    setkeyv(peaks,range_keys)
    return(peaks)
}

human_peaks <- read_peaks('human_filtered_orth_consensus_peaks.txt')
chimp_peaks <- read_peaks('chimp_filtered_orth_consensus_peaks.txt')

## convert the species-specific peaks to the other genome coordinates and then attach the granges
human_ortfilt_chimp_coord <- convert_coord(human_peaks,'hg38ToPanTro5.over.chain')

orth_peaks_overlap = foverlaps(human_ortfilt_chimp_coord,chimp_peaks,type='any')%>%na.omit()
orth_peaks_overlap = orth_peaks_overlap[,width:=i.end-i.start][width>=150 & width<=2000]
duplicated_human_chimp_peaks <- copy(orth_peaks_overlap)[,c('peakID','i.peakID')]%>%unique()

duplicated_human_chimp_peaks =duplicated_human_chimp_peaks[,duphuman:=.N,by=.(i.peakID)][,dupchimp:=.N,by=.(peakID)][duphuman == 1 & dupchimp==1]

orth_peaks_overlap = orth_peaks_overlap[peakID %in% duplicated_human_chimp_peaks$peakID]%>%setnames(old=c('peakID','i.peakID'),new=c('chimp_peakID','human_peakID'))
orth_peaks_overlap = orth_peaks_overlap[,peakID_info:=paste(chimp_peakID,human_peakID,sep='.')]

final_human_peaks = copy(human_peaks)[peakID %in%orth_peaks_overlap$human_peakID]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
final_chimp_peaks = copy(chimp_peaks)[peakID %in%orth_peaks_overlap$chimp_peakID]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

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
raw_counts_human <- regionCounts(hg38_bams, final_human_peaks,param=param)%>%asDGEList()
rownames(raw_counts_human) = orth_peaks_overlap$peakID_info
colnames(raw_counts_human) = samples_names[7:12]

human_peaks_width <- copy(final_human_peaks@ranges@width)

rpkm_counts_human <- copy(raw_counts_human)
rpkm_counts_human$counts <- rpkm(raw_counts_human$counts,gene.length=human_peaks_width)

raw_counts_chimp <- regionCounts(pantro5_bams, final_chimp_peaks,param=param)%>%asDGEList()
rownames(raw_counts_chimp) = orth_peaks_overlap$peakID_info
colnames(raw_counts_chimp) = samples_names[1:6]

chimp_peaks_width <-copy(final_chimp_peaks@ranges@width)
rpkm_counts_chimp <- copy(raw_counts_chimp)
rpkm_counts_chimp$counts <- rpkm(raw_counts_chimp$counts,gene.length=chimp_peaks_width)

## merge matrices 
merged_raw_counts = cbind(raw_counts_chimp,raw_counts_human)
merged_raw_counts$samples <- mutate(merged_raw_counts$samples,'group'=ifelse(rownames(merged_raw_counts$samples) %like% 'C','chimp','human'))

merged_rpkm_counts = cbind(rpkm_counts_chimp,rpkm_counts_human)
merged_rpkm_counts$samples <- mutate(merged_rpkm_counts$samples,'group'=ifelse(rownames(merged_rpkm_counts$samples) %like% 'C','chimp','human'))

##---------------------
## Filtering
##---------------------
prior_counts = 0.5

# cpm_counts <- copy(merged_raw_counts)
# cpm_counts$counts <- cpm(cpm_counts$counts,log=T,prior.counts = prior_counts)

## plot log2 CPM before and after filtering
pdf(paste(outplot_dir,"density_reads_in_peaks.pdf",sep=''),width=7,height=7)
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

keep <- filterByExpr(merged_raw_counts,large.n=3)
# keep <- rowMeans(cpm_counts$counts) > 2
# keep <- keep[apply(merged_raw_counts$counts,1, function(x) all(x!=0))),]
raw_counts_filtered <- copy(merged_raw_counts)[keep,, keep.lib.sizes=FALSE]
rpkm_counts_filtered <- copy(merged_rpkm_counts)[keep,, keep.lib.sizes=FALSE]
dim(merged_raw_counts)
dim(raw_counts_filtered)

# filtered_peak_width <- copy(peak_width)[names(peak_width) %in% rownames(raw_counts_filtered) ]

# rpkm_counts_filtered <- copy(raw_counts_filtered)
# rpkm_counts_filtered$counts <- rpkm(rpkm_counts_filtered$counts,gene.length = filtered_peak_width)

# df[rowSums(df > 10) >= 1, ]

pdf(paste(outplot_dir,"density_reads_in_peaks_filtered.pdf",sep=''),width=7,height=7)
plot_densities(make_dt(cpm(raw_counts_filtered$counts,log=T,prior.counts = prior_counts)))
dev.off()

##------------------
## Normalisation
##------------------
# cpm_counts_filtered = copy(raw_counts_filtered)
# cpm_counts_filtered$counts = cpm(cpm_counts_filtered$counts,log=T,prior.counts= prior_counts)

pdf(paste(outplot_dir,"distribution_unnormalised_filtered_counts.pdf",sep=''),width=7,height=7)
plot_boxes(make_dt(cpm(raw_counts_filtered$counts,log=T,prior.counts = prior_counts)))
dev.off()

## tmm normalisation
tmm_norm_counts <- calcNormFactors(copy(raw_counts_filtered),method='TMM')

pdf(paste(outplot_dir,"distribution_tmm_normalised_filtered_counts.pdf",sep=''),width=7,height=7)
plot_boxes(make_dt(cpm(tmm_norm_counts$counts,log=T,prior.counts = prior_counts)),ylab='log2 CPM')
dev.off()

## loess
loess_norm_counts <- copy(raw_counts_filtered)
loess_norm_counts$counts <- normalizeBetweenArrays(loess_norm_counts$counts, method = "cyclicloess")

pdf(paste(outplot_dir,"distribution_loess_normalised_filtered_counts.pdf",sep=''),width=7,height=7)
plot_boxes(make_dt(cpm(loess_norm_counts$counts,log=T,prior.counts = prior_counts)),ylab='log2 CPM')
dev.off()

## get RPKM and quantile normalise these
quantnorm_rpkm_counts <- copy(rpkm_counts_filtered)
quantnorm_rpkm_counts$counts <- normalizeBetweenArrays(quantnorm_rpkm_counts$counts, method = "quantile")
# quant_norm_counts <- voom(raw_counts_filtered$counts, normalize.method = "quantile")

# quantnorm_raw_counts <- copy(raw_counts_filtered)
# quantnorm_raw_counts$counts <- normalizeBetweenArrays(quantnorm_raw_counts$counts, method = "quantile")
# # quant_norm_counts <- voom(raw_counts_filtered$counts, normalize.method = "quantile")

pdf(paste(outplot_dir,"distribution_quantile_normalised_filtered_counts.pdf",sep=''),width=7,height=7)
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
    scale_color_manual(values = human_chimp_col) +
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

pdf(paste(outplot_dir,"pca_screeplot.pdf",sep=''),width = 7, height = 7)
plot(pca)
dev.off()

pdf(paste(outplot_dir,"pca_quantnorm_counts.pdf",sep=''),width = 7, height = 7)
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

# top_pc_associations <- top_pc_associations[
#     ,maximum_column :=  names(.SD)[max.col(.SD)], .SDcols = 2:5
#     ][
#         , .(Max = do.call(max, .SD)), .SDcols = 2:5, .(PC,maximum_column)
#         ][
#             ,Max:=round(as.numeric(Max),2)
#             ]
pdf(paste(outplot_dir,"pca_top_pc_association.pdf",sep=''),width=7,height=7)
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
pdf(paste(outplot_dir,"heatmap_quantnorm_pearson_corr.pdf",sep=''),width=7,height=7)
Heatmap(cor(log2(quantnorm_rpkm_counts$counts+prior_counts),method='pearson'), name = "Pearson corr", col=viridis(10)) 
dev.off()

pdf(paste(outplot_dir,"heatmap_quantnorm_spearman_corr.pdf",sep=''),width=7,height=7)
Heatmap(cor(log2(quantnorm_rpkm_counts$counts+prior_counts),method='spearman'), name = "Spearman corr", col=viridis(10)) 
dev.off()

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

y <- estimateDisp(quantnorm_rpkm_counts, model_species_sex_frip) 
fit <- glmQLFit(y, model_species_sex_frip, robust=TRUE)
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp,levels=model_species_sex_frip))
# results <- glmTreat(fit, contrast=makeContrasts(human-chimp,levels=model_species_sex_frip),lfc = log2(1.5))

## get df with final results_df
peaks_df <- copy(orth_peaks_overlap)

final_results <- data.frame(results$table)
final_results <- dplyr::mutate(final_results,'peakID'= rownames(final_results))%>%as.data.table()
final_results <- final_results[, c("chimp_peakID", "human_peakID") := tstrsplit(peakID, ".", fixed=TRUE)][,peakID:=NULL]

pval_threshold = 0.01

final_results <-  final_results[
        ,FDR:=p.adjust(final_results$PValue,method = 'fdr')
        ][
            ,DA:=as.factor(ifelse(FDR<=pval_threshold, 'da','non_da'))
            ][
                    ,da_species := ifelse(DA == 'non_da' ,'non_da',ifelse(DA=='da'& logFC <0,'chimp','human'))
] 

nrow(final_results[DA=='da'])
nrow(final_results[DA=='da'])/nrow(final_results)
nrow(final_results[DA=='da'][logFC>0])/nrow(final_results[DA=='da'])

final_results <- inner_join(
    final_results,as.data.table(final_human_peaks),by=c('human_peakID'='peakID')
    )%>%setnames(
        old=c(range_keys),new=c(paste('human_',range_keys,sep=''))
        )%>%inner_join(
            as.data.table(final_chimp_peaks),by=c('chimp_peakID'='peakID')
            )%>%setnames(
                old=c(range_keys),new=c(paste('chimp_',range_keys,sep=''))
                )%>%dplyr::select(
                    -c(contains('width'),contains('strand'))
)%>%as.data.table()


#             [
#                 peaks_df,on='peakID',nomatch=0
#             ][
#                 ,c('width','strand'):=NULL
#                 ][
#                     ,peak_species := ifelse(DA == 'non_da' ,'non_da',ifelse(DA=='da'& logFC <0,'da_chimp','da_human'))
# ] 

final_colors <- c(human_chimp_col,'grey')
names(final_colors) = c(names(human_chimp_col),'non_da')

## QCs
## check distribution Pvalues
pdf(paste(outplot_dir,"da_test_distribution_raw_pvals.pdf",sep=''),width = 7, height = 7)
ggplot(final_results,aes(x=round(PValue,4))) +geom_histogram(binwidth=0.01)+
xlab('raw pvalue')+
theme_classic()+
  theme(
    legend.position = "bottom",
    axis.ticks.x =element_blank()
  )
dev.off()

## check distribution peak lengths between da and nonda peaks
pdf(paste(outplot_dir,"da_test_peak_sizes.pdf",sep=''),width = 7, height = 7)
df <- copy(final_results)[,width:=human_end-human_start]
ggplot(df,aes(x=width,fill=da_species)) +
geom_density(alpha=0.5) +
scale_fill_manual(values=final_colors)+
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
        scale_color_manual(values = final_colors) + 
        geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess",se = F) + # smoothed loess fit; can add span=0.5 to reduce computation load/time
        geom_hline(yintercept = 0) + labs(col = NULL)+
        theme_classic()+
        theme(
            legend.position = "bottom",
            axis.ticks.x =element_blank()
            )
        return(plot)
}

pdf(paste(outplot_dir,"da_test_ma_plot.pdf",sep=''),width = 7, height = 7)
ma_plot(final_results)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x = logFC,y =-log10(FDR), col = da_species),alpha=0.2)+
    # scale_color_manual(values = da_palette) + 
    geom_hline(yintercept=-log10(pval_threshold), linetype='dashed', color='black', size=0.5)+
    geom_vline(xintercept=0, linetype='dashed', color='black', size=0.5)+
    scale_colour_manual(values=final_colors)+
    theme_classic()+
    xlim(-8,+8)+
    theme(
        legend.position = "bottom",
        axis.ticks.x =element_blank()
        )
    return(plot)
}

pdf(paste(outplot_dir,"da_test_volcano_plot.pdf",sep=''),width = 7, height = 7)
volcano_plot(final_results)
dev.off()

## export DA resutls
write.table(final_results,paste('./post_processing_analyses/',da_dir,'new_da_results.txt',sep=''),sep='\t',col.names=T,row.names = F,quote=F)














# all_liftovered_peaks = rbind(hg38_to_pantro5_back,pantro5_to_hg38_back)[,width:=end-start]
# all_liftovered_peaks = bin_distance(all_liftovered_peaks,all_liftovered_peaks$width)



# read_filtered_consensus_peaks = function(organimsDir){
#     peaks = fread(
#         paste0(organimsDir,paste(peakDir,'filtered_consensus_peaks.txt',sep='')),sep='\t',header=T)%>%
#         setorderv(c(range_keys),1)%>%unique()%>%
#         makeGRangesFromDataFrame(keep.extra.columns=T)%>%as.data.table()
#     peaks =peaks[,peakID:=paste('peak_',1:nrow(peaks),sep='')]
#     return(peaks)
# }

# human_consensus_peaks = read_filtered_consensus_peaks('./hg38/')[,species:='human']
# chimp_consensus_peaks = read_filtered_consensus_peaks('./pantro5/')[,species:='chimp']

# ## plot distrbution of consensus peak sizes
# all_consensus_peaks = rbind(human_consensus_peaks,chimp_consensus_peaks)
# all_consensus_peaks = bin_distance(all_consensus_peaks,all_consensus_peaks$width)

# binned_width_levels = c(
#     '0-49',
#     '50-150',
#     '151-300',
#     '301-450',
#     '451-600',
#     '601-750',
#     '751-900',
#     '901-1000',
#     '1001-1150',
#     '1151-1301',
#     '1301-1500',
#     '1501-2000',
#     '2001-3000',
#     '3001-4000',
#     '4001-5000',
#     '5001-6000',
#     '>6000'
# )

# plot_distr_binned_width = function(x){
#     p <- ggplot(x,aes(x=factor(binned_column,levels=binned_width_levels),fill=species))+
#     geom_bar()+
#     facet_wrap(species~.,scale='free')+
#     theme(
#         legend.position='none',
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
#     return(p)
# }

# ## plot distributon binned width peaks
# pdf(paste(plot_dir,'consensus_peak_width_distribution.pdf',sep=''),width=7,height=5)
# plot_distr_binned_width(all_consensus_peaks)
# dev.off()

# ##---------------------------------------------------------
# ## now liftover and then plot again the size distribution
# ##---------------------------------------------------------

# hg38_to_pantro5 <- convert_coord(human_consensus_peaks,'hg38ToPanTro5.over.chain')%>%unique()
# hg38_to_pantro5_back <- convert_coord(hg38_to_pantro5,'panTro5ToHg38.over.chain')%>%unique() 

# pantro5_to_hg38 <- convert_coord(chimp_consensus_peaks,'panTro5ToHg38.over.chain')[seqnames!='chrM']%>%unique()
# pantro5_to_hg38_back <- convert_coord(pantro5_to_hg38,'hg38ToPanTro5.over.chain')[seqnames!='chrM']%>%unique()

# all_liftovered_peaks = rbind(hg38_to_pantro5_back,pantro5_to_hg38_back)[,width:=end-start]
# all_liftovered_peaks = bin_distance(all_liftovered_peaks,all_liftovered_peaks$width)

# ## plot the distribution of the new lifted peaks to see if there are significant differences
# pdf(paste(plot_dir,'consensus_liftovered_peak_width_distribution.pdf',sep=''),width=7,height=5)
# plot_distr_binned_width(all_liftovered_peaks)
# dev.off()


# ##-----------------------
# ## mappability scores
# ##-----------------------
# ## read map scores
# read_map = function(species){
#     file=dir(paste(species,'/output/Genome_mappability/mappability',sep=''),pattern='.bed',full.names=T)
#     map = fread(file,sep='\t',header=T)[
#         ,seqnames:=ifelse(seqnames%like% 'chr',seqnames,paste('chr',seqnames,sep=''))
#         ]
#     return(map)
# }
 
# hg38_map = read_map('./hg38')
# pantro5_map = read_map('./pantro5')

# ## get map score for peak regions 
# human_peakmap = foverlaps(hg38_map,human_filt_orth_peaks,type='any')%>%na.omit()
# chimp_peakmap = foverlaps(pantro5_map,chimp_filt_orth_peaks,type='any')%>%na.omit()

# ## calculate and plot the mean mappability score for each peak
# human_meanmap = copy(human_peakmap)[,peakID:=paste(peakID,species,sep='.')][ , .(mean_map = mean(map_score)), by = peakID][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]
# chimp_meanmap = copy(chimp_peakmap)[,peakID:=paste(peakID,species,sep='.')][ , .(mean_map = mean(map_score)), by = peakID][, c("peakID", "species") := tstrsplit(peakID, ".", fixed=TRUE)]

# mean_map_scores=rbind(human_meanmap,chimp_meanmap)[,mean_map:=round(mean_map,1)]

# pdf(paste0(plot_dir,'peak_mean_mappability_score.pdf',sep=''),width = 10, height = 5)
# ggplot(mean_map_scores, aes(x=factor(mean_map),fill = species))+
#     geom_bar(position='dodge')+xlab(' ')+ylab('peak mean mappability score')+
#     facet_wrap(species~.,ncol=2)+
#     theme(
#     axis.line = element_blank(),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
#     )
# dev.off()

# ## only retain peaks with high mappability scores (i.e. 1)
# peaks_to_keep =copy(mean_map_scores)[mean_map==1]%>%split(by='species')

# human_peaks_to_keep=copy(peaks_to_keep[[1]])
# chimp_peaks_to_keep=copy(peaks_to_keep[[2]])

# human_peak_set = human_filt_orth_peaks[peakID%in%human_peaks_to_keep$peakID]
# chimp_peak_set = chimp_filt_orth_peaks[peakID%in%chimp_peaks_to_keep$peakID]

# ## write peaks
# ## for the moment avoid wrting headers
# write.table(human_peak_set,paste(out_dir,'human_peaks_hg38.bed',sep=''),sep='\t',col.names=F,quote=F,row.names=F)
# write.table(chimp_peak_set,paste(out_dir,'chimp_peaks_pantro5.bed',sep=''),sep='\t',col.names=F,quote=F,row.names=F)

# ## convert chimp orth peaks to hg38 coords
# chimp_peak_set_hg38 = convert_coord(chimp_peak_set,'panTro5ToHg38.over.chain')
# ## tmp file 
# write.table(chimp_peak_set_hg38,paste(out_dir,'chimp_peaks_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)


# ## check peaks that are in common btwn human/chimp 
# ## use bedtools
# get_common_peaks = "bedtools intersect -f 0.5 -r -a post_processing_analyses/output/final_peak_set/human_peaks_hg38.bed -b post_processing_analyses/output/final_peak_set/chimp_peaks_hg38.bed"
# common_peaks <- fread(cmd=get_common_peaks)

# common_peaks_hg38coord = foverlaps(human_peak_set,chimp_peak_set_hg38,type='any')%>%na.omit()%>%
#     setnames(old=c('start','end','i.start','i.end','peakID','i.peakID'),
#             new=c('chimp_start','chimp_end','human_start','human_end','chimp_peakID','human_peakID'))%>%
#             dplyr::select(-c(contains('species'))
# )

# ##-------------------
# ## common peaks QCs
# ##-------------------
# ## check how many peaks are in common between species 
# ## check how far start common peaks are from each other
# common_peaks_qc = copy(common_peaks_hg38coord)[
#     ,start_dist:=abs(chimp_start-human_start)
# ]
# common_peaks_qc = bin_distance(common_peaks_qc,common_peaks_qc$start_dist)

# ## plot distribution start distances between common peaks
# pdf(paste(plot_dir,'distance_btwn_common_peaks.pdf',sep=''),width=7,height=7)
# ggplot(common_peaks_qc,aes(x=factor(binned_column,levels=binned_width_levels)))+
#     geom_bar()+
#     theme(
#         legend.position='none',
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
# dev.off()

# ## keep as common peaks all those peaks within 300bp of distance from each other
# final_set_common_peaks = copy(common_peaks_qc)[start_dist<=300][,c('start_dist','binned_column'):=NULL]

# ## label common and species-specific peaks
# chimp_peaks_pantro5 = copy(chimp_peak_set)[
#     ,species:=ifelse(peakID %in% final_set_common_peaks$chimp_peakID,'common','chimp')
# ]
# human_peaks_hg38 = copy(human_peak_set)[
#     ,species:=ifelse(peakID %in% final_set_common_peaks$human_peakID,'common','human')
# ]

# ## keep common peaks where both peaks have mappability scores >0.5
# common_peaks_hg38coord = common_peaks_hg38coord[
#     ,chimp_map_score:=ifelse(chimp_peakID%in%chimp_peaks_to_keep$peakID,'high','low')
#     ][
#         ,human_map_score:=ifelse(human_peakID%in%human_peaks_to_keep$peakID,'high','low')
#         ][
#            (chimp_map_score=='high' & human_map_score=='high')
#            ][
#                ,c('chimp_map_score','human_map_score'):=NULL
# ]


# ## write files (these will be input of the DA analysis)
# write.table(human_peak_set,paste(out_dir,'common_and_human_specific_peaks_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
# write.table(chimp_peak_set,paste(out_dir,'common_and_chimp_specific_peaks_pantro5.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)
# write.table(common_peaks_hg38coord,paste(out_dir,'common_peaks_hg38coord.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

# ## For all other downstream analyses in which u compare the common vs species-specific peaks save also a file with the common regions (not peaks)
# ## this file is the same of the common_peaks_hg38coord file but with ranges being regions not the peaks. 
# ## this spares me from repeating the code below over and over again across scripts

# ## common regions
# common_regions_hg38 = copy(common_peaks_hg38coord)[
#     , max_chimp_end:= max(chimp_end), by=.(chimp_peakID)
#     ][
#         , min_chimp_start:= min(chimp_start), by=.(chimp_peakID)
#         ][
#             ,region_start:=ifelse(human_start<min_chimp_start,human_start,min_chimp_start)
#             ][
#                 ,region_end:=ifelse(human_end>max_chimp_end,human_end,max_chimp_end)
#                 ][
#                     ,c('seqnames','region_start','region_end','chimp_peakID','human_peakID')
# ]
# colnames(common_regions_hg38)[1:3] = range_keys

# ## plot width distribution of these regions
# common_regions_hg38_qc = copy(common_regions_hg38)[,width:=end-start]
# common_regions_hg38_qc = bin_distance(common_regions_hg38_qc,common_regions_hg38_qc$width)

# ## plot distribution start distances between common peaks
# pdf(paste(plot_dir,'common_regions_size_distribution.pdf',sep=''),width=7,height=7)
# ggplot(common_regions_hg38_qc,aes(x=factor(binned_column,levels=binned_width_levels)))+
#     geom_bar()+
#     theme(
#         legend.position='none',
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
# dev.off()

# ## write file
# write.table(common_regions_hg38,paste(out_dir,'common_regions_hg38.bed',sep=''),sep='\t',col.names=T,quote=F,row.names=F)

