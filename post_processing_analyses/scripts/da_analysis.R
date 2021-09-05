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

bamDir  = "output/Post_alignment/Files"
result_dir = './post_processing_analyses/output'
peakDir = paste(result_dir,"/final_peak_set/",sep='')
plot_dir = paste0(result_dir,'/plots/DA/',sep='')
peak_outdir = paste0(result_dir,'/DA/peaks/',sep='')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

##------------------------------------------------------------------
## Read orth common and species-specific peaks (all in hg38 coord) 
##------------------------------------------------------------------
read_peaks = function(file){
    peak = fread(paste(peakDir,file,sep=''),sep='\t',header=T)%>%unique()
    peak = peak[,peakID:=paste(peakID,species,sep='.')]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
    return(peak)
}

# human_and_common_peaks = read_peaks('common_and_human_specific_peaks_hg38.bed')
# chimp_and_common_peaks = read_peaks('common_and_chimp_specific_peaks_pantro5.bed')
# common_peaks_hg38coord = fread(paste(peakDir,'common_peaks_hg38coord.bed',sep=''),sep='\t',header=T)

human_peaks_hg38 = read_peaks('human_peaks_hg38.bed')
chimp_peaks_pantro5 = read_peaks('chimp_peaks_pantro5.bed')

## convert the species-specific peaks to the other genome coordinates and then attach the granges
human_peaks_pantro5 = convert_coord(as.data.table(human_peaks_hg38),'hg38ToPanTro5.over.chain')%>%unique()
human_peaks_pantro5 = human_peaks_pantro5[,i.strand:=NULL]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

chimp_peaks_hg38 = convert_coord(as.data.table(chimp_peaks_pantro5),'panTro5ToHg38.over.chain')%>%unique()
chimp_peaks_hg38 = chimp_peaks_hg38[,i.strand:=NULL]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

all_hg38_peaks = c(human_peaks_hg38,chimp_peaks_hg38)
all_pantro5_peaks = c(chimp_peaks_pantro5,human_peaks_pantro5)

##---------------------------------------
## Get the combined filtered BAM files
##---------------------------------------
## NB: I've removed the reads overlapping blacklisted regions in the snakemake pipeline 
## because if u combine the blacklist regions here
## u risk to remove reads overlapping blacklisted peaks in the other species

standard.chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag = 2000, pe = "both")

get_bamReads = function(organimsDir){
    bamReads = list.files(paste0(organimsDir,bamDir), 
        recursive = T,full.names = T,pattern="*_blacklist_removed.bam$")
    return(bamReads)
}

hg38_bams = get_bamReads('./hg38/')[-7]
pantro5_bams = get_bamReads('./pantro5/')[-7]

##-----------------------------------------------
## Count number of reads overlapping peaks 
##-----------------------------------------------
## Here i am not counting reads overlapping a slinding window of fixed width but rather I am using called peaks
## importantly, I have converted the set of common and species-specific peaks in both pantro5 and hg38 coordinates
## and I am counting the reads in the pantro5 or hg38 bams overlapping the respective set of peaks independently
## this is critical because even if MACS called a peak in a species but not in the other it doenst mean there are no reads overlapping that region

## this function returns a matrix containing the count for each library (column) at each region (row)
## importantly, to account for differences in peak sizes I am using RPKM rather than CPM
get_overlaps = function(bams,peaks,group){
    overlaps =  regionCounts(bams,peaks,param=param)
    
    samples = sub("\\_.*", "", sub('.*\\/', '', bams))
    samples = gsub("*.bam", "", samples)
    
    overlaps_DGE <- asDGEList(overlaps,group=rep(group,each=6))

    peak_width = copy(peaks)%>%as.data.table()%>%dplyr::pull('width')
    overlaps_DGE$counts = rpkm(overlaps_DGE$counts,peak_width) ## it passed manual check

    rownames(overlaps_DGE) = peaks$peakID
    colnames(overlaps_DGE) = samples

    return(overlaps_DGE)
} 

hg38_counts <- get_overlaps(hg38_bams, all_hg38_peaks,1)
pantro5_counts <- get_overlaps(pantro5_bams, all_pantro5_peaks, 2)

## combine the results for counts of common regions where there is 1-to-more match
## caused by liftover splitting a peak into multiple ones
## do it separately for common and species specific counts
## for species-specific simply cbind the matrices
## for the common ones follow comments below
get_counts = function(x,organism){
    count_matrix = copy(x)%>%as.data.frame()%>%mutate('peakID'=rownames(.))%>%as.data.table()
    count_matrix = count_matrix[!peakID %like% organism]
    return(count_matrix)
}

hg38_specific_counts = get_counts(hg38_counts$counts,'common')
pantro5_specific_counts = get_counts(pantro5_counts$counts,'common')

## during liftover a peak can get splitted into multiple ones
## so this function takes the species-specific peaks 
## it splits the species2 (the liftovered one) dt by peakID
## and gets avg of all numeric columns 
## then it reconstruct the dt
## it combines the counts for the species specific peaks of the two species into a single matrix
# combine_specific_counts = function(species1,species2,organism){
#     species1_specific = copy(species1)[peakID %like% organism]
#     species2_specific = copy(species2)[
#         peakID %like% organism
#         ][
#             , c("species_peakID", "species",'numb') := tstrsplit(peakID, ".", fixed=TRUE)
#             ][
#                 ,numb:=NULL
#                 ]%>% split(by='species_peakID')%>%
#                 lapply(function(y)y=y%>%summarize_if(is.numeric, mean, na.rm=TRUE)%>%as.data.table())
                
#     species2_specific = Map(mutate,species2_specific,peakID=names(species2_specific))%>%rbindlist()
#     species2_specific = species2_specific[,peakID:=paste(peakID,organism,sep='.')]

#     species_specific_combined_counts = inner_join(species1_specific,species2_specific,by='peakID')%>%
#     dplyr::select(c(contains('H'),contains('C'),'peakID'))
#     matrix_rownames = species_specific_combined_counts$peakID
#     species_specific_combined_counts = species_specific_combined_counts[,peakID:=NULL]%>%as.matrix()
    
#     rownames(species_specific_combined_counts) = matrix_rownames

#     return(species_specific_combined_counts)

# }
# ## here it is important to invert species1 and species2!!
# human_specific_counts = combine_specific_counts(hg38_specific_counts,pantro5_specific_counts,'human')
# chimp_specific_counts = combine_specific_counts(pantro5_specific_counts,hg38_specific_counts,'chimp')

# ##------------------------------------------------------------------
# ## now merge the counts for peaks in common between species
# ##------------------------------------------------------------------
# ## get common peakIDs for the 2 species
# common_peaks_humans_peakID = copy(common_peaks_hg38coord)[,c('human_peakID','chimp_peakID')]%>%setnames(old='human_peakID',new='peakID')
# common_peaks_chimp_peakID = copy(common_peaks_hg38coord)[,c('human_peakID','chimp_peakID')]%>%setnames(old='chimp_peakID',new='peakID')

# ## use this function to retain from the DGE object only the counts for the common peaks
# get_counts_common_peaks = function(x,common_peaksID){
#     common_counts = copy(x) %>%as.data.frame()%>%mutate('peakID'=rownames(.))%>%as.data.table()
#     common_counts = common_counts[peakID %like% 'common'][
#         ,peakID:=gsub("\\..*","",peakID)
#         ][
#             common_peaksID, on='peakID',nomatch=0]
#     return(common_counts)
# }

# human_common_counts = get_counts_common_peaks(hg38_counts$counts,common_peaks_humans_peakID)
# chimp_common_counts = get_counts_common_peaks(pantro5_counts$counts,common_peaks_chimp_peakID)

# ## join the two count matrices into a single dt
# common_counts = full_join(human_common_counts,chimp_common_counts,by=c('peakID'='human_peakID','chimp_peakID'='peakID'))
# common_counts[is.na(common_counts)]=0
# common_counts = common_counts[,npeaks:=.N,by=.(peakID)]

# ## sum the read counts if a peak in a species corresponds to >1  peaks in the other 
# ## for convenience just take the ID of the first peak as peakID 
# merged_peaks = copy(common_counts)[npeaks==2]
# merged_peaks_ids = copy(merged_peaks)[,c('peakID','chimp_peakID')]%>%split(by='peakID')%>%lapply(function(x)x=x[1,2])%>%rbindlist()
# merged_peaks = merged_peaks[,c(1:6):=NULL]%>%split(by='peakID')
# merged_peaks = lapply(merged_peaks,function(x)x=x[, lapply(.SD, sum), .SDcols=3:8])%>%rbindlist()
# merged_peaks = cbind(merged_peaks,merged_peaks_ids)

# peaks_to_merge = copy(common_counts)[npeaks==2][,c(1:8)]
# peaks_to_merge = peaks_to_merge[merged_peaks,on='chimp_peakID',nomatch=0]

# ## remove from the matrix the duplicated chimp peaks
# common_counts = common_counts[!peaks_to_merge,on='peakID']
# ## now add the merged counts
# common_counts_combined = rbind(common_counts[,npeaks:=NULL],peaks_to_merge)
# common_counts_rownames = paste(paste('H',common_counts_combined$peakID,sep='_'),paste('C',common_counts_combined$chimp_peakID,sep='_'),sep='.')
# common_counts_combined = common_counts_combined[,c('peakID','chimp_peakID'):=NULL]%>%as.matrix()
# rownames(common_counts_combined)=common_counts_rownames

## merge together the common and species-specific counts
# final_count_matrix = DGEList(
#     rbind(common_counts_combined,chimp_specific_counts,human_specific_counts),
#     ,group=rep(1:2,each=6)
# )

final_count_matrix = hg38_specific_counts[pantro5_specific_counts,on='peakID',nomatch=0]
peakIDs=final_count_matrix$peakID
final_count_matrix = final_count_matrix[,peakID:=NULL]%>%as.matrix()
rownames(final_count_matrix)=peakIDs
final_count_matrix = DGEList(final_count_matrix)
##----------
## Some QCs
##----------
## plot density of log-RPKM for each sample and decide cutoff for min number of reads in peaks
pdf(paste(plot_dir,"logRPKM_read_densities.pdf",sep=''),width=7,height=7)
plotDensities(log(final_count_matrix$counts+1), group=rownames(final_count_matrix$samples),col=brewer.pal(12, 'Set3')) 
dev.off()

## Retain only peaks with log-RPKM > -xxx
keep_exprs <- filterByExpr(final_count_matrix$counts,group=final_count_matrix$samples$group)
filtered_counts <- final_count_matrix[keep_exprs,, keep.lib.sizes=FALSE]

## plot density of log-RPKM for each sample
pdf(paste(plot_dir,"logRPKM_filtered_read_densities.pdf",sep=''),width=7,height=7)
plotDensities(log(filtered_counts$counts+1), group=rownames(filtered_counts$samples),col=brewer.pal(12, 'Set3')) 
dev.off()

##------------------
## Normalisation
##------------------
## because there is a bit too much varibility across samples we are using quantile normalisation (voom) rather than TMM
## look at it yourself
tmm_normd_filtered_counts <- calcNormFactors(filtered_counts,method='TMM')

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(plot_dir,"logRPKM_tmm_normalised_filtered.pdf",sep=''),width=7,height=7)
plotDensities(log(tmm_normd_filtered_counts$counts+1), group=rownames(tmm_normd_filtered_counts$samples),col=brewer.pal(12, 'Set3')) 
dev.off()

## quantile norm with voom
filtered_counts_quant_norm = copy(filtered_counts)
filtered_counts_quant_norm$counts =normalizeQuantiles(filtered_counts_quant_norm$counts, ties=TRUE)

## plot the density of log-RPKM for the tmm norm filtered counts for each sample
pdf(paste(plot_dir,"logRPKM_quant_normalised_filtered.pdf",sep=''),width=7,height=7)
plotDensities(log(filtered_counts_quant_norm$counts+1), group=rownames(filtered_counts_quant_norm$samples),col=brewer.pal(12, 'Set3')) 
dev.off()

##-------------------
## Run PCA 
##-------------------
## check with PCA if samples cluster differently
rpkm_count_Tmatrix = t(filtered_counts_quant_norm$counts)

## run PCA
pca <- pca(rpkm_count_Tmatrix,ncomp=5,scale=T) 
group = c(rep('human',6),rep('chimpanzee',6))

## plot PCA results
pdf(paste0(plot_dir,'PCA.pdf'),width = 7, height = 7)
plotIndiv(pca,group=group,legend=T,ind.names = T)
dev.off()

pdf(paste0(plot_dir,'PCA_comp13.pdf'),width = 7, height = 7)
plotIndiv(pca,group=group,comp = c(1,3),legend=T,ind.names = T)
dev.off()

pdf(paste0(plot_dir,'PCA_comp23.pdf'),width = 7, height = 7)
plotIndiv(pca,group=group,comp = c(2,3),legend=T,ind.names = T)
dev.off()

## calculate correaltion across samples for reads in peaks
## separate species-specific from common peaks
## but also make a single comprehensive heatmap
library(ComplexHeatmap)
library(viridis)

pdf(paste(plot_dir,"PearsonCor_read_in_all_peaks.pdf",sep=''),width=7,height=7)
Heatmap(cor(filtered_counts_quant_norm$counts,method='pearson'), name = "Pearson corr", col=viridis(100)) 
dev.off()

##------------------------------------------
## DIFFERENTIAL ACCESSIBILITY ANALYSIS
##------------------------------------------
## setup design matrix and fit model to data
design_matrix <- model.matrix(~0+group, data=filtered_counts_quant_norm$samples)  ## 0+ in the formula is an instruction not to include an intercept column and instead to include a column for each group
colnames(design_matrix) <- c("human", "chimp")

## model the dispersion/variance in accessibility assuming most regions arent DA  
## and stabilize dispersion estimates with empirical bayes
## edgeR determines the common dispersion across features.
## It determines the common variance across features and uses this to model the mean-variance relationship.
## This is then used to calculate a dispersion estimate per feature which is the necessary to test for DE/DA

## The variance for each feature measures the degree of inter-library variation for that feature
## edgeR calculates the common dispersion across feature which, in turn, gives an idea of the overall variability across the genome for the dataset

y <- estimateDisp(filtered_counts_quant_norm, design_matrix) 

## Fits a quasi-likelihood negative binomial generalized log-linear model to data
fit <- glmQLFit(y, design_matrix, robust=TRUE)

## test for DA
results <- glmQLFTest(fit, contrast=makeContrasts(human-chimp, levels=design_matrix))

## correct for multiple hypothesis
fdr_results = getBestTest(rownames(results$table),results$table)
fdr_results = data.table(peakID=rownames(fdr_results),FDR=fdr_results$FDR,direction=fdr_results$direction)

## Concatenate all relevant statistical data and peak info (e.g. regions and species) into a single data.table
final_results = data.frame(results$table)
final_results = dplyr::mutate(final_results,'peakID'=rownames(final_results))%>%as.data.table()
final_results = final_results[
    # genomic_regions,on='peakID',nomatch=0
    # ][
        fdr_results,on='peakID',nomatch=0
        ][
            ,significant:=ifelse(FDR<0.05,'significant','non_significant')
            ][
                ,species:=ifelse(peakID %in% rownames(common_counts_combined),'common',
                            ifelse(grepl("c",peakID),'chimp','human'))
            ][
                , c("human_peakID", "chimp_peakID") := tstrsplit(peakID, ".", fixed=TRUE)
                ][
                    ,human_peakID := ifelse(species=='common',sub(".*?_", "", human_peakID),peakID)
                    ][
                        ,chimp_peakID := ifelse(species=='common',sub(".*?_", "", chimp_peakID),peakID)
]

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

pdf(paste0(plot_dir,'rpkm_quant_norm_normalisation_plot.pdf'),width = 7, height = 7)
ma_plot(final_results)
dev.off()

## Volcano plot for DA regions
## color points with species they belong to
volcano_plot <-function(df){
    plot <- ggplot(df) + 
    geom_point(aes(x=logFC,y=-log10(FDR),col =species))+
    geom_hline(yintercept=1.3, linetype='dashed', color='black', size=1)+facet_wrap(species~.,ncol=3)
    
    return(plot)
}

pdf(paste0(plot_dir,'rpkm_quant_norm_volcano_humanVSchimp.pdf'),width = 10, height = 7)
volcano_plot(final_results)
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
