# ## script used to analyses the distribution of DA peaks nearby lncRNA TSSs
# library(dplyr)
# library(data.table)
# library(magrittr)
# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggthemes)
# library(ggplot2)
# library(ggpubr)

# ## read all and DA peaks
# ## and see if the DA peaks have different GC/CpG content 
# options(width=150)
# set.seed(2021)
# setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

# scripts_dir <- './scripts/'
# source(paste(scripts_dir,'utils.R',sep=''))

# outplot_dir <- create_dir(plot_dir,'GC')

# ## get DA peaks
# da_file <- paste(da_dir,'da_results.txt',sep='')
# da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA'))[,type:=DA][,DA:=NULL]

# ## get independent sets of promoters 
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# all_promoters <- transcripts(txdb)%>%as.data.table()
# all_promoters <- promoters[!seqnames%like% "M|Un|_"][,end:=start][,start:=start-450][,width:=end-start][,c(1:3)][,type:='upstream_txdb']

# chrom_state_dir <- '../data/iPSC_chrom_states_hg38'
# ipsc_promoters <- read_chromstate(chrom_state_dir)[chrom_state %like% 'Tss|Flnk'][,c(1:3)][,type:='ipsc_promoters']%>%unique()

# promoters <- rbind(da_results,ipsc_promoters,all_promoters)

# ## get CpG content
# cpg_content <- copy(promoters)%>%makeGRangesFromDataFrame()%>%Repitools::cpgDensityCalc(organism=Hsapiens)
# cpg_content <- cpg_content[,cpg_content:=cpg]

# ##--------------------------------------------
# ## First calculate CpG content of these peaks
# ##--------------------------------------------
# ## expectations: DA peaks should have higher CpG content than non-DA peaks as they are enriched for TssA 
# peaks_cpg <- copy(da_results)%>%split(by='DA')%>%lapply(
#   function(x){
#     cpg <- copy(x)%>%makeGRangesFromDataFrame()%>%Repitools::cpgDensityCalc(organism=Hsapiens)
#     x<-x[,cpg_content:=cpg]
#     return(x)
#   } 
# )

# ## check distribution CpGs and remove outliers (only for plotting purpouses)
# lapply(peaks_cpg,function(x)summary(x$cpg_content))
# peaks_cpg <- lapply(peaks_cpg,function(x)x<-x[!cpg_content>100])
# # common_peaks = common_peaks[!cpg_content>6000] ## single outlier
# # all_peaks = list(chimp_peaks_pantro5,common_peaks,human_peaks_hg38)

# all_peaks_cpg = copy(peaks_cpg)%>%lapply(function(x)x=x[,c('cpg_content','DA')])%>%rbindlist()

# ## test with permutation whether differences between da and non da peaks are larger 
# ## than what expected by chance. calculate z score for the observed value, i.e.
# ## (observed difference - mean(permuted differences))/ std(permuted differences)

# cpg_permutation <- permute_data(all_peaks_cpg$cpg_content,all_peaks_cpg$DA)

# ## permutation plot
# pdf(paste0(outplot_dir,'permuted_cpg_differences.pdf',sep=''),width = 7, height = 7)
# ggplot(cpg_permutation$permuted_scores,aes(x=permuted_scores))+
# geom_histogram(binwidth=0.05)+
# geom_vline(xintercept=cpg_permutation$observed_value,linetype='dashed',color='red')+
# annotate("text", x=cpg_permutation$observed_value-0.5, y=600, label = paste('zscore = ',round(cpg_permutation$stat_results[[1]],2),' ', '-log10 pval = ',round(-log10(cpg_permutation$stat_results[[2]]),2),sep='')) + 
# xlab('mean difference CpG content non-da vs da peaks')+ylab('counts')+theme_bw()
# dev.off()

# ## plot CpG densities
# pdf(paste0(outplot_dir,'peaks_CpG_content.pdf',sep=''),width = 7, height = 7)
# ggplot(all_peaks_cpg,aes(x=DA,y=cpg_content,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# scale_fill_manual(values=da_palette)+
# xlab('')+ylab('Normalised CpG content')+
# geom_hline(yintercept=median(all_peaks_cpg[DA=='non_da']$cpg_content),linetype='dashed')+
# theme_bw()+
# theme(
#   legend.position = 'none',
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()


# ##-------------------------
# ## Now look at GC content
# ##-------------------------
# peaks_cg <- copy(da_results)%>%split(by='DA')%>%lapply(
#   function(x) 
#   x<-x[
#     ,sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
#     ][
#       !sequence%like%'N' ## do this because u need to remove sequences with Ns 
#       ][
#         ,gc_content:=letterFrequency(DNAStringSet(sequence),"GC")/width(DNAStringSet(sequence))
#         ][
#           ,sequence:=NULL
#           ]
# )

# all_peaks_gc <- copy(peaks_cg)%>%rbindlist()%>%dplyr::select(c('gc_content','DA'))

# gc_permutation <- permute_data(all_peaks_gc$gc_content,all_peaks_gc$DA)

# ## permutation plot
# pdf(paste0(outplot_dir,'permuted_gc_differences.pdf',sep=''),width = 7, height = 7)
# ggplot(gc_permutation$permuted_scores,aes(x=permuted_scores))+
# geom_histogram(binwidth=0.0005)+
# geom_vline(xintercept=gc_permutation$observed_value,linetype='dashed',color='red')+
# annotate("text", x=0, y=1500, label = paste('zscore = ',round(gc_permutation$stat_results[[1]],2),' ', '-log10 pval = ',round(-log10(gc_permutation$stat_results[[2]]),2),sep='')) + 
# xlab('mean difference GC content non-da vs da peaks')+ylab('counts')+theme_bw()
# dev.off()

# ## plot GC content
# pdf(paste0(outplot_dir,'peaks_GC_content.pdf',sep=''),width = 7, height = 7)
# ggplot(all_peaks_gc,aes(x=DA,y=gc_content,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# scale_fill_manual(values=da_palette)+
# geom_hline(yintercept=median(all_peaks_gc[DA=='non_da']$gc_content),linetype='dashed')+
# xlab('')+ylab('Normalised GC content')+
# theme_bw()+
# theme(
#   legend.position = 'none',
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()


# # # random_cres <- sample_n(cres,10000)%>%as.data.table()
# # random_cres <- copy(cres)
# # random_cres <- random_cres[,peakID:=paste('bgpeak_',1:nrow(random_cres),sep='')][,c(..range_keys,'peakID')]


