## use this script to plot regions with Gviz
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(Gviz)
library(biomaRt)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)
 

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

chrom_state_dir = '../data/iPSC_chrom_states_hg38'
plot_dir = './output/plots/GenomeBrowser/'
tmp_files_dir = './output/temp_files/'
tads_dir =  './output/TADs/'
peakDir = './output/DA/peaks/'
target_genes_dir <- './output/files/GO_enrich'

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# ##-------------------
# ## Orthologous TADs
# ##-------------------
# orth_tads <- read_tads('human_chimp_orth_tads.txt')[,c(1:3,7)]
# colnames(orth_tads)[1:3]= range_keys
# setkeyv(orth_tads,range_keys)

##---------------
## bigWigs
##---------------
read_files <- function(dir, patterns){
    files <-  list.files(
        paste('../',genome,dir,sep=''),
        recursive=T,full.names=T,pattern= patterns)
    return(files)
}

bigwigs <- read_files('/output/PeakCalling/Files','fc.signal.bigwig')

##----------------------------
## DA peaks w DE genes
##----------------------------
da_peaks_w_de_genes <- fread(paste(tmp_files_dir,'da_peaks_w_de_genes_gviz.txt',sep=''),sep='\t',header=T)
da_peaks_w_de_genes <- da_peaks_w_de_genes[
    ,concordance:=ifelse(logFC<0 & de_logFC<0,'y',ifelse(logFC>0 & de_logFC>0,'y','n'))
    ][
        concordance=='y'
        ][
            ,species:=ifelse(logFC<0,'chimp','human')
]
# da_peaks_w_de_genes_same_tad <- foverlaps(copy(da_peaks_w_de_genes), copy(orth_tads),type='within')%>%na.omit()
# da_peaks_w_de_genes_same_tad <- da_peaks_w_de_genes_same_tad[
#   ,same_tad:=ifelse(distTSS<0 & start- (i.start+distTSS)>0,'y',
#   ifelse(distTSS>0 & end - (i.end+distTSS)>0,'y','n')
#   )
#   ][
#     same_tad=='y'
#     ][
#         ,same_species:=ifelse(species=='Chimp' & peak_species=='da_chimp','y',ifelse(species=='Human' & peak_species=='da_human','y','n'))
#         ][
#             same_species=='y'
#             ][
#                 ,c('same_tad','go_signif','pleiotropy','same_species'):=NULL
# ]%>%unique()%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
# colnames(da_peaks_w_de_genes_same_tad)[2:3] = c('tad_start','tad_end')


# ## get reads and filter for those with highest difference between chimp and humans
# library(csaw)
# library(edgeR)
# peak_file <- paste('output/final_peak_set/',genome,'_all_orthologous_peaks.txt',sep='')
# peaks <- fread(peak_file,sep='\t',header=T)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

# ## list bam
# standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
# param <- readParam(pe = "both",restrict=standard_chr,max.frag=1000)

# get_bams <- function(species){
#     bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
#     return(bams)
# }

# bams <- get_bams(genome)

# reads_in_peaks <- regionCounts(bams, peaks,param=param)

# ## remove peaks with low counts
# abundances <- aveLogCPM(asDGEList(reads_in_peaks))
# keep <- abundances > 0

# reads_in_peaks_filtered <- normFactors(reads_in_peaks[keep,])

# reads_in_peaks_dge <- asDGEList(reads_in_peaks_filtered)
# rownames(reads_in_peaks_dge) <- rowData(reads_in_peaks_filtered)$peakID
# colnames(reads_in_peaks_dge) <- samples_names
# reads_in_peaks_dge$samples <- mutate(reads_in_peaks_dge$samples,'group'=ifelse(rownames(reads_in_peaks_dge$samples) %like% 'C','chimp','human'))

# peak_width <- copy(peaks)%>%as.data.table()
# peak_width <- peak_width[peakID %in% rowData(reads_in_peaks_filtered)$peakID]%>%dplyr::pull('width')
# reads_in_peaks_dge$counts = rpkm(reads_in_peaks_dge$counts,peak_width) 

# ## check which peak has highest rpkm difference between species 
# test <- copy(reads_in_peaks_dge$counts)
# test <- subset(test, rownames(test) %in% da_peaks_w_de_genes_same_tad$peakID)
# cmean <- rowSums(test[,c(1:6)])
# hmean <- rowSums(test[,c(7:12)])

# ## when -ve = higher in human, vice versa when +ve
# diff_means <- data.table(
#     peakID=names(cmean),
#     c_h_mean_diff=cmean-hmean
# )
# top_differences <- copy(diff_means)[abs(c_h_mean_diff)]%>%unique()%>%setorderv('c_h_mean_diff',1)

# top_differences[c_h_mean_diff>0]

# top_different_peaks <- list('peak_218171','peak_196193')

# # test[rownames(test) %in% top_differences[c_h_mean_diff>0]$peakID]
# ## peak_101698,peak_189586,peak_99152,peak_218171 more accessible in chimp. target gene is expressed in developing (but not adult) brain
# ## peak_196193 more accessible in human. target gene involved in development of human germ cells

# peak_gene_of_interest <-copy(da_peaks_w_de_genes_same_tad)%>%split(by='species')
# peak_gene_of_interest <- purrr::map2(peak_gene_of_interest,top_different_peaks,function(x,y){
#     x<-x[peakID %in% y][,tad_start:=min(tad_start)][,tad_end:=max(tad_end)]%>%unique()
#     return(x)
# })

top_logFC<- copy(da_peaks_w_de_genes)[
    ,c('logFC','SYMBOL','peakID','species')
][,rank:=round(logFC/max(abs(logFC)),1)]%>%split(by='species')



peaks = c('peak_226301','peak_118593') # h,c
peak_of_interest <-copy(da_peaks_w_de_genes)[peakID %in%peaks]%>%setorderv('species',1)%>%split(by='species')%>%lapply(function(x){x<-x[,c(..range_keys)]})
gene_of_interest <- copy(da_peaks_w_de_genes)[peakID %in%peaks]%>%setorderv('species',1)%>%split(by='species')%>%lapply(function(x)x$SYMBOL)

## get region of interest
region_of_interest <- copy(da_peaks_w_de_genes)%>%setorderv('species',1)%>%split(by='species')
region_of_interest <- purrr::map2(region_of_interest,peak_of_interest,function(x,y){
    x<-x[y,on=c(range_keys),nomatch=0][
            ,start:=start-10000
            ][
                ,end:=end+10000
                ][
                    ,c(..range_keys,'species')
                    ]%>%unique()
    setkeyv(x,range_keys)
    }
)
# ## chromstate for region to plot
# ipsc_chromstate <- read_chromstate(chrom_state_dir) ## these contain info for sex chr and are in hg38 coord
# setkeyv(ipsc_chromstate,range_keys)

# chromstate_of_interest <- copy(region_of_interest)%>%lapply(
#     function(x){
#         setkeyv(x,range_keys)
#         x<-foverlaps(x,ipsc_chromstate,type='any')%>%na.omit()
#         x <- x[,c(..range_keys,'chrom_state','species')]%>%unique()
#         x <- x[, .SD[which.min(start)], by=.(end)][, .SD[which.max(end)], by=.(start)]
#         x <- x[
#             ,state_numb:=as.numeric(gsub('\\_.*','',chrom_state))
#             ][
#                 ,c(..range_keys,'state_numb')
#                 ]%>%unique()%>%makeGRangesFromDataFrame(keep.extra.columns=T)
#         return(x)
#     }
# )

# ## phastcons score
# phastCons = phastCons7way.UCSC.hg38

# expanded_region <- copy(region_of_interest)%>%lapply(function(x)
#     x<-x[
#         ,list(start = seq(start, end)),by=.(seqnames)
#         ][
#             ,end:=start
#             ]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
# )

# conservation_scores <- lapply(expanded_region,function(x){
#     df <-copy(x)
#     df <-gscores(phastCons,df)%>%as.data.table()
#     df$default[is.na(df$default)] <- 0
#     df <- df[,c('width','strand'):=NULL]
#     return(df)
#     }
# )
##---------------
## read bigWigs
##---------------
read_files <- function(patterns){
    files <-  list.files(
        paste('../',genome,'/output/PeakCalling/Files',sep=''),
        recursive=T,full.names=T,pattern= patterns)
    return(files)
}
bigwigs <- read_files('fc.signal.bigwig')
##----------------------
## Genome Browser Tracks
##----------------------
biomart_track <- copy(region_of_interest)%>%lapply(function(x){
    track<-BiomartGeneRegionTrack(
    chromosome = x$seqnames, 
    start = x$start,
    end = x$end,
    biomart = ensembl,
    genome='hg38',
    name='Target gene',
    size=5,
    min.height = 1)

    return(track)
    }
)

biomart_track <-purrr::map2(biomart_track,gene_of_interest,function(x,y){
    x<-x[x@range$symbol %in% y]
    return(x)
    }
)

## plot region that spans the longest gene transcript
region_to_plot <- copy(biomart_track)%>%lapply(function(x){
    expand_range = 5000
    x <- as.data.table(x@range)
    x <- x[,start:=min(start)-expand_range][,end:=max(end)+expand_range][,c(..range_keys)]%>%unique()
})

gtrack <- lapply(region_to_plot,function(x){
    x<-copy(x)
    track<-GenomeAxisTrack(
    chromosome = x$seqnames,
    start = x$start,
    end =  x$end,
    size = 3)
    return(track)
    }
)

# tad_track <- lapply(tad_of_interest,function(x){
#     x<- AnnotationTrack(
#     chromosome = x$seqnames,
#     start = x$start,
#     end = x$end,
#     genome = genome,
#     fill='#ffecd1',
#     col = 'black', 
#     size=3,
#     name = "TAD")
#     }
# )

# conservation_track <- lapply(conservation_scores,function(x){
#     gr <- makeGRangesFromDataFrame(copy(x),keep.extra.columns=T)
#     track <- DataTrack(
#         range = gr,
#         name = "phastCons7way",
#         type='hist',
#         col.histogram = "#8cb369", 
#         fill.histogram = "#8cb369",
#         size = 3, 
#         genome = genome)
#         return(track)
#         }
# )

# ## chromstate track
# ipsc_chromstate_track <- lapply(chromstate_of_interest,function(x){
#     track <- DataTrack(
#             range = x,
#             genome = genome,
#             chromosome = x$seqnames,start=x$start, end=x$end,
#             type = "heatmap",
#             name = "iPSCs chrom state",
#             gradient = chrom_state_colors,
#             size = 3,
#             min.height = 1)
#     return(track)
#     }
# )
## bigwigs
get_track <- function(file,region,colors){
    track_colors = as.list(c(rep(colors[[1]],6),rep(colors[[2]],6)))
    
    files = copy(file)
    sample_names = lapply(files,function(z)gsub('\\_.*','',basename(z)))
    sample_names = lapply(sample_names,function(x)gsub('\\..*','',x))

   track = purrr::pmap(
        list(files,sample_names,track_colors),
        function(x,y,z)
        DataTrack(
            range = x, 
            genome = genome,
            type = "h", 
            chromosome = unique(as.character(region$seqnames)),
            start = unique(region$start),
            end = unique(region$end),
            name = y,
            col=z,
            col.histogram=z
            )
        )
    return(track)
}

da_chimp_bigwigs_track <- get_track(bigwigs,region_to_plot[[1]],c('#f79824','#33a1fd'))
da_human_bigwigs_track <- get_track(bigwigs,region_to_plot[[2]],c('#f79824','#33a1fd'))

ht <- purrr::map2(peak_of_interest,list(da_chimp_bigwigs_track,da_human_bigwigs_track),function(x,y){
    gr <- copy(x)
    track <- HighlightTrack(
        trackList = y,
        genome=genome,
        start = gr$start,end=gr$end,
        chromosome = gr$seqnames)
        return(track)
        }
)

chimp_tracks <- c(gtrack[[1]],biomart_track[[1]],ht[[1]])
human_tracks <- c(gtrack[[2]],biomart_track[[2]],ht[[2]])

## plot tracks
pdf(paste0(plot_dir,'chimp_da_genome_plot.pdf',sep=''),width = 20, height = 10)
plotTracks(
    chimp_tracks,
    from = region_to_plot[[1]]$start,
    to =  region_to_plot[[1]]$end,
    collapseTranscripts='longest', 
    transcriptAnnotation = 'symbol'
)
dev.off()

pdf(paste0(plot_dir,'human_da_genome_plot.pdf',sep=''),width = 20, height = 10)
plotTracks(
    human_tracks,
    from = region_to_plot[[2]]$start,
    to =  region_to_plot[[2]]$end,
    collapseTranscripts='longest', 
    transcriptAnnotation = 'symbol'
)
dev.off()



