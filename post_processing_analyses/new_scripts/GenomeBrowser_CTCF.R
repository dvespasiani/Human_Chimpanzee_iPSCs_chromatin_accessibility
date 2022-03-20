## use this script to plot CTCF overlap at TAD boundaries using Gviz
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

# ctcf_dir <- './output/homer/homer_output/ctcf/'
outplot_dir <- create_dir(plot_dir,'GenomeBrowser')
peakDir = './output/DA/peaks/'
tads_dir =  './output/TADs/'
files_dir <- './output/files/'
homer_ctcf_dir <-'output/homer/homer_output/ctcf'

## read common TADs
common_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][species=='Common'][,label:=NULL]%>%setorderv(c(range_keys),1)

common_tads_hg38 <- common_tads_hg38[
    ,tadID:=paste('tad_',seq(1:nrow(common_tads_hg38)),sep='')
    ][
        ,species:=NULL
        ][
            ,dist_upstreamTAD:=start-lag(end),by=.(seqnames) ## calculate distance between consecutive TADs
            ][
                ,dist_downstreamTAD:=lead(start)-end,by=.(seqnames)
                ][
                    ,dist_upstreamTAD:=ifelse(dist_upstreamTAD %in% NA,0,dist_upstreamTAD) 
                    ][
                        ,dist_downstreamTAD:=ifelse(dist_downstreamTAD %in% NA,0,dist_downstreamTAD) 
]

# names(all_tads) = species_names[-2]

## first select TADs that are located within 100kb from other up/downstream TADs
close_tads <- copy(common_tads_hg38)[abs(dist_upstreamTAD)<=100000 & abs(dist_downstreamTAD)<=100000]

## then extend these regions to 300kb (because u want to see whether CTCFs density is higher at boundaries or not)
## get TAD boundaries
boundary <-  100000
boundaries <- copy(close_tads)[
        ,upstream_start:= start-boundary
        ][
            ,upstream_end:= start
            ][
                ,downstream_start:= end 
                ][
                    ,downstream_end:= end+ boundary
]

upstream_boundaries <- copy(boundaries)%>%
    dplyr::select(c('seqnames',contains('upstream'),'tadID'))%>%
    setnames(old=c('upstream_start','upstream_end','dist_upstreamTAD'),new=c(range_keys[-1],'distTAD'))%>%
    mutate('boundary'='upstream')


downstream_boundaries = copy(boundaries)%>%dplyr::select(c('seqnames',contains('downstream'),'tadID'))%>%
    setnames(old=c('downstream_start','downstream_end','dist_downstreamTAD'),new=c(range_keys[c(2,3)],'distTAD'))%>%
    mutate('boundary'='downstream')

tad_boundaries <- rbind(upstream_boundaries,downstream_boundaries)
setkeyv(tad_boundaries,range_keys)

##-----------------------
## read CTCF results
##-----------------------
deepbind_output_file <- dir(paste(files_dir,'deepbind_output',sep=''),recursive=F,full.names=T)[[2]]
ctcf_peak_files <- dir(homer_ctcf_dir,recursive=T,full.names=T,pattern='all')[[2]]

## add DA info and genomic coordinates
da_file <- paste(da_dir,'da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'logFC','DA','peakID','peak_species'))

ctcf_peaks <- fread(ctcf_peak_files,header=T,sep='\t')[,'Motif Name':=NULL]%>%setnames(old=c('FASTA ID','Sequence'),new=c('peakID','sequence'))
ctcf_peaks <- ctcf_peaks[da_results,on='peakID',nomatch=0]

ctcf_affinities <-fread(deepbind_output_file,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(1,2,3),new=c('peakID','sequence','prediction'))

ctcf_peaks <- ctcf_peaks[copy(ctcf_affinities),on=c('peakID','sequence'),nomatch=0]
setkeyv(ctcf_peaks,range_keys)

# ##--------------------
# ## read CTCF peaks
# ##--------------------
# deepbind_output_file <- dir(paste(files_dir,'deepbind_predictions',sep=''),recursive=F,full.names=T,pattern='fimo')
# ctcf_peak_file <- dir(paste(files_dir,'fimo_ctcf',sep=''),recursive=F,full.names=T,pattern='fimo')

# ctcf_peaks <- fread(ctcf_peak_file,header=T,sep='\t')
# ctcf_affinities <- fread(deepbind_output_file,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(2,3),new=c('sequence','prediction'))

# ctcf_peaks <- ctcf_peaks[
#     ctcf_affinities,on=c('peakID','sequence'),nomatch=0
#     ][
#             DA=='non_da'
#             ][
#                 ,c('peak_species','sequence'):=NULL
# ]
# setkeyv(ctcf_peaks,range_keys)


## find the tad with a high ctcf overlaps at its boundary
## separately overalp the ctcf file with the up/downstream boundaries
## then get those regions with max number of overlaps (cherry picking for plotting)
## then select those with minimum distTAD
ctcf_boundary_overlap <- foverlaps(ctcf_peaks,tad_boundaries,type='within')%>%
    na.omit()%>%split(by='boundary')%>%
    lapply(function(y){
        y <- y[
            ,numbpeaks:=.N,by=.(tadID)
            ][
                numbpeaks==max(numbpeaks)
                ][
                    ,tad_row:=as.numeric(gsub('.*_','',tadID))
                    ]
    }
)%>%rbindlist()

## now get the region spanning the TADs u want to plot
## beware that this region will likely be too huge, so you'll have to zoom into it later
# tads_of_interest = purrr::map2(
#     all_tads,ctcf_boundary_overlap,
#     function(x,y)x[c(y$tad_row-1,y$tad_row,y$tad_row+1),]%>%unique()
# )

ctcf_peak_of_interests <- copy(ctcf_boundary_overlap)[tadID %in% 'tad_10']
tads_of_interest <- copy(common_tads_hg38)[
        c(ctcf_peak_of_interests$tad_row-1,ctcf_peak_of_interests$tad_row),
]
## get the region to plot and then all the CTCF peaks that overlap in here
region_to_plot <- copy(tads_of_interest)[tadID %in% 'tad_10'][
            ,start:=start-60000
            ][
                ,end:=end-100000
                ][
                    ,c(..range_keys)
]%>%unique()
setkeyv(region_to_plot,range_keys)

## all the ctcf overlapping within the region
ctcf_overlap_region <- foverlaps(ctcf_peaks,region_to_plot,type='within')%>%na.omit()
ctcf_overlap_region <- ctcf_overlap_region[,c('start','end'):=NULL]%>%setnames(old=c('i.start','i.end'),new=range_keys[-1])

##--------------------------------------------------
## get phastcons score only the region around CTCFs
##--------------------------------------------------
## NB: convert the coordinates first from pantro5 to hg38, get the scores for those regions 
## and then convert them back to pantro5
phastCons = phastCons7way.UCSC.hg38

# expanded_region <- copy(region_to_plot)[
#         ,c(..range_keys,'peakID')
#         ][
#             ,start:=min(start),by=.(peakID)
#             ][
#                 ,end:=max(end),by=.(peakID)
# ]%>%unique()
expanded_region <- copy(region_to_plot)[
        ,list(start = seq(start, end)),by=.(seqnames)
        ][
            ,end:=start
]%>%makeGRangesFromDataFrame(keep.extra.columns=T)


conservation_scores <- gscores(phastCons,expanded_region)%>%as.data.table()
# conservation_scores = lapply(conservation_scores,function(x)x=x[,peakID:=paste('peak_',1:nrow(x),sep='')])

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
# chimp_bigwigs = read_files('./pantro5','fc.signal.bigwig')## pantro5 coords

##-------------------------
## generate genomic tracks 
##-------------------------
# ideoTrack <- purrr::map2(
#     region_to_plot,list('panTro5','hg38'),
#     function(x,y)
#     IdeogramTrack(
#     genome = y,
#     chromosome = unique(x$seqnames),
#     size = 2
#     )
# )

gtrack <- GenomeAxisTrack(
    chromosome = region_to_plot$seqnames,
    start = region_to_plot$start,
    end =  region_to_plot$end,
    size = 3
)

# tad_track <- purrr::map2(
#     tads_of_interest,list('panTro5','hg38'),
#     function(x,y)
#     AnnotationTrack(
#     chromosome = as.character(unique(x$seqnames)),
#     start = as.numeric(x$start),
#     end = as.numeric(x$end),
#     genome = y,
#     name = "TAD"
#     )
# )

tad_track <- AnnotationTrack(
    chromosome = as.character(unique(tads_of_interest$seqnames)),
    start = unique(tads_of_interest$start),
    end = unique(tads_of_interest$end),
    genome = genome,
    fill='#ffecd1',
    col = 'black', 
    size=3,
    name = "TAD"
)

## list of CTCF peaks
ctcf_peak_track <- copy(ctcf_overlap_region)[
        ,col:='darkgreen'
        ][
            ,c(..range_keys,'col')
]%>%setorderv(range_keys,1)%>%makeGRangesFromDataFrame(keep.extra.columns=T)

ctcf_peaks_track <- AnnotationTrack(
        range=ctcf_peak_track,
        genome=genome,
        type='l',
        fill=c(ctcf_peak_track@elementMetadata$col),
        col='white',
        # groups = c(ctcf_peaks[[2]]@elementMetadata$species),
        name = "CTCF peaks"
)
## affinities
ctcf_affinity <- copy(ctcf_overlap_region)[,c(..range_keys,'prediction')]%>%makeGRangesFromDataFrame(keep.extra.columns=T)

ctcf_affinity_track <-DataTrack(
        ctcf_affinity,
        genome=genome,
        fill='#ef233c',
        col='#ef233c',
        type='histogram',
        lwd = 0.1,
        min.width = 3,
        name = "CTCF affinity"
)

## phatcons score
conservation_track <- DataTrack(
        range = makeGRangesFromDataFrame(copy(conservation_scores),keep.extra.columns=T),
        name = "phastCons7way",
        type='hist',
        col.histogram = "#8cb369", 
        fill.histogram = "#8cb369",
        size = 3, 
        genome = genome
)

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

bigwigs_track <- get_track(bigwigs,region_to_plot,c('#f79824','#33a1fd'))

human_tracks <- c(gtrack,tad_track,conservation_track,ctcf_affinity_track,bigwigs_track)

## plot tracks
pdf(paste0(outplot_dir,'human_ctcf_genome_plot.pdf',sep=''),width = 10, height = 10)
plotTracks(
    human_tracks,
    from = region_to_plot$start+6000,
    to =  region_to_plot$end-3000
)
dev.off()

