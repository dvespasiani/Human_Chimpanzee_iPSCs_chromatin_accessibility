## use this script to plot CTCF overlap at TAD boundaries using Gviz
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(liftOver)
library(rtracklayer)
library(Gviz)
library(biomaRt)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)
 

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

chrom_state_dir = './data/iPSC_chrom_states_hg38'
plot_dir = './post_processing_analyses/output/plots/GenomeBrowser/'
tmp_files_dir = './post_processing_analyses/output/temp_files/'
tads_dir =  './post_processing_analyses/output/TADs/'
peakDir = './post_processing_analyses/output/DA/peaks/'

##------------
## read TADs
##------------
human_tads = read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]
chimp_tads = read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

all_tads = list(copy(chimp_tads),copy(human_tads))%>%lapply(function(x)setorderv(x,c(range_keys[-3]),1))
all_tads = lapply(all_tads,function(x)x=x[
    ,tadID:=paste('tad_',seq(1:nrow(x)),sep='')
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
    
) 
names(all_tads) = species_names[-2]

## first select TADs that are located within 100kb from other up/downstream TADs
close_tads = copy(all_tads)%>%lapply(function(x)x=x[abs(dist_upstreamTAD)<=100000 & abs(dist_downstreamTAD)<=100000])

## then extend these regions to 300kb (because u want to see whether CTCFs density is higher at boundaries or not)
## get TAD boundaries
boundary =  100000
boundaries = copy(all_tads)%>%
lapply(
    function(x)x=x[
        ,upstream_start:= start-boundary
        ][
            ,upstream_end:= start
            ][
                ,downstream_start:= end 
                ][
                    ,downstream_end:= end+ boundary
                    ]
)

upstream_boundaries = copy(boundaries)%>%
lapply(
    function(x)x=x%>%dplyr::select(c('seqnames',contains('upstream'),'tadID'))%>%
    setnames(old=c('upstream_start','upstream_end','dist_upstreamTAD'),new=c(range_keys[-1],'distTAD'))%>%
    mutate('boundary'='upstream')
)
# lapply(upstream_boundaries,function(x)setkeyv(x,range_keys))

downstream_boundaries = copy(boundaries)%>%
lapply(
    function(x)x=x%>%dplyr::select(c('seqnames',contains('downstream'),'tadID'))%>%
    setnames(old=c('downstream_start','downstream_end','dist_downstreamTAD'),new=c(range_keys[c(2,3)],'distTAD'))%>%
    mutate('boundary'='downstream')
)
# lapply(downstream_boundaries,function(x)setkeyv(x,range_keys))

tad_boundaries = purrr::map2(upstream_boundaries,downstream_boundaries,rbind)
lapply(tad_boundaries,function(x)setkeyv(x,range_keys))

##--------------------
## read CTCF peaks
##--------------------
ctcf_peaks =list.files(tmp_files_dir,recursive=F,full.names=T,pattern='ctcf')%>%
lapply(function(x)fread(x,sep='\t',header=T)[
    ,c('seqnames','ctcf_peak_start','ctcf_peak_end','MotifScore','prediction','peakID','species')
    ]%>%setnames(old=c(1:3),new=range_keys)
)
names(ctcf_peaks)=paste(species_names,'ctcf',sep='_')

## make 2 files by combining common and species specific
## in their respective coordinates
chimp_peaks = rbind(ctcf_peaks[[1]],copy(ctcf_peaks[[2]]%>%convert_coord('hg38ToPanTro5.over.chain')))
human_peaks = rbind(ctcf_peaks[[3]],ctcf_peaks[[2]])

all_ctcf_peaks = list(chimp_peaks,human_peaks)
names(all_ctcf_peaks) = species_names[-2]

## to find the tad with highest ctcf overlaps at its boundary
## separately overalp the ctcf file with the up/downstream boundaries
## then get those regions with max number of overlaps (cherry picking for plotting)
## then select those with minimum distTAD

## PS this is a bit convoluted, i'll find a neater way to code it
ctcf_boundary_overlap = purrr::map2(all_ctcf_peaks,tad_boundaries,function(x,y)foverlaps(x,y,type='within')%>%na.omit())
ctcf_boundary_overlap = 
lapply(
    ctcf_boundary_overlap,function(x)x%>%split(by='boundary')%>%
    lapply(
        function(y)y=y[
            ,numbpeaks:=.N,by=.(tadID)
        ]%>%split(by='tadID')%>%
        lapply(
            function(z)z=z[
                ,regulation:=length(levels(as.factor(y$species)))
                ][
                    regulation>1 
                    ][
                        ,regulation:=NULL
                        ])%>%rbindlist()
                        )%>%rbindlist()
)

ctcf_boundary_overlap= 
lapply(ctcf_boundary_overlap,function(x)x=x[
    numbpeaks>5
    ][
        distTAD==min(abs(distTAD))
        ]%>%split(by='tadID')%>%
        lapply(function(y)y=y[
            ,numbpeaks_species:=.N,by=.(species)
            ][
                ,numb_common_peaks:=numbpeaks-min(numbpeaks_species)
            ]
        )%>%rbindlist()
)
## these are all the CTCF peaks that overlap at around TAD boundaries
ctcf_boundary_overlap =lapply(ctcf_boundary_overlap,function(x)x=x[numb_common_peaks==max(numb_common_peaks)][,tad_row:=as.numeric(gsub('.*_','',tadID))]%>%dplyr::select(-c(contains('numb'))))

## now get the region spanning the TADs u want to plot
## beware that this region will likely be too huge, so you'll have to zoom into it later
tads_of_interest = purrr::map2(
    all_tads,ctcf_boundary_overlap,
    function(x,y)x[c(y$tad_row-1,y$tad_row,y$tad_row+1),]%>%unique()
)

## get the region to plot and then all the CTCF peaks that overlap in here
region_to_plot = copy(tads_of_interest)%>%
lapply(function(x)x=x[
            ,start:=min(start)-50000
            ][
                ,end:=max(end)+50000
                ][
                    ,c(..range_keys)
                    ]%>%unique()
)
lapply(region_to_plot,function(x)setkeyv(x,range_keys))

## all the ctcf overlapping within the region
ctcf_overlap_region = purrr::map2(all_ctcf_peaks,region_to_plot,function(x,y)foverlaps(x,y,type='within')%>%na.omit())
ctcf_overlap_region =lapply(ctcf_overlap_region,function(x)x=x[,c('start','end'):=NULL]%>%setnames(old=c(1:3),new=range_keys))

##--------------------------------------------------
## get phastcons score only the region around CTCFs
##--------------------------------------------------
## NB: convert the coordinates first from pantro5 to hg38, get the scores for those regions 
## and then convert them back to pantro5
phastCons = phastCons7way.UCSC.hg38

region_cons_score = copy(ctcf_overlap_region)
chimp_region_hg38 = copy(region_cons_score[[1]])%>%convert_coord('panTro5ToHg38.over.chain')

region_hg38 = list(chimp_region_hg38,region_cons_score[[2]])%>%lapply(
    function(x)x=x[
        ,c(..range_keys,'peakID')
        ][
            ,start:=min(start),by=.(peakID)
            ][
                ,end:=max(end),by=.(peakID)]%>%unique()
)
expanded_region_hg38 =copy(region_hg38)%>%
lapply(function(x)x=x[,c(..range_keys,'peakID')]%>%split(by='peakID')%>%
lapply(function(y)y=y[
        ,list(start = seq(start, end)),by=peakID
        ][
            ,end:=start
        ])%>%rbindlist()
)

conservation_scores=purrr::map2(
    expanded_region_hg38,ctcf_overlap_region,function(x,y)
    x=x[,seqnames:=unique(y$seqnames)]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
)

conservation_scores = lapply(conservation_scores,function(x)gscores(phastCons,x)%>%as.data.table())
conservation_scores = lapply(conservation_scores,function(x)x=x[,peakID:=paste('peak_',1:nrow(x),sep='')])

chimp_conscore_pantro5 = copy(conservation_scores[[1]])%>%convert_coord('hg38ToPanTro5.over.chain')
chimp_conscore_pantro5 = chimp_conscore_pantro5[seqnames %in% region_cons_score[[1]]$seqnames]
region_conservation_scores = list(chimp_conscore_pantro5,conservation_scores[[2]])%>%
lapply(function(x)x=x[
    ,c(..range_keys,'default')
    ]%>%makeGRangesFromDataFrame(keep.extra.columns=T)
)
##---------------
## read bigWigs
##---------------
read_files = function(species,patterns){
    files =  list.files(
        paste(species,'/output/PeakCalling/Files',sep=''),
        recursive=T,full.names=T,pattern= patterns)
    return(files)
}

human_bigwigs = read_files('./hg38','fc.signal.bigwig')
chimp_bigwigs = read_files('./pantro5','fc.signal.bigwig')## pantro5 coords
bigwigs = c(chimp_bigwigs,human_bigwigs)

##-------------------------
## generate genomic tracks 
##-------------------------
ideoTrack <- purrr::map2(
    region_to_plot,list('panTro5','hg38'),
    function(x,y)
    IdeogramTrack(
    genome = y,
    chromosome = unique(x$seqnames),
    size = 2
    )
)

gtrack = lapply(
    region_to_plot,function(x)
    GenomeAxisTrack(
    chromosome = unique(x$seqnames),
    start = min(x$start),
    end =  max(x$end),
    size = 2
    )
)

tad_track = purrr::map2(
    tads_of_interest,list('panTro5','hg38'),
    function(x,y)
    AnnotationTrack(
    chromosome = as.character(unique(x$seqnames)),
    start = as.numeric(x$start),
    end = as.numeric(x$end),
    genome = y,
    name = "TAD"
    )
)

## list of CTCF peaks
ctcf_peaks = copy(ctcf_overlap_region)%>%
lapply(
    function(x)x=x[
        ,col:=ifelse(species=='common','darkred','darkgreen')
        ][
        ,c(..range_keys,'species','col')
        ]%>%setorderv(range_keys,1)%>%makeGRangesFromDataFrame(keep.extra.columns=T)
)

ctcf_peaks_track = purrr::map2(
    ctcf_peaks,list('panTro5','hg38'),
    function(x,y)
    AnnotationTrack(
        range=x,
        genome=y,
        type='l',
        fill=c(x@elementMetadata$col),
        col='white',
        # groups = c(ctcf_peaks[[2]]@elementMetadata$species),
        name = "CTCF peaks"
    )
)
## affinities
ctcf_affinity = copy(ctcf_overlap_region)%>%lapply(function(x) makeGRangesFromDataFrame(x[,c(..range_keys,'prediction')],keep.extra.columns=T))

ctcf_affinity_track = purrr::map2(
    ctcf_affinity,list('panTro5','hg38'),
    function(x,y)
    DataTrack(
        x,
        genome=y,
        type='histogram',
        name = "CTCF affinity")
)

## motif scores
ctcf_score = copy(ctcf_overlap_region)%>%lapply(function(x) makeGRangesFromDataFrame(x[,c(..range_keys,'MotifScore')],keep.extra.columns=T))

ctcf_motifscore_track = purrr::map2(
    ctcf_score,list('panTro5','hg38'),
    function(x,y)
    DataTrack(
        x,
        genome=y,
        type='histogram',
        name = "CTCF score")
)

## phatcons score
conservation_track = lapply(
    region_conservation_scores,function(x)
    DataTrack(
        range = x,
        name = "phastCons7way",
        type='hist',
        col.histogram = "darkblue", 
        fill.histogram = "darkblue",
        size = 3
    )
)

## bigwigs
get_track = function(file,region,colors){
    track_colors = as.list(c(rep(colors[[1]],6),rep(colors[[2]],6)))
    
    files = copy(file)
    sample_names = lapply(files,function(z)gsub('\\_.*','',basename(z)))
    sample_names = lapply(sample_names,function(x)gsub('\\..*','',x))

   track = purrr::pmap(
        list(files,sample_names,track_colors),
        function(x,y,z)
        DataTrack(
            range = x, 
            genome = "hg38",
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
bigwigs_track = get_track(bigwigs,region_to_plot[[2]],c('deepskyblue1','dodgerblue3'))

chimp_tracks = c(ideoTrack[[1]],gtrack[[1]],tad_track[[1]],conservation_track[[1]],ctcf_peaks_track[[1]],ctcf_affinity_track[[1]],ctcf_motifscore_track[[1]])
human_tracks = c(ideoTrack[[2]],gtrack[[2]],tad_track[[2]],conservation_track[[2]],ctcf_peaks_track[[2]],ctcf_affinity_track[[2]],ctcf_motifscore_track[[2]])


## plot tracks
pdf(paste0(plot_dir,'chimp_ctcf_genome_plot.pdf',sep=''),width = 7, height = 7)
plotTracks(
    chimp_tracks,
    from = min(region_to_plot[[1]]$start),
    to =  max(region_to_plot[[1]]$end)

)
dev.off()

pdf(paste0(plot_dir,'human_ctcf_genome_plot.pdf',sep=''),width = 7, height = 7)
plotTracks(
    human_tracks,
    from = min(region_to_plot[[2]]$start),
    to =  max(region_to_plot[[2]]$end)
)
dev.off()

