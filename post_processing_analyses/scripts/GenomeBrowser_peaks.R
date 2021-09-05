## use this script to plot regions with Gviz
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
chimp_tads=setorderv(chimp_tads,range_keys,1)
chimp_tads = chimp_tads[,peakID:=paste('peak_',seq(1:nrow(chimp_tads)),sep='')] ## this is temporary, I should change the function

## convert chimp TADs into hg38
chimp_tads_hg38 = convert_coord(chimp_tads,'panTro5ToHg38.over.chain')
chimp_tads_hg38 = chimp_tads_hg38[,peakID:=NULL]%>%unique()

all_tads = list(chimp_tads_hg38,human_tads)
all_tads = lapply(all_tads,function(x) x=x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]) ## use this column only to split dataframes
lapply(all_tads,function(x)setkeyv(x,range_keys))
names(all_tads) = species_names[-2]

##-------------------------
## read iPSCs chrom states
##-------------------------
ipsc_chromstate = read_chromstate(chrom_state_dir)

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
chimp_bigwigs = read_files('./pantro5','fc.signal.hg38.bigwig')

##---------------
## read peaks
##---------------
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))

common_da_results = read_da_peaks('common_da_results.txt',c('significant','non_significant'))
common_da_results = common_da_results[
  ,peakID:=paste(paste('C_',chimp_peakID,sep=''),paste('H_',human_peakID,sep=''),sep='.')
]
human_common_peaks = rbind(
    human_peaks[,c(..range_keys,'peakID')],
    common_da_results%>%dplyr::select(c(contains(paste('human_',range_keys,sep='')),'peakID'))%>%
    setnames(old=c(1:3),new=range_keys)
)
chimp_common_peaks = rbind(
    chimp_peaks[,c(..range_keys,'peakID')],
    common_da_results%>%dplyr::select(c(contains(paste('chimp_',range_keys,sep='')),'peakID'))%>%
    setnames(old=c(1:3),new=range_keys)
)%>%convert_coord('panTro5ToHg38.over.chain')

da_peaks = list(chimp_common_peaks,human_common_peaks)
names(da_peaks) = species_names[-2]

##--------------------------------
## read peaks w genes in same TAD
##--------------------------------
peaks_w_genes_same_tad = fread(paste(tmp_files_dir,'list_peaks_w_genes_same_tad.txt',sep=''),header=T,sep='\t')%>%split(by='file')
peaks_w_genes_same_tad = lapply(
    peaks_w_genes_same_tad,function(x)x=x[
        ,c(..range_keys,'peakID','peak_species','tadID','gene','distTSS','EnsemblID','numbpeaks_per_gene_same_tad')
    ]%>%split(by='gene')%>%
    lapply(function(y)y=y[
    ,regulation:=ifelse(
        length(unique(y$peak_species)) == 1, unique(y$peak_species),'common_and_specific')
        ]
      )%>%rbindlist()%>%dplyr::select(-c(all_of(range_keys))) ## these will come back with the next join (see below)
)

## correct here a potential mistake and change it back in the da_peaks_in_tads.R script
correct = copy(peaks_w_genes_same_tad)%>%lapply(function(x)x=x[,c(1,2,5,7,8):=NULL]%>%unique())
correct = lapply(correct,function(x)x=x[,numbtads_per_gene:=.N,by=.(gene)][numbtads_per_gene==1])

peaks_w_genes_same_tad  = purrr::map2(peaks_w_genes_same_tad,correct,function(x,y)x[tadID %in% y$tadID])

## function to get genes targeted either by species-specific or common peaks
get_gene = function(df,regulation_type){
    peak = copy(df)%>%
    lapply(function(x)x=x[
        numbpeaks_per_gene_same_tad>1 & numbpeaks_per_gene_same_tad<= 4 # i set this one for simplicity
        ][
            regulation %in% regulation_type
            ]
    )
    ## here range_keys come back in the right coordinates (also common peaks will be back as peaks not regions!)
    ## which is important for plotting
    peak = purrr::map2(peak,da_peaks,function(x,y)x[y,on='peakID',nomatch=0])

    if(regulation_type =='human' || regulation_type=='chimp'){
        all_genes = copy(peak)%>%lapply(function(x)levels(as.factor(x$gene)))
        peak_w_specific_genes = copy(peak)
        peak_w_specific_genes = purrr::map2(
            peak_w_specific_genes, all_genes,
            function(x,y)
            x=x[gene==sample(y,1,replace=F)]
        )
        return(peak_w_specific_genes)
    }else{
        common_genes = copy(peak)
        common_genes = merge(common_genes[[1]],common_genes[[2]],by='gene',allow.cartesian=TRUE)%>%dplyr::pull('gene')%>%unique()
        common_genes =sample(common_genes,1,replace=F)
        peak_w_common_gene = copy(peak)%>%lapply(
            function(x)x=x[gene%in%common_genes]
            )
        return(peak_w_common_gene)

    }
} 

common_specific_gene = get_gene(peaks_w_genes_same_tad,c('common_and_specific'))
common_specific_gene = Map(mutate,common_specific_gene,file=names(common_specific_gene))%>%rbindlist()
species_specific_genes = get_gene(peaks_w_genes_same_tad,c('human','chimp'))

## need to specify chimp and human marts
human_mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl",host='www.ensembl.org')

## decide how many bp to show around peak(s)
## this region must include the entire transcript 
## thus, first get gene start/end using this track generator
get_region = function(peak){
    df = copy(peak)
   gene_coord = BiomartGeneRegionTrack(
        chromosome = unique(as.character(df$seqnames)),
        filter=list(ensembl_gene_id=unique(as.character(df$EnsemblID))),
        biomart = human_mart,
        name='Target gene')
    gene_coord = as.data.table(gene_coord@range)[,start:=min(start)][,end:=max(end)][,c(..range_keys,'symbol')]%>%unique()
    
    region = copy(gene_coord)[
        ,start:=min(c(df$start,gene_coord$start))-20000
        ][
            ,end:=max(c(df$end,gene_coord$end))+20000
            ][
                ,c(..range_keys)
                ]%>%unique()

    return(region)
}

## then get the peak-surrounding base pairs
## these are determined as 1kb after gene end and 1kb before peak start
## note that this wont make the figure centered on the peak of interest (wonder who cares)
common_specific_region = get_region(common_specific_gene)
species_specific_region = lapply(species_specific_genes,function(x)get_region(x))

## function to get coordinate of TAD of interest
get_tad = function(peak_file,tad_file){
    tads =  copy(tad_file)
    peak = copy(peak_file)[,start:=min(start)][,end:=max(end)][,c(..range_keys)]%>%unique()
    tad_of_interest = foverlaps(peak,tads,type='within')%>%na.omit()%>%dplyr::select(c(all_of(range_keys)))
    tad_of_interest = tad_of_interest[,start:=min(start)][,end:=max(end)]%>%unique()
    return(tad_of_interest)

}
common_specific_tads = get_tad(common_specific_gene,all_tads[[2]])
species_specific_tads = purrr::map2(species_specific_genes,all_tads,function(x,y)get_tad(x,y))

## get chrom states for the TAD region across all iPSCs
## I deciced to combine all cell types
get_chromstate = function(region){
    grange = copy(region)
    chromstates = copy(ipsc_chromstate)[
        seqnames %in% unique(grange$seqnames)
        ][
            start>=unique(grange$start)
            ][
                end<=unique(grange$end)
                ][
                    ,state_numb:=as.numeric(gsub('\\_.*','',chrom_state))
                    ][
                        ,c('chrom_state','cell_type'):=NULL
                        ]%>%unique()%>%makeGRangesFromDataFrame(keep.extra.columns=T)
    return(chromstates)

}

common_specific_chromstate_tad = get_chromstate(common_specific_region)
species_specific_chromstate_tad = lapply(species_specific_region,function(x)get_chromstate(x))

## get conservation score for the region (peak region + surroundings)
get_cons_score = function(grange){
    tads = copy(grange)[,c(..range_keys)]
    tads=tads[,tadID:=paste('tad_',1:nrow(tads),sep='')]
    expanded_tads = copy(tads)[
        ,list(start = seq(start, end)),by=tadID
        ][
            ,seqnames:=unique(tads$seqnames)
            ][
                ,end:=start
                ]%>%makeGRangesFromDataFrame()
    
    phastCons = phastCons7way.UCSC.hg38

    conservation_scores = gscores(phastCons,expanded_tads)
    

   return(conservation_scores)
}

common_specific_cons_score = get_cons_score(common_specific_region)
species_specific_cons_score = lapply(species_specific_region,function(x)get_cons_score(x))

## test and polish single track before expanding for all 
ideoTrack =IdeogramTrack(
    genome = 'hg38',
    chromosome = unique(as.character(common_specific_tads$seqnames)),
    size = 2
)

gtrack = GenomeAxisTrack(
    chromosome = unique(as.character(common_specific_tads$seqnames)),
    start = min(common_specific_tads$start),
    end =  max(common_specific_tads$end),
    size = 5
)

chromstate_track = DataTrack(
        common_specific_chromstate_tad,
        type = "heatmap",
        name = "chrom state",
        gradient = chrom_state_colors,
        size = 5
)

biomart_track = BiomartGeneRegionTrack(
    chromosome = as.character(unique(common_specific_gene$seqnames)),
    start = as.numeric(common_specific_gene$start),
    end = as.numeric(common_specific_gene$end),
    biomart =human_mart,
    name='Target gene'
)

tad_track = AnnotationTrack(
    chromosome = as.character(unique(common_specific_tads$seqnames)),
    start = as.numeric(common_specific_tads$start),
    end = as.numeric(common_specific_tads$end),
    genome = 'hg38',
    fill='slategray',
    col='slategray',
    name = "TAD"
)

conservation_track = DataTrack(
        range = common_specific_cons_score,
        name = "phastCons7way",
        type='hist',
        col.histogram = "orangered2", 
        fill.histogram = "orangered2",
        size = 3
)

bigwigs = c(chimp_bigwigs,human_bigwigs)

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
bigwigs_track = get_track(bigwigs,common_specific_tads,c('deepskyblue1','dodgerblue3'))

# all_tracks = c(ideoTrack,gtrack,tad_track,biomart_track,conservation_track,bigwigs_track,chromstate_track)
all_tracks = c(biomart_track,conservation_track,bigwigs_track,chromstate_track)

common_region = copy(common_specific_gene)[peak_species=='common'][,start:=min(start)][,end:=max(end)][,c(..range_keys)]%>%unique()
sp_specific_peaks = copy(common_specific_gene)[peak_species!='common'][,c(..range_keys)]%>%unique()

highlight_region = rbind(common_region,sp_specific_peaks)
highlight_region_track = HighlightTrack(
    trackList = all_tracks,
    chromosome = unique(as.character(highlight_region$seqnames)),
    start = c(as.character(highlight_region$start)),
    end =  c(as.character(highlight_region$end)),
    col=c('coral2','cornsilk3','darkgoldenrod2'),
    fill=c('coral2','cornsilk3','darkgoldenrod2')
)

human_all_tracks = c(ideoTrack,gtrack,tad_track,highlight_region_track)

pdf(paste0(plot_dir,'human_common.pdf',sep=''),width = 10, height = 7)
plotTracks(
    human_all_tracks,
    from = common_specific_region$start,
    to =  common_specific_region$end,
    collapseTranscripts='longest', 
    transcriptAnnotation = 'symbol'
)
dev.off()



