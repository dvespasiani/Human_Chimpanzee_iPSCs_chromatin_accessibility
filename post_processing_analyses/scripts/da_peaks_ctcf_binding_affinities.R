## use this script to analyse the predicted CTCF binding profile 
## in particular compare the CTCF binding affinities between peaks containing the ctcf bs and all the genome-wide bs fo ctcf
## then do the analysis by comparing the affinities btwn tad associated and non-tad associated ctcf
## then look if there is an enrichment of CTCF bs near DC TADs 

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(liftOver)
library(rtracklayer)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

ctcf_dir = './post_processing_analyses/output/homer/homer_output/ctcf/'
plot_dir = './post_processing_analyses/output/plots/ctcf/'
peakDir = './post_processing_analyses/output/DA/peaks/'
tads_dir =  './post_processing_analyses/output/TADs/'
tmp_files_dir = './post_processing_analyses/output/temp_files/'

##-----------------------
## read deepBind results
##-----------------------
ctcf_affinities = list.files(
    paste(ctcf_dir,'deepbind_predictions/',sep=''),recursive=F,full.names=T,pattern='affinity')%>%
    lapply(
        function(x)fread(x,sep='\t',header=F,col.names=c('sequence','prediction'))[-1,] ## remove 1st row containing module name only
) 
names(ctcf_affinities) =  gsub("\\..*","",list.files(paste(ctcf_dir,'deepbind_predictions/',sep=''),recursive=F,full.names=F,pattern='affinity'))

## this files contains the peakIDs for each predicted CTCF binding affinity
ctcf_peakID = fread(paste(ctcf_dir,'ctcf_peak_sequences/','ctcf_match_file.txt',sep=''),sep='\t',header=T)

## assign metadata to predictions
## beware that u will have duplicates because homer returns multiple CTCFs for a single peak which are offset by few bp from each other
## to counteract this simply take the sequence with the highest absolute binding affinity for each single peak
ctcf_affinities_peakIDs = copy(ctcf_affinities[[3]])[
    ctcf_peakID,on='sequence',nomatch=0
    ][
        ,prediction:=as.numeric(prediction)
        ][ 
            ,.SD[which.max(abs(prediction))], by=.(peakID,file)
]%>%split(by='file')

## Add genomic coordinates to CTCF-peaks
ctcf_peak_coords = list.files(paste(ctcf_dir,'ctcf_peaks/',sep=''),recursive=F,full.names=T)%>%lapply(function(x)fread(x,sep='\t',header=T)%>%unique())
names(ctcf_peak_coords) =  gsub("\\..*","",list.files(paste(ctcf_dir,'ctcf_peaks/',sep=''),recursive=F,full.names=F))

ctcf_peak_coords =
lapply(ctcf_peak_coords,function(z)z=z[
        ,peak_center:= round((end-start)/2,0)
        ][
            ,ctcf_peak_start:=start+(peak_center+Offset) # this returns the coord of the start of ctcf peak
            ][
                ,ctcf_peak_end:=ctcf_peak_start+nchar(sequence)  # this returns the coord of the end of ctcf peak (based on length of motif sequence)
                ][
                    ,c('peak_center','sequence'):=NULL

        ]
)

## add motif score and strand information
ctcf_score_strand = list.files(paste(ctcf_dir,'target_sequences/',sep=''),recursive=F,full.names=T)%>%
lapply(
    function(x)fread(x,sep='\t',header=T,drop=c('Sequence','Offset','Motif Name'))[
        ,.SD[which.max(abs(MotifScore))],by=.(`FASTA ID`) ## if a ctcf has two scores (1 per strand) keep the highest one
        ]%>%setnames(old='FASTA ID',new='peakID')%>%unique() 
) 

##-------------------------------------------------------------
## create single ctcf file with all the info gathered till now
##-------------------------------------------------------------
ctcf_peaks = purrr::map2(ctcf_peak_coords,ctcf_affinities_peakIDs,function(x,y)x[y,on='peakID',nomatch=0])
ctcf_peaks = purrr::map2(ctcf_peaks,ctcf_score_strand,function(x,y)x[y,on='peakID',nomatch=0])
ctcf_peaks =lapply(
    ctcf_peaks,function(x)x=x[
        ,c('sequence','file'):=NULL
        ][
            , .SD[which.min(abs(Offset))],by=.(peakID,MotifScore,prediction) ## if a peak has multiple CTCFs with same prediction and motifscore, keep the one closest to peak center
            ][
                ,numbpeaks:=.N
                ]
)
lapply(ctcf_peaks,function(x)setkeyv(x,c('seqnames','ctcf_peak_start','ctcf_peak_end')))

## write this file for Gviz 
ctcf_out_filenames = paste(tmp_files_dir,names(ctcf_peaks),'_allinfo.txt',sep='')
mapply(write.table,ctcf_peaks, file = ctcf_out_filenames,col.names = T, row.names = F, sep = "\t", quote = F)

# ##-----------------------------------------------
# ## Compare the predicted binding affinities 
# ## and Homer MotifScores differ between species
# ##-----------------------------------------------
# all_ctcf_peaks = copy(ctcf_peaks)
# all_ctcf_peaks = Map(mutate,all_ctcf_peaks,species=species_names)%>%rbindlist()

# stat_test_ctcf_affinities = 
# compare_means(
#   prediction~species,
#   all_ctcf_peaks,
#   method = "wilcox.test",
#   paired = FALSE,
#   ref.group = 'common',
#   p.adjust.method = "fdr"
# ) %>% as.data.table()%>%setnames(old='group2',new='species')

# # stat_test_ctcf_affinities =Map(mutate,stat_test_ctcf_affinities,file=species_names[-2])%>%rbindlist()

# stat_test_MotifScore = 
# compare_means(
#   MotifScore~species,
#   all_ctcf_peaks,
#   method = "wilcox.test",
#   paired = FALSE,
#   ref.group = 'common',
#   p.adjust.method = "fdr"
# ) %>% as.data.table()%>%setnames(old='group2',new='species')

# # stat_test_MotifScore =Map(mutate,stat_test_MotifScore,file=species_names[-2])%>%rbindlist()

# # all_ctcf_peaks = Map(mutate,all_ctcf_peaks,file=species_names[-2])%>%rbindlist()

# ##plot
# pdf(paste(plot_dir,'ctcf_binding_affinities.pdf',sep=''),width=10,height = 7)
# ggplot(all_ctcf_peaks, aes(x=species, y=prediction,fill=species)) +
# geom_boxplot(notch=T) +ylab('binding affinity prediction score')+
# geom_text(data=stat_test_ctcf_affinities, aes(x=species, y=max(all_ctcf_peaks$prediction+0.01), label=p.signif), col='black', size=7)+
# theme(
#     axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#     )
# dev.off()

# pdf(paste(plot_dir,'ctcf_motifscores.pdf',sep=''),width=10,height = 7)
# ggplot(all_ctcf_peaks, aes(x=species, y=MotifScore,fill=species)) +
# geom_boxplot(notch=T) +ylab('Homer predicted MotifScore')+
# geom_text(data=stat_test_ctcf_affinities, aes(x=species, y=max(all_ctcf_peaks$MotifScore+0.01), label=p.signif), col='black', size=7)+
# theme(
#     axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#     )
# dev.off()

##-------------------------------------------------
## look for overlaps between ctcf at TAD boundaries
##--------------------------------------------------
## do it separately for TADs and TAD boundaries 
## for CTCFs overlapping TAD boundaries 

## read TADs
human_tads = read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]
chimp_tads = read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

common_hg38_tads = copy(human_tads)[species %in% 'Common'] ## because numb of TADs is the same and because common peaks are in hg38
human_hg38_tads = copy(human_tads)[species %in% 'Human']
chimp_pantro5_tads = copy(chimp_tads)[species %in% 'Chimp']

all_tads = list(chimp_pantro5_tads,common_hg38_tads,human_hg38_tads)
lapply(all_tads,function(x)setkeyv(x,range_keys))
names(all_tads) = species_names

## define TAD boundaries
boundary =  1000000
boundaries = copy(all_tads)%>%
lapply(
    function(x)x=x[
    ,tadID:=paste('tad_',seq(1:nrow(x)),sep='')
    ][
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
    function(x)x=x%>%dplyr::select(c('seqnames',contains('upstream'),'species','tadID'))%>%
    setnames(old=c('upstream_start','upstream_end'),new=range_keys[-1])%>%
    mutate('boundary'='upstream')
)
downstream_boundaries = copy(boundaries)%>%
lapply(
    function(x)x=x%>%dplyr::select(c('seqnames',contains('downstream'),'species','tadID'))%>%
    setnames(old=c('downstream_start','downstream_end'),new=range_keys[c(2,3)])%>%
    mutate('boundary'='downstream')
)

tad_boundaries = purrr::map2(upstream_boundaries,downstream_boundaries,rbind)
lapply(tad_boundaries,function(x)setkeyv(x,range_keys))

## look at ctcf overlap at tad boundaries
ctcf_boundary_overlap = purrr::map2(ctcf_peaks,tad_boundaries,function(x,y)foverlaps(x,y,type='within')%>%na.omit())
ctcf_boundary_overlap = lapply(ctcf_boundary_overlap,
    function(x)x=x[
        ,distance_ctcf_tad:= ifelse(boundary=='upstream',ctcf_peak_end-end,ctcf_peak_start-start)
    ]
)

## plot distribution distances between CTCF and TAD boundaries
## if peak is associated with multiple tads keep closest one
distance_ctcf_boundaries = copy(ctcf_boundary_overlap)%>%
lapply(function(x)x=x[
    ,c('peakID','prediction','species','distance_ctcf_tad','MotifScore','Strand')
    ][
        , .SD[which.min(abs(distance_ctcf_tad))], by=.(peakID,species)
        ][
            ,log10dist:=log10(abs(distance_ctcf_tad))
            ][
                ,log10dist:=ifelse(distance_ctcf_tad<0, -log10dist,log10dist)
                ][
                    ,associated:=ifelse(abs(distance_ctcf_tad)<=50000,'tad_associated','non_associated')
                    ]
)
distance_ctcf_boundaries = Map(mutate,distance_ctcf_boundaries,species=species_names)%>%rbindlist()

pdf(paste(plot_dir,'ctcf_tad_distances.pdf',sep=''),width=10,height = 7)
ggplot(distance_ctcf_boundaries, aes(x=distance_ctcf_tad,fill=species)) +
geom_density(alpha=0.5) +ylab('Density CTCF peaks')+ xlab('CTCF-TAD Distance')+
theme(
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()

## get some % CTCF peaks within tab boundaries 
numb_ctcf_peaks = copy(ctcf_peaks)%>%lapply(function(x)x[,c(..range_keys)]%>%unique()%>%nrow())
numb_ctcf_peaks_tad_overlap_1mb = copy(ctcf_boundary_overlap)%>%lapply(function(x)x[,c('seqnames','i.start','i.end')]%>%unique()%>%nrow())
numb_ctcf_peaks_tad_overlap_50kb = copy(ctcf_boundary_overlap)%>%lapply(function(x)x[abs(distance_ctcf_tad)<=50000][,c('seqnames','i.start','i.end')]%>%unique()%>%nrow())

prop_ctcf_peaks_tad_overlap_1mb = purrr::map2(numb_ctcf_peaks_tad_overlap_1mb,numb_ctcf_peaks,function(x,y)x/y*100)
prop_ctcf_peaks_tad_overlap_50kb = purrr::map2(numb_ctcf_peaks_tad_overlap_50kb,numb_ctcf_peaks,function(x,y)x/y*100)

##------------------------------------------------
## Get conservation score 
## of common/species specific CTCF peaks near TADs 
##------------------------------------------------
ctcf_peaks_cons_score = copy(ctcf_peaks)%>%lapply(function(x)x[,c('seqnames','ctcf_peak_start','ctcf_peak_end','peakID')]%>%unique())

ctcf_chimp_peaks_pantro5 = copy(ctcf_peaks_cons_score[[1]])
colnames(ctcf_chimp_peaks_pantro5)[1:3] = range_keys
ctcf_chimp_peaks_hg38 = convert_coord(ctcf_chimp_peaks_pantro5,'panTro5ToHg38.over.chain')

all_ctcf_peaks_hg38 = list(ctcf_chimp_peaks_hg38,ctcf_peaks_cons_score[[2]],ctcf_peaks_cons_score[[3]])%>%
lapply(function(x)makeGRangesFromDataFrame(x,keep.extra.columns=T))
names(all_ctcf_peaks_hg38) = species_names

## get average conservation score
phastCons = phastCons7way.UCSC.hg38
ctcf_peak_scores = lapply(all_ctcf_peaks_hg38,function(x)gscores(phastCons, x)%>%as.data.table())

ctcf_boundary_overlap_50kb = lapply(ctcf_boundary_overlap,function(x)x=x[abs(distance_ctcf_tad)<=50000])

ctcf_peak_scores = purrr::map2(ctcf_peak_scores,ctcf_boundary_overlap_50kb,function(x,y)
x=x[
    ,tad_overlap:=ifelse(
        peakID %in% y$peakID,'yes','no'
        )
        ][
            , default := nafill(default,fill=0)
            ]
)

ctcf_peak_scores = Map(mutate,ctcf_peak_scores,species=names(ctcf_peak_scores))%>%rbindlist()

## plot function
plot_ctcf_tad = function(x,column,ylab,xlab){
    df=copy(x)[
        ,column_to_plot:=column
    ]
    p = ggplot(df, aes(x=tad_overlap, y=column_to_plot,fill=tad_overlap)) +
    geom_boxplot(notch=T) +
    ylab(ylab)+xlab(xlab)+
    facet_wrap(species~.,ncol=3)+
    stat_compare_means(
        method = "wilcox.test",
        # label = "p.signif",
        label.y = max(df$column_to_plot+0.1),
        size=5
        )+
    theme(
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
    return(p)
}

## plot cons scores
pdf(paste(plot_dir,'ctcf_peaks_PhastCons_score.pdf',sep=''),width= 10,height = 8)
plot_ctcf_tad(ctcf_peak_scores,ctcf_peak_scores$default,ylab='PhastCons score',xlab=' ') 
dev.off()

## now plot motif binding affinities for CTCFs within 50kb
ctcf_peak_affinities = purrr::map2(
    ctcf_peaks,ctcf_boundary_overlap_50kb,function(x,y)
    x=x[
        ,tad_overlap:=ifelse(peakID %in% y$peakID,'yes','no')
        ]
)
ctcf_peak_affinities = Map(mutate,ctcf_peak_affinities,species=names(ctcf_peak_affinities))%>%rbindlist()

## plot Deepbind CTCF affinities
pdf(paste(plot_dir,'ctcf_peaks_Deepbind_affinities.pdf',sep=''),width= 10,height = 8)
plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
dev.off()

## plot Motifscores
pdf(paste(plot_dir,'ctcf_peaks_Motifscores.pdf',sep=''),width= 10,height = 8)
plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$MotifScore,ylab='CTCF HOMER motif score',xlab=' ') 
dev.off()

