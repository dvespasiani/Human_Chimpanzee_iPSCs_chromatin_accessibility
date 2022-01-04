
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(phastCons7way.UCSC.hg38)
library(GenomicScores)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

# ctcf_dir <- './output/homer/homer_output/ctcf/'
outplot_dir <- create_dir(plot_dir,'ctcf')
peakDir = './output/DA/peaks/'
tads_dir =  './output/TADs/'
files_dir <- './output/files/'

##-----------------------
## read deepBind results
##-----------------------
## read deepbind input with all the info 
deepbind_output_file <- dir(paste(files_dir,'deepbind_predictions',sep=''),recursive=F,full.names=T,pattern='fimo')
ctcf_peak_file <- dir(paste(files_dir,'fimo_ctcf',sep=''),recursive=F,full.names=T,pattern='fimo')

ctcf_peaks <- fread(ctcf_peak_file,header=T,sep='\t')
ctcf_affinities <- fread(deepbind_output_file,sep='\t',header=T)%>%dplyr::select(c(1:2),contains('D00328.003'))%>%setnames(old=c(2,3),new=c('sequence','prediction'))

ctcf_peaks <- ctcf_peaks[
    ctcf_affinities,on=c('peakID','sequence'),nomatch=0
    ][
        ,peak_species:=ifelse(peak_species=='non_da',"common",ifelse(peak_species=='da_chimp','chimp','human'))
]
setkeyv(ctcf_peaks,range_keys)

ctcf_peaks <- setorderv(ctcf_peaks,'peak_species')%>%split(by='peak_species')

##-------------------------------------------------
## look for overlaps between ctcf at TAD boundaries
##--------------------------------------------------
orth_tads <- read_tads('human_chimp_orth_tads.txt')
chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL][species=='Chimp']
human_tads_hg38 <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]

chimp_tads_hg38 <- copy(orth_tads)[species=='Chimp']
colnames(chimp_tads_hg38)[4:6] = range_keys
chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'species'),nomatch=0][
  ,c(range_keys):=NULL
]%>%setnames(old=c(1:3),new=range_keys)

common_tads_hg38 <- copy(human_tads_hg38)[species=='Common']

all_tads <- list(chimp_tads_hg38,common_tads_hg38,human_tads_hg38[species!='Common'])

all_tads <- lapply(all_tads,function(x){
  x <- setorderv(x,range_keys,1)
  x<-x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]
  setkeyv(x,range_keys)
  return(x)
  }
) 
names(all_tads) = species_names

## define TAD boundaries
boundary =  500000
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
ctcf_boundary_overlap <- purrr::map2(ctcf_peaks,tad_boundaries,function(x,y)foverlaps(x,y,type='within')%>%na.omit())
ctcf_boundary_overlap <- lapply(ctcf_boundary_overlap,
    function(x)x=x[
        ,distance_ctcf_tad:= ifelse(boundary=='upstream',i.end-end,i.start-start)
    ]
)

## plot distribution distances between CTCF and TAD boundaries
## if peak is associated with multiple tads keep closest one
distance_ctcf_boundaries <- copy(ctcf_boundary_overlap)%>%
lapply(function(x)x=x[
    ,c('peakID','prediction','species','distance_ctcf_tad')
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
distance_ctcf_boundaries <- Map(mutate,distance_ctcf_boundaries,peak_species=species_names)%>%rbindlist()

pdf(paste(outplot_dir,'fimo_ctcf_tad_distances.pdf',sep=''),width=10,height = 7)
ggplot(distance_ctcf_boundaries, aes(x=distance_ctcf_tad,fill=peak_species)) +
geom_density(alpha=0.5) +ylab('Density CTCF peaks')+ xlab('CTCF-TAD Distance')+
theme(
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
)
dev.off()

## get some % CTCF peaks within tab boundaries 
numb_ctcfs_associated <- copy(distance_ctcf_boundaries)[associated=='tad_associated'][,numb_ctcfs:=.N,by=.(peak_species)][,c('numb_ctcfs','peak_species')]%>%unique()%>%split(by='peak_species')
prop_ctcfs_associated <- purrr::map2(ctcf_peaks,numb_ctcfs_associated,function(x,y){
    x<-copy(x)
    prop<-round((y$numb_ctcfs/nrow(x))*100,2)
    return(prop)
})

##-------------------------------------------------------------------
## work only with common (non-da) peaks with CTCF binding sites 
##-------------------------------------------------------------------
## Get conservation score 
## of common/species specific CTCF peaks near TADs 

common_ctcf_peaks <- copy(ctcf_peaks[[2]])[distance_ctcf_boundaries,on=c('peakID','peak_species','prediction'),nomatch=0]%>%dplyr::select(-c(contains('dist')))

## get average conservation score
phastCons <- phastCons7way.UCSC.hg38

common_ctcf_peaks_conScore <- copy(common_ctcf_peaks)%>%makeGRangesFromDataFrame(keep.extra.columns=T)
common_ctcf_peaks_conScore <- gscores(phastCons,common_ctcf_peaks_conScore)%>%as.data.table()

common_ctcf_peaks_conScore <- common_ctcf_peaks_conScore[, default := nafill(default,fill=0)]

## plot function
plot_ctcf_tad = function(x,column,ylab,xlab){
    df=copy(x)[
        ,column_to_plot:=column
    ]
    p = ggplot(df, aes(x=associated, y=column_to_plot,fill=associated)) +
    geom_violin(trim=T,scale = "width")+
    geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+ 
    ylab(ylab)+xlab(xlab)+
    facet_wrap(peak_species~.,ncol=3)+
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
pdf(paste(outplot_dir,'fimo_ctcf_peaks_PhastCons_score.pdf',sep=''),width= 10,height = 8)
plot_ctcf_tad(common_ctcf_peaks_conScore,common_ctcf_peaks_conScore$default,ylab='PhastCons score',xlab=' ') 
dev.off()

## plot Deepbind CTCF affinities
pdf(paste(outplot_dir,'fimo_ctcf_peaks_Deepbind_affinities.pdf',sep=''),width= 10,height = 8)
plot_ctcf_tad(common_ctcf_peaks_conScore,common_ctcf_peaks_conScore$prediction,ylab='CTCF DeepBind predicted affinity',xlab=' ') 
dev.off()

# ## plot Motifscores
# pdf(paste(plot_dir,'ctcf_peaks_Motifscores.pdf',sep=''),width= 10,height = 8)
# plot_ctcf_tad(ctcf_peak_affinities,ctcf_peak_affinities$MotifScore,ylab='CTCF HOMER motif score',xlab=' ') 
# dev.off()






