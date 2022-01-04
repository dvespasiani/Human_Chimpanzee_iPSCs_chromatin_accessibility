library(dplyr)
library(data.table)
library(magrittr)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(viridis)
library(viridislite)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

fimo_dir <- './output/fimo_meme'
outplot_dir <- create_dir(plot_dir,'tfbs')

outfile_dir <- create_dir(outdir,'files/fimo_ctcf')

## read da peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))%>%setorderv('DA',1)%>%split(by='DA')

## read fimo results
fimo_outfiles <- list.files(fimo_dir,recursive=T,full.names=T,pattern='fimo.tsv')

fimo_out <- lapply(fimo_outfiles,function(x)x<-fread(x,sep='\t',header=T))
names(fimo_out) = c('da','non_da')

fimo_out <- lapply(
  fimo_out,function(x){
    x<-setnames(x,old=c('p-value','q-value'),new=c('p_value','q_value'))[
    q_value<=0.05
    ][
      ,log10_adjP:=-log10(q_value)
      ][
        ,c('q_value','p_value','motif_alt_id','strand'):=NULL
        ]
  }
)

## count numb peaks with at least 1 motif and prop of tot numb peaks 
purrr::map2(fimo_out,da_results,function(x,y){
  seqid_motif<-copy(x)[,c('sequence_name')]%>%unique()%>%nrow()
  seqid_allpeaks<-copy(y)[,c('peakID')]%>%unique()%>%nrow()

  prop<-round((seqid_motif/seqid_allpeaks)*100,2)
  return <- list(seqid_motif,seqid_allpeaks,prop)
  names(return)= c('numb_peaks_w_motif','numb_all_peaks','prop_peaks_w_motif')
  return(return)
})

## add cluster info
motif_cluster <- fread('/data/projects/punim0586/dvespasiani/Annotation_and_other_files/Motifs_clusters/motifs_clusters.txt',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

fimo_out_cluster <- copy(fimo_out)%>%lapply(
  function(x){
    x <- x[
    motif_cluster,on='motif_id',nomatch=0
    ][
      ,c('sequence_name','start','stop','motif_id','Name','log10_adjP')
    ]
  }
)

## count numb peaks associated to cluster 
purrr::map2(fimo_out_cluster,fimo_out,function(x,y){
  seqid_motif_cluster<-copy(x)[,c('sequence_name')]%>%unique()%>%nrow()
  seqid_motif<-copy(y)[,c('sequence_name')]%>%unique()%>%nrow()

  prop<-round((seqid_motif_cluster/seqid_motif)*100,2)
  return <- list(seqid_motif_cluster,seqid_motif,prop)
  names(return)= c('numb_peaks_w_motif_cluster','numb_peaks_w_motif','prop_peaks_w_motif_cluster')
  return(return)
})


##--------
## QCs
##--------
## number of motifs per peak (and avg)
numb_motif_peak <- copy(fimo_out_cluster)%>%lapply(function(x){
  x <- x[,c('sequence_name','Name')]%>%unique()
  x <- x[
    ,numb_motifs_per_peak:=.N,by=.(sequence_name)
    ][
      ,avg_motif_per_peak:=round(mean(numb_motifs_per_peak),2)
      ]
  return(x)
})

## distribution number motifs per peak
pdf(paste(outplot_dir,'qc_number_motifs_per_peak.pdf',sep=''),width = 10,height = 7)
ggplot(corr_numb_motif_peak_size,aes(x=numb_motifs_per_peak,fill=DA))+
geom_bar()+facet_wrap(DA~.,scale='free')
dev.off()

## check correlation with peak sizes
corr_numb_motif_peak_size <- purrr::map2(numb_motif_peak,da_results,function(x,y){
  y<-copy(y)[,width:=end-start][,peak_species:=NULL]
  x<-copy(x)[,peakID:=sequence_name][,sequence_name:=NULL][y,on='peakID',nomatch=0]
  return(x)
})%>%rbindlist()

cor.test(corr_numb_motif_peak_size$width,corr_numb_motif_peak_size$numb_motifs_per_peak)

# ## distance between consecutive motifs in same peak
# dist_btwn_near_motif <-  purrr::map2(numb_motif_peak,da_results,function(x,y){
#   y<-copy(y)[,width:=end-start][,c('start','end','peak_species'):=NULL]
#   x<-copy(x)[
#     ,peakID:=sequence_name
#     ][
#       ,c('motif_alt_id','strand','sequence_name'):=NULL
#       ][
#         y,on='peakID',nomatch=0
#         ]%>%setorderv(c('peakID','start','stop'),1)
#   x <- x[,dist_motifs:=start - data.table::shift(start),by=.(peakID)]%>%na.omit()
#   x <- x[,avg_dist_motifs:=round(mean(dist_motifs),1)]
#   return(x)
# })

# pdf(paste(outplot_dir,'qc_distance_between_motifs_per_peak.pdf',sep=''),width = 7,height = 7)
# ggplot(qc2[distance_motifs<=50],aes(x=distance_motifs,fill=DA))+
# geom_bar()+facet_wrap(DA~.,scale='free_y')
# dev.off()




# ## count peaks per cluster in each files
# fimo_results_cluster <- lapply(
#   fimo_results_cluster,function(x){
#     x<-x[
#       ,numb_motifs_cluster:=.N,by=.(Name)
#       ][
#         ,numb_motifs:=.N
#       ][
#         ,numb_motifs_noncluster:=numb_motifs-numb_motifs_cluster
#       ]
#   }
# )

## proceed with odds ratio of tfs per cluster
fimo_out_cluster_top_motif <- copy(fimo_out_cluster)%>%lapply(
  function(x)
  x<-x[
    ,.SD[which.max(log10_adjP)], by=.(sequence_name)
    ][
      ,numb_motifs_cluster:=.N,by=.(Name)
      ][
        ,numb_motifs:=.N
        ][
          ,numb_motifs_noncluster:=numb_motifs-numb_motifs_cluster
          ]
)

## look overlap motif clusters between non-da da peaks 
library(VennDiagram)

set_of_clusters <- copy(fimo_out_cluster_top_motif)%>%lapply(function(x)x=x[,Name]%>%unique())

venn.diagram(
    x = set_of_clusters,
    category.names = c('da','non_da'),
    filename = paste(outplot_dir,'set_of_motif_clusters_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col=da_palette,
    fill = c(alpha(da_palette[[1]],0.3), alpha(da_palette[[2]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-30, 30),
    cat.dist = c(0.1, 0.1),
    cat.fontfamily = "sans",
    cat.col = da_palette
)


counts_motifs_clusters <- copy(fimo_out_cluster_top_motif)%>%lapply(
  function(x){x<-x[,c('Name','numb_motifs_cluster','numb_motifs_noncluster')]%>%unique()}
)

nonda_vs_da_or <- calculate_or(counts_motifs_clusters[[2]],counts_motifs_clusters[[1]],'Name')%>%setnames(old='elements',new='Name')
nonda_vs_da_or <- nonda_vs_da_or[
  counts_motifs_clusters[[2]],on='Name',nomatch=0
  ][
    ,numb_motifs_noncluster:=NULL
    ][
      ,log10_numb_snps_tf:=log10(numb_motifs_cluster)
      ][
        ,log10_padj:= -log10(p.adj)
        ][
          ,log10_padj:=ifelse(log10_padj=='Inf',0,log10_padj)
          ][
            ,replace_p:=ifelse(p.signif!=' ' & log10_padj==0,'y','n')
]

## plot
tf_enrich_plot=function(x){
  df <- copy(x)[,'Log10 total number motifs per cluster':=log10_numb_snps_tf]
  gradient <- scale_colour_viridis(aes(`Log10 total number motifs per cluster`),option="inferno",discrete = F)
  text <- ifelse(!df$p.signif%in%' ',df$Name,'')
  
  ggplot(df,aes(x=log2(odds_ratio),log10_padj,label = text,col=log10_numb_snps_tf))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.2)+
    geom_text_repel(size = 5,color='black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines"),
                    max.overlaps = 50
    )+
    gradient+
    geom_hline(yintercept=1.3,linetype='dashed',size=0.5)+
    xlab('\n Log2 non da vs da OR \n')+
    ylab('-Log10 (adj-P)')+
    xlim(-8,8)+
    theme_bw()+ 
    theme(
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", colour = "black"),
        axis.line = element_blank())
}

pdf(paste(outplot_dir,'nonda_vs_da_or_fimo.pdf',sep=''),width = 10,height = 10)
tf_enrich_plot(nonda_vs_da_or)
dev.off()

##--------------
## CTCF peaks
##--------------
## counts number of CTCF peaks
count_motif <- function(tfbs_peaks,da_peaks,motif){

  peaks_w_motif <- copy(tfbs_peaks)[Name %like% motif]

  prop_peaks_w_motif <- copy(da_peaks)[
    ,peak_motif:=ifelse(peakID %in% peaks_w_motif$sequence_name,motif,'other')
    ][
      ,numb_motifs:=.N,by=.(peak_motif,DA)
      ][
        ,all_motifs:=.N,by=.(DA)
        ][
          ,prop_motifs:=round((numb_motifs/all_motifs)*100,2)
          ][
            peak_motif==motif
            ]
  return(prop_peaks_w_motif)
}

observed_ctcf_motifs <- purrr::map2(fimo_out_cluster,da_results,function(x,y)count_motif(x,y,'CTCF'))


## determine whether CTCF enrichment is larger than expected by chance

number_20bp_peaks <- lapply(observed_ctcf_motifs,function(p){
    p <- p%>%dplyr::select(-c(contains('peak')))
    p <- p[
      ,binned_width:=round(width/20,0)
      ][
        ,tot_numb_binned_width:=sum(binned_width)
        ]%>%dplyr::select('tot_numb_binned_width','DA',contains('motifs'))%>%unique()
})

cluster_names <- copy(fimo_out_cluster)%>%lapply(function(c)unique(c$Name))

## this takes too long 
motif_permutations <- purrr::map2(number_20bp_peaks,cluster_names,function(x,y){

    score_distribution = list()

  for(i in 1:10000){
    dt <- data.table(
      peakID = paste('peak_',1:x$tot_numb_binned_width,sep=''),
      Name = sample(y, x$tot_numb_binned_width, replace = T),
      DA = rep(x$DA,x$tot_numb_binned_width)
      )
    dt_tfbs_peaks <- copy(dt)%>%setnames(old='peakID',new='sequence_name')
    permuted_ctcf <- count_motif(dt_tfbs_peaks,dt,'CTCF')[,c('DA','numb_motifs','prop_motifs')]%>%unique()
    
    score_distribution[[i]] <- permuted_ctcf # add it to your list
  }
  
  score_distribution <- rbindlist(score_distribution)

  observed_value <- copy(x$numb_motifs)

  zscore <- (observed_value-mean(score_distribution$numb_motifs))/sd(score_distribution$numb_motifs)
  pvalue <- 2*pnorm(q=abs(zscore), lower.tail=FALSE)

  return <- score_distribution[,zscore:=zscore][,pvalue:=pvalue][,observed_value:=observed_value]
  return(score_distribution)
  }
)

test <- copy(motif_permutations)%>%rbindlist()
test_stat <- copy(test)[,c('DA','zscore','pvalue','observed_value')]%>%unique()
test_stat <-test_stat[,x:=observed_value][,y:=c(10,10)
  ][
    ,stat:=paste('zscore = ',round(zscore,2),'\n', '-log10 pval = ',round(-log10(pvalue),2),sep='')
]

## permutation plot
pdf(paste0(outplot_dir,'permuted_number_ctcfs.pdf',sep=''),width = 15, height = 10)
p <- ggplot(test,aes(x=numb_motifs,fill=DA))+
geom_histogram(binwidth=4)+
geom_vline(data = test, aes(xintercept=observed_value),linetype='dashed',color='red')+
facet_wrap(DA~.,ncol=2,scale='free')+
scale_fill_manual(values=da_palette)+
xlab('permutation number CTCFs')+ylab('counts')+
theme_bw()+
theme(legend.position='none')

p + geom_text(data = test_stat,aes(x=x-120,y=y,label=stat),size=5)
dev.off()


## export CTCF peaks 
ctcf_peaks <- copy(fimo_out_cluster)%>%lapply(
  function(x){ 
  x <- x[
    ,c('sequence_name','Name','start','stop')
    ][
      Name %like%'CTCF'
      ][
        ,start:=as.numeric(start)
      ]%>%
      setnames(old=c(1,3),new=c('peakID','begin'))%>%
      setorderv('peakID')
  }
)

ctcf_peaks <- purrr::map2(
  ctcf_peaks,da_results,
  function(x,y){
    z<-merge(x,y,by='peakID') 
    return(z)
    }
)%>%rbindlist()

ctcf_peaks <- ctcf_peaks[
      ,.SD[which.min(begin)], by=.(peakID)
      ][
        ,.SD[which.max(stop)], by=.(peakID)
        ][
          ,motif_start:=start+begin
          ][
            ,motif_end:=start+stop
] ## checked that these are the same of peaks with ctcf  

## get dna seq 
library(BSgenome.Hsapiens.UCSC.hg38)

ctcf_peaks<-ctcf_peaks[,c('begin','stop','start','end'):=NULL]
ctcf_peaks <- ctcf_peaks%>%
makeGRangesFromDataFrame(
  start.field="motif_start",end.field="motif_end",keep.extra.columns=T)%>%
  resize(width=20, fix = "center")%>%
as.data.table()

ctcf_peaks <- ctcf_peaks[
        ,sequence:=as.character(getSeq(Hsapiens, seqnames,start, end))
        ][
          ,c('width','strand'):=NULL
]

write.table(ctcf_peaks,file=paste(outfile_dir,'fimo_ctcf_peaks.txt',sep=''),sep='\t',row.names = F,quote=F,col.names=T)