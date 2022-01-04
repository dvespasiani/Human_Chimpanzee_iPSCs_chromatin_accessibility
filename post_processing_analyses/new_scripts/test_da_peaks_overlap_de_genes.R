## script used to look for overlap between DA peaks and DE genes
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(VennDiagram)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'rna_seq')
target_genes_dir <- './output/files/GO_enrich'
# tmp_files_dir <- './output/temp_files/'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'logFC','DA','peakID','peak_species'))
setkeyv(da_results,range_keys)
da_results <- da_results[
  ,c(range_keys):=NULL
]%>%setorderv('DA',1)%>%split(by='DA')

## read target genes
target_genes <- list.files(target_genes_dir,full.names= T,recursive=F,pattern='target')%>%lapply(
  function(x)fread(x,sep='\t',header=T))%>%rbindlist()%>%split(by='DA')
  
peaks_w_genes <- purrr::map2(da_results,target_genes,function(x,y)x[y,on=c('peakID','peak_species','DA'),nomatch=0])%>%rbindlist()

##-----------------
## read DE genes 
##-----------------
## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
de_genes <- fread(
  "../rna_seq/de_output/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",
  sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B')
  )[
    ,logFC:=-logFC
    ][
        ,DE:=ifelse(adj.P.Val <= 0.01,'de','non_de')
        ][
          ,c('P.Value','t','B','genes'):=NULL
]

## keep genes with ensembl ids from v 86
ensembl_86_ids <- fread("../rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

peaks_w_genes_ensembl <- copy(peaks_w_genes)[
  ensembl_86_ids,on='gene',nomatch=0
  ][
    de_genes,on='EnsemblID',nomatch=0
    ][
      ,c('width','strand','go_signif'):=NULL
]%>%unique()


human_common_peaks <- copy(peaks_w_genes_ensembl[peak_species!='chimp'])
chimp_common_peaks <- copy(peaks_w_genes_ensembl[peak_species!='human'])

all_peaks <- list(chimp_common_peaks,human_common_peaks)
names(all_peaks) = species_names[-2]
lapply(all_peaks,function(x)setkeyv(x,range_keys))

##--------------------
## peaks within TADs
##--------------------
## read TADs
orth_tads <- read_tads('human_chimp_orth_tads.txt')%>%setnames(old='species','tad_species')
chimp_tads_pantro5 <- read_tads('pantro5_tads.bed')[,tad_species:=ifelse(label==1,'Common',species)][,c('label','species'):=NULL]
human_tads_hg38 <- read_tads('hg38_tads.bed')[,tad_species:=ifelse(label==1,'Common',species)][,c('label','species'):=NULL]

chimp_tads_hg38 <- copy(orth_tads)[tad_species=='Chimp']
colnames(chimp_tads_hg38)[4:6] = range_keys
chimp_tads_hg38 <- chimp_tads_hg38[chimp_tads_pantro5,on=c(range_keys,'tad_species'),nomatch=0][
  ,c(range_keys):=NULL
]%>%setnames(old=c(1:3),new=range_keys)

common_tads <- copy(human_tads_hg38)[tad_species=='Common']
human_tads <- copy(human_tads_hg38)[tad_species!='Common']
chimp_tads <- copy(chimp_tads_hg38)

all_tads <- list(rbind(chimp_tads,common_tads),rbind(human_tads,common_tads))
names(all_tads) = species_names[-2]
lapply(all_tads,function(x)setkeyv(x,range_keys))

## overlaps
peaks_in_tads <- purrr::map2(all_peaks,all_tads,function(x,y){
  x<-copy(x)
  y<-copy(y)
  z<-foverlaps(x,y,type='within')%>%na.omit()
  # z <- z[
  #     ,gene_start := ifelse(distTSS<0,i.start+distTSS,i.end+distTSS)
  #     ][
  #         ,same_tad := ifelse((gene_start-start>=0),'y','n')
  #         ]
  return(z)
  }
)

## test enrichment of DA peaks in species-specific TADs 
qc <-copy(peaks_in_tads)%>%lapply(function(x){
  x<-x[,numb_peak_type_in_tad_type:=.N,by=.(tad_species,DA)][,c('numb_peak_type_in_tad_type','tad_species','DA')]%>%unique()
})

## test whether there are more DE genes in human/chimp tads compared to common ones 


## keep peaks w genes in same tad
peaks_w_genes_same_tad <- copy(peaks_w_genes_tads)[
  same_tad=='y'
  ][
    ,c('same_tad','gene_start'):=NULL
]

## calculate proportion of peaks with genes in the same tad
prop_peaks_w_genes_same_tad <- length(unique(peaks_w_genes_ensembl[peakID %in% peaks_w_genes_same_tad$peakID]$peakID))/length(unique(peaks_w_genes_ensembl$peakID))*100
## 83.3%

## if peak is associated to multiple tads keep the widest one 
peaks_w_genes_same_tad <- peaks_w_genes_same_tad[
  ,tad_width:=end-start
  ][
  ,.SD[which.max(tad_width)], by=.(peakID)
]


## check how many peaks each gene has associated and stratify that by de/non-de  
peaks_w_genes_same_tad <- peaks_w_genes_same_tad[,numb_peaks_per_gene:=.N,by=.(gene,DE)]

pdf(paste0(outplot_dir,'numb_peaks_per_gene.pdf',sep=''),width = 10, height = 8)
ggplot(peaks_w_genes_same_tad,aes(x=DE,y=numb_peaks_per_gene,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

## compare proportion da peaks per de/non de gene
peaks_w_genes_same_tad <- peaks_w_genes_same_tad[
  ,numb_peak_type_per_gene:=.N,by=.(gene,DE,DA)
  ][
    ,prop_peaktype_per_gene:=numb_peak_type_per_gene/numb_peaks_per_gene
]

pdf(paste0(outplot_dir,'prop_da_peaks_per_gene.pdf',sep=''),width = 10, height = 8)
df <- copy(peaks_w_genes_same_tad)[numb_peaks_per_gene>1][,c('prop_peaktype_per_gene','gene','DA','DE')][DA=='da']%>%unique()
ggplot(df,aes(x=DE,y=prop_peaktype_per_gene,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
ylab('proportion DA peaks per gene')+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

##-------
## QCs
##-------
## compare distribution distances between peaks associtated to same gene for de non de genes
qc <- copy(peaks_w_genes_same_tad)[
  ,min_start:=min(start),by=.(gene)
  ][
    ,max_end:=max(end),by=.(gene)
    ][
      ,regulatory_region:=log10(max_end-min_start)
]

pdf(paste0(outplot_dir,'span_peaks_regulatory_region_per_gene.pdf',sep=''),width = 10, height = 8)
ggplot(qc,aes(x=DE,y=regulatory_region,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
ylab('log10 bp regulatory region')+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

## difference btw max and min atac logfc per gene
qc <- qc[,diff_logFC:=ifelse(numb_peaks_per_gene>1,round(max(logFC)-min(logFC),2),NA),by=.(gene)]%>%na.omit()

pdf(paste0(outplot_dir,'diff_atac_logFC_per_gene.pdf',sep=''),width = 10, height = 8)
ggplot(qc,aes(x=DE,y=diff_logFC,fill=DE))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_de',
  size=5
  )
dev.off()

##--------------------------------------------------------------------------------------------------
## for each gene keep the peak with highest logFC (this would bias towards DA peaks)
##--------------------------------------------------------------------------------------------------
one_peak_per_gene <- copy(peaks_w_genes_ensembl)[
  ,.SD[which.max(abs(logFC))], by=.(EnsemblID)
]

# one_peak_per_gene <- copy(peaks_w_genes_same_tad)[
#   ,peak:=ifelse(DA=='non_da',0,1),by=.(gene)
#   ][
#     ,da_peaks:=sum(peak),by=.(gene)
# ]
# gene_da_peaks <- copy(one_peak_per_gene)[da_peaks>0][DA=='da'][
#       ,.SD[which.max(logFC)], by=.(EnsemblID)
# ]
# gene_nonda_peaks <- copy(one_peak_per_gene)[da_peaks==0][
#       ,.SD[which.max(logFC)], by=.(EnsemblID)
# ]

# one_peak_per_gene <- rbind(gene_da_peaks,gene_nonda_peaks)


##------------------------------------------------------
## look at avg expression genes by peaks (DA vs non-DA)
##------------------------------------------------------
all_genes_expr <- copy(de_genes)[,c('AveExpr')][,DA:='all_genes']

comparisons = list(
  c('non_da','da'),
  c('da','all_genes'),
  c('non_da','all_genes')
)

plot_de_genes <- function(df,xaxis,yaxis,ylab){
  df_plot <- copy(df)[
    ,xcol:=xaxis
    ][
      ,ycol:=yaxis
      ]

  p <- ggplot(df_plot,aes(x=xcol,y=ycol,fill=xcol))+
  geom_violin(trim=T,scale = "width")+
  geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
  xlab(' ')+ylab(ylab)+
  stat_compare_means(
  method = "wilcox.test",
  comparisons = comparisons,
  size=5
  )+
    theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
  return(p)
}

pdf(paste0(outplot_dir,'avg_expr_genes_by_peak_significance.pdf',sep=''),width = 10, height = 8)
df <- rbind(one_peak_per_gene[,c('AveExpr','DA')],all_genes_expr)
plot_de_genes(df,df$DA,df$AveExpr,'AveExpr')
dev.off()


## get loeuf score for these genes
gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
colnames(gnomad_loeuf)[1] = 'EnsemblID'

## add loeuf score to genes
one_peak_per_gene_loeuf <- copy(one_peak_per_gene)[gnomad_loeuf,on='EnsemblID',nomatch=0]

all_genes_loeuf <- copy(gnomad_loeuf)[,DA:='all_genes']
peak_genes_loeuf <- copy(one_peak_per_gene_loeuf)%>%dplyr::select(c(all_of(colnames(all_genes_loeuf))))

pdf(paste0(outplot_dir,'loeuf_score_genes_by_peak.pdf',sep=''),width = 10, height = 8)
df <- rbind(all_genes_loeuf,peak_genes_loeuf)
plot_de_genes(df,df$DA,df$oe_lof_upper,'oe_lof_upper')
dev.off()


##-----------------------------------------
## ovelap between DE genes and DA peaks
##-----------------------------------------

significant_da_genes <- copy(one_peak_per_gene)[DA=='da'][,EnsemblID]%>%unique()
significant_de_genes <- copy(one_peak_per_gene)[DE=='de']$EnsemblID%>%unique()

da_de_overlap <-  list(significant_da_genes,significant_de_genes)
venn.diagram(
    x = da_de_overlap,
    category.names = c("DA","DE"),
    filename = paste(outplot_dir,'da_de_overlap.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
   col=c('#fde725ff',"#440154ff"),
    fill = c(alpha('#fde725ff',0.3),alpha("#440154ff",0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = c('#fde725ff',"#440154ff")
)

## test significance overlap
numb_da_de_genes <- Reduce(intersect,da_de_overlap)%>%length()
numb_nonda_de_genes <- copy(one_peak_per_gene)[DA!='da'][DE=='de'][,EnsemblID]%>%unique()%>%length()

numb_da_nonde_genes <- copy(one_peak_per_gene)[DA=='da'][DE!='de'][,EnsemblID]%>%unique()%>%length()
numb_nonda_nonde_genes <- copy(one_peak_per_gene)[DA!='da'][DE!='de'][,EnsemblID]%>%unique()%>%length()

matrix <- data.frame(
  as.numeric(numb_da_de_genes),
  as.numeric(numb_nonda_de_genes), 
  as.numeric(numb_da_nonde_genes),
  as.numeric(numb_nonda_nonde_genes)
)%>%unlist()%>%matrix(nrow=2,byrow=T)

enrichment_dade_overlap <-fisher.test(matrix) # 3.56e-05
## focus on genes that are DE and DA 
## and look at the direction of their respective logFC
da_de_genes <- significant_da_genes[significant_da_genes %in% significant_de_genes]

logfc_da_de_genes <- copy(one_peak_per_gene)[EnsemblID %in% da_de_genes][DA =='da']%>%unique()

count_genes_per_quadrant = function(x){
  df = copy(x)[
    ,quadrant:=ifelse((logFC>0 & i.logFC<0),'1',
      ifelse((logFC>0 & i.logFC>0),'2',
        ifelse((logFC<0 & i.logFC<0),'3','4'
        )))
        ][
          ,maxDAfc:=max(logFC)
          ][
            ,minDAfc:=min(logFC)
            ][
              ,maxDEfc:=max(i.logFC)
              ][
                ,minDEfc:=min(i.logFC)
                ][
                  ,c("EnsemblID",'quadrant','maxDAfc','minDAfc','maxDEfc','minDEfc')
                  ]%>%unique()
  df = df[
    ,numbgenes:=.N,by=.(quadrant)
    ][
      ,EnsemblID:=NULL
      ]%>%unique()%>%setorderv('quadrant',1)

  df = df[
    ,y:=ifelse((quadrant=='1'|quadrant=='2'),maxDAfc+0.5,minDAfc-0.5)
    ][
      ,x:=ifelse((quadrant=='1'|quadrant=='3'),minDEfc+0.5,maxDEfc-0.5)
      ]
  return(df)
}

logfc_plot = function(x){
  df=copy(x)
  numb_genes = count_genes_per_quadrant(x)
  p = ggplot(df, aes(x=i.logFC, y=logFC)) +
    geom_point()+
    geom_text(data = numb_genes, aes(label = numb_genes$numbgenes,y = numb_genes$y, x = numb_genes$x))+
    geom_vline(xintercept=0,linetype='dashed')+
    geom_hline(yintercept=0,linetype='dashed')+
    ylab('DA logFC')+xlab('DE logFC')+
    # facet_wrap(peak_species~.,ncol=3)+
    theme()
  
  return(p)
}

pdf(paste(outplot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 8,height = 5)
logfc_plot(logfc_da_de_genes)
dev.off()
