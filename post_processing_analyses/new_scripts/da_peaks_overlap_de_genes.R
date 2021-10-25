## script used to look for overlap between DA peaks and DE genes

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)

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
]%>%setorderv('peak_species',1)%>%split(by='peak_species')

## read target genes
target_genes <- list.files(target_genes_dir,full.names= T,recursive=F)%>%lapply(function(x)fread(x,sep='\t',header=T))

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
        ,de_significant:=ifelse(adj.P.Val<= 0.05,'de','non_de')
]

# all_peaks = list(chimp_peaks,common_da_peaks,human_peaks)%>%lapply(function(x)x=x[,c('logFC','logCPM','FDR','significant','peakID')]%>%unique())

# ## read list peaks w genes same TAD + regulation info
# peaks_w_genes_same_tad <- fread(paste(tmp_files_dir,'list_peaks_w_genes_same_tad.txt',sep=''),header=T,sep='\t')[,species:=peak_species]%>%setorderv('species',1)%>%split(by='species')
# lapply(peaks_w_genes_same_tad,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

# ## add logFC into to these peaks and bind the lists
# peaks_w_genes_same_tad = copy(peaks_w_genes_same_tad)%>%
# purrr::map2(
#   all_peaks,function(x,y)
#   x[
#     y,on=c('significant','peakID'),nomatch=0
#     ][
#         ,c('tadID','numbpeaks_per_gene_same_tad','peak_species'):=NULL]%>%unique()
# )%>%rbindlist()

## retain only genes with ensembl ids from v 86
ensembl_86_ids <- fread("../rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

peaks_w_genes_ensembl <- copy(peaks_w_genes)[ensembl_86_ids,on='gene',nomatch=0]

## get the prop of genes retained after this step
# length(unique(peaks_w_genes_same_tad_ensembl$gene))/length(unique(peaks_w_genes_same_tad$gene))

## now add DE info to the remaining genes
peaks_w_de_genes_ensembl <- copy(peaks_w_genes_ensembl)[de_genes,on='EnsemblID',nomatch=0]

all_genes_expr <- rbind(
  peaks_w_de_genes_ensembl[,c('AveExpr','peak_species')],
  de_genes[,c('AveExpr')][,peak_species:='all_genes']
)
##------------------------------------------------------
## look at avg expression genes by peaks (DA vs non-DA)
##------------------------------------------------------
comparisons = list(
  c('common','chimp'),
  c('common','human'),
  c('common','all_genes')

)

plot_de_genes <- function(df,xaxis,yaxis,fill){
  # labely <- copy(df)%>%dplyr::pull(yaxis)
  # labely <- max(labely)
  p <- ggboxplot(
    df, x = xaxis, y = yaxis,fill = xaxis,
    notch = T,order=c('all_genes','common','chimp','human')
    )+ 
    xlab(' ')+ylab(yaxis)+
    stat_compare_means(
     method = "wilcox.test",
    #  label.y = labely,
     size=5,
     ref.group='common'
  )+
  theme(
  legend.position = "none",
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )
  return(p)
}

pdf(paste0(outplot_dir,'logFC_genes_by_peak_significance.pdf',sep=''),width = 10, height = 8)
plot_de_genes(all_genes_expr,'peak_species','AveExpr')
dev.off()


## look at ovelap between DE genes and DA peaks
library(VennDiagram)
 
significant_da_genes <- copy(peaks_w_de_genes_ensembl)[DA=='da'][,EnsemblID]%>%unique()
significant_de_genes <- copy(de_genes)[de_significant=='de']$EnsemblID

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

## focus on genes that are DE and DA 
## and look at the direction of their respective logFC
da_de_genes <- significant_da_genes[significant_da_genes %in% significant_de_genes]

logfc_da_de_genes <- copy(peaks_w_de_genes_ensembl)[ EnsemblID %in% da_de_genes][
  DA =='da'
  ][
    , .SD[which.min(abs(distTSS))], by=.(EnsemblID)
]


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
    ylab('DA logFC')+xlab('DE logFC')+
    # facet_wrap(peak_species~.,ncol=3)+
    theme()
  
  return(p)
}

pdf(paste(outplot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 8,height = 5)
logfc_plot(logfc_da_de_genes)
dev.off()

# ## now look what happens if u filter for peaks nearby genes (5kb)
# pdf(paste(plot_dir,'/proximal_peaks_logFC_DA_DE_direction.pdf',sep=''),width= 8,height = 5)
# logfc_plot(logfc_da_de_genes[abs(distTSS)<=5000])
# dev.off()

## loeuf score for genes
gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene_id','oe_lof_upper'))
colnames(gnomad_loeuf)[1] = 'EnsemblID'

## add loeuf score to genes
peaks_w_de_genes_ensembl_loeuf <- copy(peaks_w_de_genes_ensembl)[gnomad_loeuf,on='EnsemblID',nomatch=0]

all_genes_loeuf <- copy(gnomad_loeuf)[,peak_species:='all_genes']
peak_genes_loeuf <- copy(peaks_w_de_genes_ensembl_loeuf)%>%dplyr::select(c(all_of(colnames(all_genes_loeuf))))
combined_loeuf <- rbind(all_genes_loeuf,peak_genes_loeuf)

pdf(paste0(outplot_dir,'loeuf_score_genes_by_peak.pdf',sep=''),width = 10, height = 8)
plot_de_genes(combined_loeuf,'peak_species','oe_lof_upper')
dev.off()
