## script used to perform GO enrichment of common and species specific DA peaks using GREAT
library(dplyr)
library(data.table)
library(magrittr)
library(rGREAT)
library(GenomicRanges)
library(openxlsx)
library(ggthemes)
library(ggplot2)
library(ggpubr)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'GO_enrich')
outfile_dir <- create_dir(outdir,'files/GO_enrich')

# table_dir <- './post_processing_analyses/output/tables/GO_enrich/'
# target_genes_dir <- './post_processing_analyses/output/GO/target_genes/'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))
setkeyv(da_results,range_keys)
da_results <- da_results%>%setorderv('peak_species',1)

all_peaks <- split(da_results,by='DA')

## GREAT enrichments
background <- copy(da_results)

get_enrichment = function(gr,bg){
    enrich = submitGreatJob(
    gr=gr,
    bg=bg,
    species = genome,
    rule= "twoClosest", 
    adv_twoDistance = 1000,  
    includeCuratedRegDoms = T,
    request_interval=10)
    
    go_enrich = getEnrichmentTables(enrich,ontology='GO Biological Process')%>%rbindlist()
    go_enrich = go_enrich[Hyper_Adjp_BH <= 0.05]
    
    # get list of putatve target genes 
    target_genes = plotRegionGeneAssociationGraphs(enrich,type=1,request_interval = 10)%>%as.data.table() %>% na.omit() %>% unique()
    great_enrich = list(go_enrich,target_genes)
    return(great_enrich)
}

go_enrichment <- lapply(all_peaks,function(x)get_enrichment(x,background))

go_tables <- list(go_enrichment[[1]][[1]],go_enrichment[[2]][[1]])
names(go_tables) = names(go_enrichment)

write.xlsx(go_tables,paste(outfile_dir,'GO_enrich_terms_da_non_da.xlsx',sep=''),append=T,overwrite=T)

## plot go bp first 30 terms for simplicity of visualization
go_bp <- Map(mutate,go_tables,DA = as.factor(names(go_tables)))
go_bp <- lapply(go_bp,function(x) x=x[1:15,c(1,2,13,14)])%>%rbindlist()%>%na.omit()
go_bp <- go_bp[
  ,log10_adj_p := -log10(Hyper_Adjp_BH)
]%>%setorderv('log10_adj_p',1)%>%setorderv('DA',1)

go_enrichment_plot = function(df){
  p <- ggplot(df, aes(x=factor(df$name,levels=df$name), y=log10_adj_p,fill=DA)) +
  geom_bar(stat = 'identity',position = 'dodge',col='black')+
  scale_fill_manual(values=da_palette)+
  xlab(" ") +ylab("\n -Log10 (P adj.) \n ") +
  theme(
    legend.position='bottom',
    legend.background=element_rect(),
    axis.text.x=element_text(angle=0, hjust=1.10),
    axis.text.y=element_text(angle=0, vjust=0.8),
    axis.title=element_text(),
    axis.line = element_line(color = "black",size = 0, linetype = "solid"),
    panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    title=element_text(),
    strip.background = element_blank()
    ) +
  coord_flip()
  return(p)
  
}

pdf(paste(outplot_dir,'go_barplot.pdf',sep=''),width=10,height = 15)
go_enrichment_plot(go_bp)
dev.off()

## check ontology index of significant GO terms 
library(ontologyIndex)
data(go)

significant_gos <- copy(go_tables)%>%lapply(function(x)x[,'ID']%>%unique()%>%split(by='ID'))


ancestor_terms <- copy(significant_gos)%>%lapply(
  function(x) x<-lapply(x,function(y)
  y<-length(get_ancestors(go,as.character(y))))
)

ancestor_terms <- lapply(ancestor_terms,function(x){
  goids <- names(x)
  table <- lapply(x,function(y)y=data.table(ancestor_terms=y))%>%rbindlist()
  table <- table[,goids:=goids]
  return(table)
})

ancestor_terms <- Map(mutate,ancestor_terms,DA=names(ancestor_terms))%>%rbindlist()

pdf(paste0(outplot_dir,'ancestor_terms.pdf',sep=''),width = 10, height = 8)
ggplot(ancestor_terms,aes(x=DA,y=ancestor_terms,fill=DA))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
scale_fill_manual(values=da_palette)+
stat_compare_means(
  method = "wilcox.test",
  ref.group = 'non_da',
  size=5
)+
ylab('number of ancestor terms per GO')+
theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

##-------------------------
## check the target genes
##-------------------------
target_genes = list(go_enrichment[[1]][[2]],go_enrichment[[2]][[2]])
target_genes = purrr::map2(target_genes,all_peaks,function(x,y)x[y,on=range_keys,nomatch=0])
names(target_genes)=c('da','non_da')

## add ensembl ID info
library(biomaRt)

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl") 

# genes_assoc_w_signf_go <- lapply(
#   go_tables,function(x){ 
#     go_id_genes <-copy(x)[,'ID']%>%split(by='ID')%>%
#     lapply(
#       function(y)y<-getBM(
#         attributes = c('hgnc_symbol'), 
#         filters = 'go', 
#         values = x$ID, 
#         mart = ensembl
#         )
#         )%>%rbindlist()
#   }
# )

## get loeuf score genes targeted by da, non-da and da + non-da peaks 
gnomad_loeuf <- fread('../../../../punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T,select = c('gene','gene_id','oe_lof_upper'))
colnames(gnomad_loeuf)[1] = 'EnsemblID'

target_genes_w_loeuf <- copy(target_genes)%>%lapply(function(x){
  x<-x[gnomad_loeuf,on='gene',nomatch=0]%>%unique()%>%na.omit()
  }
)%>%rbindlist()

## divide genes by peak 
target_genes_w_loeuf <- target_genes_w_loeuf[
  ,c('gene','oe_lof_upper','DA')
]%>%unique()

target_genes_w_loeuf <-target_genes_w_loeuf[
  ,peaktype_assoc:=.N,by=.(gene)
  ][
    ,regulation:=ifelse(peaktype_assoc>1,'nonda_da',DA)
    ][
      ,DA:=NULL
]

all_genes_loeuf <- copy(gnomad_loeuf)[,gene_id:=NULL][,regulation:='all_genes']%>%na.omit()


comparisons = list(
  # c('non_da','da'),
  c('da','all_genes'),
  c('non_da','all_genes'),
  c('nonda_da','all_genes')
  # c('nonda_da','non_da'),
  # c('nonda_da','da')
)

pdf(paste0(outplot_dir,'loeuf_score_target_genes.pdf',sep=''),width = 10, height = 8)
df <- rbind(all_genes_loeuf,target_genes_w_loeuf%>%dplyr::select(c(all_of(colnames(all_genes_loeuf)))))
ggplot(df,aes(x=regulation,y=oe_lof_upper,fill=regulation))+
geom_violin(trim=T,scale = "width")+
geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
geom_hline(yintercept=median(all_genes_loeuf$oe_lof_upper),linetype='dashed',size=0.5)+
stat_compare_means(
  method = "wilcox.test",
  comparisons = comparisons,
  size=5
  )+
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()



## get number and % peaks with target genes
numb_peaks_w_gene <- copy(target_genes)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
numb_peaks <- copy(all_peaks)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())

purrr::map2(numb_peaks_w_gene,numb_peaks,`/`)

numb_genes <- copy(target_genes)%>%lapply(function(x)x=x[,'gene']%>%unique()%>%nrow())

filenames <- paste0(outfile_dir,paste(names(target_genes),'target_genes.txt',sep='_'),sep='')
mapply(write.table,target_genes, file = filenames,col.names = T, row.names = F, sep = "\t", quote = F)

## plot venn diagramm of shared genes between peaks
library(VennDiagram)
 
genes <- copy(target_genes)%>%lapply(function(x)x=x[,gene]%>%unique())
venn.diagram(
    x = genes,
    category.names = c('da','non_da'),
    filename = paste(outplot_dir,'target_genes_venn.png',sep=''),
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
    cat.pos = c(-27, 27),
    cat.dist = c(0.055, 0.055),
    cat.fontfamily = "sans",
    cat.col = da_palette
)

##  plot distance between peak and target genes
target_genes <- lapply(target_genes,function(x)x=x[
    ,log10_abs_dist:= log10(abs(distTSS)+1)
    ]
)%>%rbindlist()

# target_genes <- Map(mutate,target_genes,peak_type=names(target_genes))%>%rbindlist()

# pdf(paste(outplot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
# ggplot(target_genes,aes(x=DA,y=log10_abs_dist,fill=DA))+
# geom_violin(trim=T,scale = "width")+
# geom_boxplot(width=.1, position =  position_dodge(width = 0.4),outlier.size=0.2,fill='white',notch=T)+
# xlab('')+ylab('Log10 bp dist from TSS')+
# scale_fill_manual(values=da_palette)+
# stat_compare_means(
#   method = "wilcox.test",
#   label.y = (max(target_genes$log10_abs_dist)+0.5),
#   size=5
#   )+
# theme(
#   legend.key = element_rect(fill = "white", colour = "black"),
#   axis.line = element_blank()
# )
# dev.off()

pdf(paste(outplot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
df <- copy(target_genes)[,log10_dist:=ifelse(distTSS<0,-log10_abs_dist,log10_abs_dist)]
ggplot(df, aes(x=log10_dist,fill=DA)) +
    geom_density(alpha=0.5)+scale_fill_manual(values=da_palette)
dev.off()


