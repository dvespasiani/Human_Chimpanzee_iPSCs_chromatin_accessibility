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
genome <- 'hg38'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC'))
setkeyv(da_results,range_keys)

da_peaks <- copy(da_results)[DA=='da']
nonda_peaks <- copy(da_results)[DA=='non_da']

all_peaks = list(da_peaks,nonda_peaks)%>%lapply(function(x)x<-x[,c(..range_keys,'peakID')]%>%unique())
names(all_peaks) = c('da','non_da')

## GREAT enrichments
background <- copy(all_peaks)%>%rbindlist()

get_enrichment = function(test,background){
    enrich = submitGreatJob(
    gr=test,
    bg=background,
    species = "hg38",
    rule= "oneClosest", 
    adv_oneDistance = 1000,  
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
names(go_tables) = c('da','non_da')

write.xlsx(go_tables,paste(outfile_dir,'GO_enrich_terms_da_non_da.xlsx',sep=''),append=T,overwrite=T)

## plot go bp first 30 terms for simplicity of visualization
go_bp <- Map(mutate,go_tables,peak_set = as.factor(names(go_tables)))
go_bp <- lapply(go_bp,function(x) x=x[1:15,c(1,2,13,14)])%>%rbindlist()
go_bp <- go_bp[
  ,log10_adj_p := -log10(Hyper_Adjp_BH)
]%>%setorderv('log10_adj_p',1)%>%setorderv('peak_set',1)

go_enrichment_plot = function(df){
  p <- ggplot(df, aes(x=factor(df$name,levels=df$name), y=log10_adj_p,fill=peak_set)) +
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

pdf(paste(outplot_dir,'da_nonda_go.pdf',sep=''),width=10,height = 15)
go_enrichment_plot(go_bp)
dev.off()

## now check the target genes
target_genes = list(go_enrichment[[1]][[2]],go_enrichment[[2]][[2]])
target_genes = purrr::map2(target_genes,all_peaks,function(x,y)x[y,on=range_keys,nomatch=0])
names(target_genes)=c('da','non_da')

## get number and % peaks with target genes
numb_peaks_w_gene <- copy(target_genes)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
numb_peaks <- copy(all_peaks)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())

# purrr::map2(numb_peaks_w_gene,numb_peaks,`/`)

numb_genes <- copy(target_genes)%>%lapply(function(x)x=x[,'gene']%>%unique()%>%nrow())

filenames <- paste0(outfile_dir,paste(names(target_genes),'target_genes.txt',sep='_'),sep='')
mapply(write.table,target_genes, file = filenames,col.names = T, row.names = F, sep = "\t", quote = F)

## plot venn diagramm of shared genes between peaks
library(VennDiagram)
 
genes <- copy(target_genes)%>%lapply(function(x)x=x[,gene]%>%unique())
venn.diagram(
    x = genes,
    category.names = c("da", "non_da"),
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
)

target_genes <- Map(mutate,target_genes,peak_type=names(target_genes))%>%rbindlist()

# pdf(paste(plot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
# ggplot(target_genes,aes(x=peak_type,y=log10_abs_dist,fill=peak_type))+
# geom_boxplot(notch=T)+
# xlab('')+ylab('Log10 bp dist from TSS')+
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
ggplot(target_genes, aes(x=log10_abs_dist,fill=peak_type)) +
    geom_density(alpha=0.5)+scale_fill_manual(values=da_palette)
dev.off()


