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
library(liftOver)
library(rtracklayer)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))


chrom_state_dir <- './data/iPSC_chrom_states_hg38'
plot_dir <- './post_processing_analyses/output/plots/GO_enrich/'
da_dir <- './post_processing_analyses/output/DA/'
table_dir <- './post_processing_analyses/output/tables/GO_enrich/'
target_genes_dir <- './post_processing_analyses/output/GO/target_genes/'

## get DA peaks
da_file <- paste(da_dir,genome,'/',genome,'_da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC'))
setkeyv(da_results,range_keys)


species_colors <- c("#F8766D","#7CAE00","#00BFC4")

## using peaks after DA analysis (both da and non DA)
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))

chimp_peaks_hg38 = convert_coord(chimp_peaks,'panTro5ToHg38.over.chain')


all_peaks = list(chimp_peaks_hg38,common_regions,human_peaks)%>%lapply(function(x)x=x[,c(..range_keys,'peakID')]%>%unique())
names(all_peaks)=species_names

## GREAT enrichments
background_common = copy(all_peaks)%>%rbindlist()
background_chimp_specific = copy(all_peaks[c(1:2)])%>%rbindlist()
background_human_specific = copy(all_peaks[c(2,3)])%>%rbindlist()

get_enrichment = function(test,background){
    enrich = submitGreatJob(
    gr=test,
    bg=background,
    species = "hg38",
    rule= "oneClosest", 
    adv_oneDistance=1000,  
    includeCuratedRegDoms = T,
    request_interval=10)
    
    go_enrich = getEnrichmentTables(enrich,ontology='GO Biological Process')%>%rbindlist()
    go_enrich = go_enrich[Hyper_Adjp_BH <= 0.05]
    
    # get list of putatve target genes 
    target_genes = plotRegionGeneAssociationGraphs(enrich,type=1,request_interval = 10)%>%as.data.table() %>% na.omit() %>% unique()
    great_enrich = list(go_enrich,target_genes)
    return(great_enrich)
}

chimp_enrich = get_enrichment(all_peaks[[1]],background_chimp_specific)
common_enrich = get_enrichment(all_peaks[[2]],background_common)
human_enrich = get_enrichment(all_peaks[[3]],background_human_specific)

go_tables = list(chimp_enrich[[1]],common_enrich[[1]],human_enrich[[1]])
names(go_tables)=species_names

write.xlsx(go_tables,paste(table_dir,'all_da_peaks_signif_GO_enrich_terms.xlsx',sep=''),append=T)

## plot go bp first 30 terms for simplicity of visualization
go_bp=Map(mutate,go_tables,'peak_set'=names(go_tables))
go_bp=lapply(go_bp,function(x) x=x[1:10,c(1,2,13,14)])

go_enrichment_plot=function(x,y){
  ggplot(x, aes(x=reorder(name,-log10(Hyper_Adjp_BH)), y=-log10(Hyper_Adjp_BH),fill=peak_set)) +
    geom_bar(stat = 'identity',position = 'dodge')+
    xlab(" ") +
    ylab("\n -Log10 (adj. P) \n ") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(x$Hyper_Adjp_BH)), by = 2), 1)) +
    scale_fill_manual(name= " ",values=y)+
   theme(
      legend.position='none',
      legend.background=element_rect(),
      axis.text.x=element_text(angle=0, hjust=1.10),
      axis.text.y=element_text(angle=0, vjust=0.8),
      axis.title=element_text(),
      axis.line = element_line(color = "black",size = 0, linetype = "solid"),
     panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      title=element_text()) +
    coord_flip()
  
}

pdf(paste(plot_dir,'chimp_go_plot.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(go_bp$chimp,species_colors[[1]])
dev.off()

pdf(paste(plot_dir,'common_go_plot.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(go_bp$common,species_colors[[2]])
dev.off()

pdf(paste(plot_dir,'human_go_plot.pdf',sep=''),width=10,height = 7)
go_enrichment_plot(go_bp$human[1:10,],species_colors[[3]])
dev.off()

## now check the target genes
target_genes = list(chimp_enrich[[2]],common_enrich[[2]],human_enrich[[2]])
target_genes = purrr::map2(target_genes,all_peaks,function(x,y)x[y,on=range_keys,nomatch=0])
names(target_genes)=species_names

## get number and % peaks with target genes
numb_peaks_w_gene=copy(target_genes)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
numb_peaks=copy(all_peaks)%>%lapply(function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())

purrr::map2(numb_peaks_w_gene,numb_peaks,`/`)

numb_genes=copy(target_genes)%>%lapply(function(x)x=x[,'gene']%>%unique()%>%nrow())


filenames=paste0(target_genes_dir,paste(names(target_genes),'target_genes.txt',sep='_'),sep='')
mapply(write.table,target_genes, file = filenames,col.names = T, row.names = F, sep = "\t", quote = F)

## plot venn diagramm of shared genes between peaks
library(VennDiagram)
 
genes = copy(target_genes)%>%lapply(function(x)x=x[,gene]%>%unique())
venn.diagram(
    x = genes,
    category.names = c("chimp",'common', "human"),
    filename = paste(plot_dir,'target_genes_venn.png',sep=''),
    output = TRUE ,
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 400,
    lwd = 1,
    col=species_colors,
    fill = c(alpha(species_colors[[1]],0.3), alpha(species_colors[[2]],0.3), alpha(species_colors[[3]],0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = species_colors
)

##  plot distance between peak and target genes
target_genes = lapply(target_genes,function(x)x=x[
    ,log10_abs_dist:= log10(abs(distTSS)+1)
    ]
)

target_genes = Map(mutate,target_genes,file=names(target_genes))%>%rbindlist()

pdf(paste(plot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
ggplot(target_genes,aes(x=file,y=log10_abs_dist,fill=file))+
geom_boxplot(notch=T)+
xlab('')+ylab('Log10 bp dist from TSS')+
stat_compare_means(
  method = "wilcox.test",
  label.y = (max(target_genes$log10_abs_dist)+0.5),
  size=5,
  ref.group='common'
  )+
theme(
  legend.key = element_rect(fill = "white", colour = "black"),
  axis.line = element_blank()
)
dev.off()

pdf(paste(plot_dir,'dist_peak_target_genes.pdf',sep=''),width=10,height = 7)
ggplot(target_genes, aes(x=log10_abs_dist,fill=file)) +
    geom_density(alpha=0.5)
dev.off()


