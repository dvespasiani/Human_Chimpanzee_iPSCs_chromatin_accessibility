## script used to look for overlap between DA peaks and DE genes

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(liftOver)
library(rtracklayer)
library(biomaRt)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

peakDir = './post_processing_analyses/output/DA/peaks/'
target_genes_dir = './post_processing_analyses/output/GO/target_genes/'
tmp_files_dir = './post_processing_analyses/output/temp_files/'

plot_dir = './post_processing_analyses/output/plots/rna_seq/'

##-----------------
## read DE genes 
##-----------------
## NB: because the DE was tested chimp vs human whereas DA was human vs chimp simply revert the sign of the DE logFC 
de_genes = fread("./rna_seq/de_output/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B'))
de_genes = de_genes[
    ,logFC:=-logFC
    ][
        ,de_significant:=ifelse(adj.P.Val<= 0.05,'significant','non_significant')
]

## read peaks
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))
common_da_peaks = read_da_peaks('common_da_results.txt',c('significant','non_significant'))

common_da_peaks = common_da_peaks[,peakID:=paste(paste('C_',chimp_peakID,sep=''),paste('H_',human_peakID,sep=''),sep='.')]
# common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))[,species:='common']

all_peaks = list(chimp_peaks,common_da_peaks,human_peaks)%>%lapply(function(x)x=x[,c('logFC','logCPM','FDR','significant','peakID')]%>%unique())

## read list peaks w genes same TAD + regulation info
peaks_w_genes_same_tad=fread(paste(tmp_files_dir,'list_peaks_w_genes_same_tad.txt',sep=''),header=T,sep='\t')[,species:=peak_species]%>%setorderv('species',1)%>%split(by='species')
lapply(peaks_w_genes_same_tad,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

## add logFC into to these peaks and bind the lists
peaks_w_genes_same_tad = copy(peaks_w_genes_same_tad)%>%
purrr::map2(
  all_peaks,function(x,y)
  x[
    y,on=c('significant','peakID'),nomatch=0
    ][
        ,c('tadID','numbpeaks_per_gene_same_tad','peak_species'):=NULL]%>%unique()
)%>%rbindlist()

## retain only genes with ensembl ids from v 86
ensembl_86_ids = fread("./rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

peaks_w_genes_same_tad_ensembl = copy(peaks_w_genes_same_tad)[
    ensembl_86_ids,on='gene',nomatch=0
]
## get the prop of genes retained after this step
length(unique(peaks_w_genes_same_tad_ensembl$gene))/length(unique(peaks_w_genes_same_tad$gene))

## now add DE info to the remaining genes
peaks_w_genes_same_tad_ensembl = peaks_w_genes_same_tad_ensembl[de_genes,on='EnsemblID',nomatch=0]

##------------------------------------------------------
## look at avg expression genes by peaks (DA vs non-DA)
##------------------------------------------------------

comparisons = list(
  c('significant','non_significant')
)

pdf(paste0(plot_dir,'/logFC_genes_by_peak_significance.pdf',sep=''),width = 7, height = 8)
p <- ggboxplot(
  peaks_w_genes_same_tad_ensembl, x = "significant", y = "AveExpr",
  fill = "significant",notch = T,short.panel.labs = FALSE
)+xlab(' ')+ylab('Avg gene expression')+ scale_fill_discrete(name="ATAC peaks")
p + stat_compare_means(
  comparisons = comparisons,
  show.legend=FALSE,
  method='wilcox.test', 
  p.adjust.method = "fdr",
  color = "black", size = 7
  ) +
theme(
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)

dev.off()

## look at ovelap between DE genes and DA peaks
library(VennDiagram)
 
significant_da_genes = copy(da_peaks_w_genes_same_tad_ensembl)[significant=='significant'][,EnsemblID]%>%unique()
significant_de_genes = copy(de_genes)[de_significant=='significant']$EnsemblID

da_de_overlap= list(significant_da_genes,significant_de_genes)
venn.diagram(
    x = da_de_overlap,
    category.names = c("DA","DE"),
    filename = paste(plot_dir,'da_de_overlap.png',sep=''),
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
da_de_genes = significant_da_genes[significant_da_genes %in% significant_de_genes]

logfc_da_de_genes = copy(da_peaks_w_genes_same_tad_ensembl)[ EnsemblID %in% da_de_genes][significant=='significant']
logfc_da_de_genes = logfc_da_de_genes[, .SD[which.min(abs(distTSS))], by=.(EnsemblID)]
# ## look direction logFC atac and genes 
# ## remove the non significant ones (both DE and DA only)
# accessibility_expression = copy(da_peaks_w_de_genes_same_tad)%>%lapply(
#     function(x)x=x[
#             ,c(..range_keys,'EnsemblID','logFC',"i.logFC",'significant','de_significant',"peak_species",'distTSS','peakID')
#             ][
#               , .SD[which.min(abs(distTSS))], by=.(EnsemblID)
#                   ]
# )%>%rbindlist()
# # [,logFC:=ifelse(peak_species=='chimp',-logFC,logFC) ## for chimp revert the sign of DA and DE
# # ][i.logFC:=ifelse(peak_species=='chimp',-i.logFC,i.logFC) ## for chimp revert the sign of DA and DE
# # ]
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

pdf(paste(plot_dir,'logFC_DA_DE_direction.pdf',sep=''),width= 8,height = 5)
logfc_plot(logfc_da_de_genes[regulation=='specific'])
dev.off()

## now look what happens if u filter for peaks nearby genes (5kb)
pdf(paste(plot_dir,'/proximal_peaks_logFC_DA_DE_direction.pdf',sep=''),width= 8,height = 5)
logfc_plot(logfc_da_de_genes[abs(distTSS)<=5000])
dev.off()


cor.test(logfc_da_de_genes$logFC,logfc_da_de_genes$i.logFC)






# # common_regions = common_regions[
# #     ,human_peakID:=paste('peak_',gsub('.*_','',peakID),sep='')
# # ]

# # ## get logFC for the common peaks
# # common_da_results = read_da_peaks('common_da_results.txt',c('significant','non_significant'))

# # common_regions = common_regions[
# #     common_da_results,on=c('human_peakID','significant'),nomatch=0
# #     ][
# #         ,c(paste('human',range_keys,sep='_')):=NULL
# # ]


# # all_peaks = list(chimp_peaks,common_regions,human_peaks)

# # lapply(all_peaks,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

# # target_genes = list.files(target_genes_dir,full.names=T,recursive=F)%>%
# # lapply(
# #     function(x)fread(x,sep='\t',header=T)
# # )
# # names(target_genes) = gsub('\\..*','',list.files(target_genes_dir,full.names=F,recursive=F))

# ## liftover back to pantro5 the chimp peaks with target genes
# chain_path <- './data/LiftOver_chains/'

# chimp_peaks_w_targets = copy(target_genes[[1]])
# chimp_peaks_w_targets = chimp_peaks_w_targets[,peakID:=paste('peak_',1:nrow(chimp_peaks_w_targets),sep='')]

# chimp_peaks_w_target_pantro5 = convert_coord(chimp_peaks_w_targets,'hg38ToPanTro5.over.chain')
# chimp_peaks_w_target_pantro5 = chimp_peaks_w_target_pantro5[,peakID:=NULL]

# new_target_genes = list(chimp_peaks_w_target_pantro5,target_genes[[2]],target_genes[[3]])

# lapply(new_target_genes,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())

# peaks_w_target_genes = purrr::map2(all_peaks,new_target_genes,function(x,y)x[y,on=c(range_keys),nomatch=0])
# names(peaks_w_target_genes) = species_names

# lapply(peaks_w_target_genes,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())


# # ## now get the ensembl id for those genes
# # pantro_mart <- useMart(biomart="ensembl", dataset="ptroglodytes_gene_ensembl")
# # human_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")


# # orth_genes = getBM(
# #     attributes = c('ptroglodytes_homolog_ensembl_gene','ensembl_gene_id'), ## get all orth genes
# #     mart = human_mart
# # )%>%as.data.table()

# # orth_genes = orth_genes[,orth:=ifelse(ptroglodytes_homolog_ensembl_gene=='','no_ort','orth')]

# # human_genes_with_orth = orth_genes[orth=='orth']$ensembl_gene_id
# # chimp_genes_with_orth = orth_genes[orth=='orth']$ptroglodytes_homolog_ensembl_gene


# ##-----------------------------------
# ## as long as ensembl server is off 
# ##-----------------------------------
# ## read all ensembl ids and gene names from most recent ensembl version (104)
# ## this is done to match GREAT assigned gene names
# # ensembl_104_ids = fread("./rna_seq/de_output/all_ensembl_gene_ids.txt",sep='\t',header=T)

# # ## assign these ensembl ids to the gene names
# # ## merge ensemblID info
# # peaks_w_target_genes = lapply(peaks_w_target_genes,function(x)x[
# #     ensembl_104_ids,on='gene',nomatch=0
# #     ]
# # )
# # lapply(peaks_w_target_genes,function(x)x[,c(..range_keys)]%>%unique()%>%nrow()) ## in this way 99% of the genes get assigned


# ## now read ensembl gene ids from archive 87
# ensembl_86_ids = fread("./rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

# peaks_w_target_genes = lapply(peaks_w_target_genes,function(x)x[
#     ensembl_86_ids,on='gene',nomatch=0
#     ]
# )
# lapply(peaks_w_target_genes,function(x)x[,c(..range_keys)]%>%unique()%>%nrow())


# # ensembl_id_pantro_orthologous = fread("./rna_seq/de_output/ensembl_id_pantro_orthologous.txt",sep='\t',header=T) 
# # ensembl_id_pantro_orthologous = ensembl_id_pantro_orthologous[
# #     ,seqnames:=paste('chr',seqnames,sep='')
# # ]

# ## merge ensemblID info
# peaks_w_target_genes = lapply(peaks_w_target_genes,function(x)x[
#     ensembl_id_hugo_symbol,on='gene',nomatch=0
#     ]
# )

# ## first get the % of genes targeted by DA peaks
# ## for which we have expression data
# peaks_w_genes_w_expr_data = lapply(peaks_w_target_genes,function(x)merge(x,de_genes,by='EnsemblID'))

# lapply(peaks_w_genes_w_expr_data,function(x)length(unique(x$EnsemblID)))
# length(de_genes$EnsemblID)

# prop_target_genes_w_expr_data = copy(peaks_w_genes_w_expr_data)%>%
# lapply(
#     function(x)
#     x=length(unique(x$EnsemblID))/length(de_genes$EnsemblID)
# )
# prop_target_genes_wo_expr_data = copy(prop_target_genes_w_expr_data)%>%lapply(function(x)1-x)

# prop_target_genes = purrr::map2(
#     prop_target_genes_w_expr_data,prop_target_genes_wo_expr_data,
#     function(x,y)
#     z=rbind(x,y)%>%as.data.table()%>%
#     mutate('with_expression'=c('yes','no'))%>%
#     setnames(old=1,new='prop')
# )

# prop_target_genes = Map(mutate,prop_target_genes,species=names(prop_target_genes))%>%rbindlist()

# pdf(paste(plot_dir,'prop_target_gene_expr_data.pdf',sep=''),width= 7,height = 7)
# ggplot(prop_target_genes, aes(x=species, y=prop,fill=with_expression)) +
#     geom_bar(stat='identity')+
#     ylab('proportion target genes \n with expression data')+xlab('species')+
#     theme(
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
# dev.off()

# ## now check the proportion of target genes that are DE
# prop_target_de =  copy(peaks_w_genes_w_expr_data)%>%lapply(
#     function(x)x=x[,c('EnsemblID','species','de_significant')]%>%unique()
# )%>%rbindlist()
# prop_target_de = prop_target_de[
#     ,numb_genes:=.N,by=.(species)
#     ][
#         ,numb_de_genes:=.N,by=.(species,de_significant)
#         ][
#             ,prop_de_genes:=numb_de_genes/numb_genes
#             ][
#                 ,c('species','prop_de_genes','de_significant')
# ]%>%unique()

# pdf(paste(plot_dir,'prop_target_gene_de.pdf',sep=''),width= 7,height = 7)
# ggplot(prop_target_de, aes(x=species, y=prop_de_genes,fill=de_significant)) +
#     geom_bar(stat='identity')+
#     ylab('proportion target genes')+xlab('species')+
#     theme(
#         axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
#         )
# dev.off()

