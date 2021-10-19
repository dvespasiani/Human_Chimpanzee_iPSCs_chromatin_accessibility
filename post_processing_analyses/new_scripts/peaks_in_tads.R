## script used to look for overlaps between peaks and TADs (eres et al 2019)
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes)
library(ggplot2)
library(ggpubr)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

go_dir <- './output/GO/'
tads_dir <-  './output/TADs/'
peakDir <- './output/DA/peaks/'
tmp_files_dir <- './output/temp_files/'

outplot_dir <- create_dir(plot_dir,'TADs')
genome <- 'hg38'

## get DA peaks
da_file <- paste(da_dir,genome,'/','da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','logFC'))
setkeyv(da_results,range_keys)


##------------
## read TADs
##------------
hg38_tads <- read_tads('hg38_tads.bed')[,species:=ifelse(label==1,'Common',species)][,label:=NULL]


all_tads = list(chimp_tads,human_tads)
all_tads = lapply(all_tads,function(x) x=x[,tadID:=paste('tad_',seq(1:nrow(x)),sep='')]) ## use this column only to split dataframes
lapply(all_tads,function(x)setkeyv(x,range_keys))
names(all_tads) = species_names[-2]

##---------------
## read peaks
##---------------
human_peaks = read_da_peaks('human_da_results.txt',c('significant','non_significant'))[,c(..range_keys,'peakID','significant','species')]
chimp_peaks = read_da_peaks('chimp_da_results.txt',c('significant','non_significant'))[,c(..range_keys,'peakID','significant','species')]

common_regions = read_da_peaks('common_regions.txt',c('significant','non_significant'))[,species:='common'][,c(..range_keys,'peakID','significant','species')]
chimp_common_regions = copy(common_regions)%>%convert_coord('hg38ToPanTro5.over.chain')

## make 2 objects with species-specific + common peaks all in hg38 coordinates
human_common_specific = rbind(human_peaks,common_regions)
chimp_common_specific = rbind(chimp_peaks,chimp_common_regions)

all_peaks = list(chimp_common_specific,human_common_specific)
lapply(all_peaks,function(x)setkeyv(x,range_keys))
names(all_peaks) = species_names[-2]

##---------------------------------
## overlaps of peaks within TADs
##---------------------------------
peak_tads_overlap = purrr::map2(
  all_peaks,all_tads,
  function(x,y)foverlaps(x,y,type='within')%>%
  na.omit()%>%
  setnames(
    old=c(paste('i.',range_keys[-1],sep=''),'i.species'),
    new=c(paste('peak_',range_keys[-1],sep=''),'peak_species')
    )
)

## some numbers and %
lapply(peak_tads_overlap,function(x)x[,c('seqnames','peak_start','peak_end')]%>%unique()%>%nrow())
# $chimp 20290
# $human 24030

fraction_peaks_overalp = purrr::map2(
  peak_tads_overlap,all_peaks,
  function(x,y)
  x[,c(1,6,7)]%>%unique()%>%nrow()
  /
  y[,c(..range_keys)]%>%unique()%>%nrow() *100
)
# $chimp 72.0
# $human 77.0

## calculate the number of peaks per TAD and then plot:
## 1.numb peaks per tad
## 2.corr numb peaks tad width
## 3. distribution TAD sizes
peak_tads_overlap = lapply(
  peak_tads_overlap,function(x)x=split(x,by='tadID')%>%
  lapply(
    function(y)y=y[
      ,numbpeaks_per_tad:=.N
  ]
  )%>%rbindlist()
)

## some QCs for peak-TAD overlap
qc_peaks_in_tads = copy(peak_tads_overlap)%>%
lapply(
  function(x)x=x[,c(..range_keys,'numbpeaks_per_tad')]%>%unique()
)
qc_peaks_in_tads = Map(mutate,qc_peaks_in_tads,species=names(qc_peaks_in_tads))%>%rbindlist()
qc_peaks_in_tads = qc_peaks_in_tads[,tad_width:=end-start][,log10_tads_width:=log10(tad_width)]

pdf(paste0(plot_dir,'/TADs/distribution_numb_peaks_per_tad.pdf',sep=''),width = 7, height = 7)
ggplot(qc_peaks_in_tads,aes(x=numbpeaks_per_tad,fill=species))+
  geom_bar(position='dodge')+
  facet_wrap(species~.,ncol=1)+
  xlab('number peaks per tad')+ylab('Counts')+
  theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()
 
pdf(paste0(plot_dir,'/TADs/pearson_corr_numb_peaks_tad_size.pdf',sep=''),width = 7, height = 7)
qc_peaks_in_tads%>%split(by='species')%>%
lapply(function(x)
  ggscatter(x,
    x = "tad_width",
    y = "numbpeaks_per_tad",
    add = "reg.line",                                 
    conf.int = TRUE,                                
    add.params = list(color = "blue",fill = "lightgray"))+
    stat_cor(method = "pearson", label.x = 3, label.y = 30)
)
dev.off()

## QC distribution TAD sizes between human and chimp
## this could be the reason why we see many more human peaks in TADs
pdf(paste0(plot_dir,'/TADs/distribution_tad_sizes.pdf',sep=''),width = 10, height = 7)
ggplot(qc_peaks_in_tads,aes(x=log10_tads_width,fill=species))+
  geom_density(alpha=0.7)+
  facet_wrap(species~.,ncol=2)+
  xlab('TAD size')+ylab('Density')+
  theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()

##-------------------------------------------
## now look whether peaks and target genes 
## are in the same TAD
##-------------------------------------------
## read target genes
target_genes = list.files(paste(go_dir,'target_genes',sep=''),full.names=T,recursive=F)%>%
lapply(
  function(x)fread(x,sep='\t',header=T)[,c('peakID','gene','distTSS')])
names(target_genes) = gsub("\\_.*","",list.files(paste(go_dir,'target_genes',sep=''),full.names=F,recursive=F))
target_genes = Map(mutate,target_genes,target_gene_species=names(target_genes))

combined_target_genes = list(
  rbindlist(target_genes[c(1:2)]),
  rbindlist(target_genes[c(2:3)])
)
names(combined_target_genes) = species_names[-2]

peak_in_tads_w_genes = purrr::map2(
  peak_tads_overlap,combined_target_genes,
  function(x,y) inner_join(x,y,by=c('peakID','peak_species'='target_gene_species'))
)

## Not all peaks have target genes associated (~ 1%)
## you can use these to calculate the %
lapply(peak_tads_overlap,function(x)x[,c('seqnames','peak_start','peak_end')]%>%unique()%>%nrow())
lapply(peak_in_tads_w_genes,function(x)x[,c('seqnames','peak_start','peak_end')]%>%unique()%>%nrow())

## check if peaks and target genes are in the same TAD
peak_in_tads_w_genes = copy(peak_in_tads_w_genes)%>% ##add this copy otherwise it complains
lapply(function(x)
  x=x[
      ,gene_start := ifelse(distTSS<0,peak_start+distTSS,peak_end+distTSS)
      ][
          ,distTAD := gene_start-start
            ][,.SD[which.min(abs(distTAD))], by=.(peakID)][
              ,peak_genes_same_tad := ifelse(distTAD>=0,'y','n')
              ]
)

## count numb peaks and target genes are in same tad + % of all peaks in tads
qc_peaks_w_genes_in_tads = copy(peak_in_tads_w_genes)%>%
lapply(
  function(x)x=x%>%dplyr::select(c('seqnames',contains('peak_'),'gene','significant','distTSS','distTAD'))%>%unique()
)
qc_peaks_w_genes_in_tads = 
lapply(qc_peaks_w_genes_in_tads,
  function(x)x=x[
    ,numbpeaks_w_genes_in_out_tad:=.N,by=.(peak_genes_same_tad)
    ][
      ,numbpeaks_in_tad:=.N
      ][
        ,prop_peaks_w_genes_in_out_tad:=numbpeaks_w_genes_in_out_tad/numbpeaks_in_tad
      ]
)
qc_peaks_w_genes_in_tads = Map(mutate,qc_peaks_w_genes_in_tads,file=names(qc_peaks_w_genes_in_tads))%>%rbindlist()

pdf(paste0(plot_dir,'/TADs/prop_peaks_w_gene_in_out_tad.pdf',sep=''),width = 7, height = 7)
ggplot(
  qc_peaks_w_genes_in_tads[,c('prop_peaks_w_genes_in_out_tad','peak_genes_same_tad','file')]%>%unique(),
  aes(x=file,y=prop_peaks_w_genes_in_out_tad,fill=peak_genes_same_tad))+
geom_bar(stat='identity')+
xlab(' ')+ylab('Proportion peaks with gene')+
theme(
  legend.key = element_rect(fill = "white", colour = "black"),
  axis.line = element_blank()
)
dev.off()

##----------------------------------------------------------
## Now keep only peaks with target genes in the same TAD
##----------------------------------------------------------
peak_w_genes_same_tad = copy(peak_in_tads_w_genes)%>%
lapply(
  function(x)
  x=x[
    peak_genes_same_tad=='y'
    ][
      ,c('seqnames','peak_start','peak_end','peakID','significant','peak_species','tadID','gene','distTSS')
      ]%>%unique()
      )%>%lapply(
            function(y)y=y[
              ,numbpeaks_per_gene_same_tad:=.N,by=.(gene,tadID)
              ]%>%setnames(old=c(1:3),new=range_keys)%>%unique()
)
numb_peaks_w_genes_same_tad =lapply(peak_w_genes_same_tad,function(x)x[,c(1:3)]%>%unique()%>%nrow())
# $chimp 18189
# $human 21624

numb_genes_w_peaks_same_tad=lapply(peak_w_genes_same_tad,function(x)x[,c('gene')]%>%unique()%>%nrow())
# $chimp 8619
# $human 10107

peak_w_genes_same_tad = Map(mutate,peak_w_genes_same_tad,file=names(peak_w_genes_same_tad))%>%rbindlist()

pdf(paste0(plot_dir,'/TADs/numb_peaks_per_gene_same_tad.pdf',sep=''),width = 7, height = 7)
ggplot(peak_w_genes_same_tad,aes(x=file,y=numbpeaks_per_gene_same_tad,fill=file))+
geom_boxplot(notch=T)+
xlab(' ')+ylab('number peaks per gene in the same TAD')+
stat_compare_means(
  method = "wilcox.test",
  label = "p.signif",
  label.y = max(peak_w_genes_same_tad$numbpeaks_per_gene_same_tad+1),
  size=7
  )+
  theme(
    legend.key = element_rect(fill = "white", colour = "black"),
    axis.line = element_blank()
)
dev.off()

## look at the type of regulation of each gene (i.e. common, common +specific and/or specific)
peak_w_genes_same_tad = peak_w_genes_same_tad%>%split(by='file')%>%
lapply(
  function(x)x%>%split(by='gene')%>%
    lapply(
      function(y)y[
        ,regulation:=ifelse(
          length(unique(y$peak_species)) == 1, unique(y$peak_species),'common_and_specific')
      ]
    )%>%rbindlist()%>%unique()
)

## get the % of genes for each of these regulation groups
peak_w_genes_same_tad = lapply(
  peak_w_genes_same_tad,function(x)x=x[
    ,regulation:=ifelse(regulation %in% c('human','chimp'),'specific',regulation)
    ][
      ,numb_peaks:=.N
      ][
        ,numb_peaks_per_regulation:=.N,by=.(regulation)
        ][
          ,prop_regulation:=numb_peaks_per_regulation/numb_peaks
        ]
)%>%rbindlist()

## plot proportion genes per regulation type
pdf(paste0(plot_dir,'/TADs/prop_genes_per_regulation_type.pdf',sep=''),width = 7, height = 7)
ggplot(
  peak_w_genes_same_tad[,c('file','regulation','prop_regulation')]%>%unique(),
  aes(x=file,y=prop_regulation,fill=regulation))+
geom_bar(stat='identity',position='stack')+
xlab('')+ylab('Proportion genes')+
theme(
  legend.key = element_rect(fill = "white", colour = "black"),
  axis.line = element_blank()
)
dev.off()

## for those set of genes with common and specific regulation
## look at the distribution of distances between common and species-specific peaks
## importantly, there are no genes in common between the sets (within each file of course) 
common_sp_specific =copy(peak_w_genes_same_tad)[regulation=='common_and_specific'][,c('peak_species','gene','distTSS','file')]
common_sp_specific = common_sp_specific[,peak_species:=ifelse(peak_species=='common',peak_species,'specific')]%>%unique()

pdf(paste0(plot_dir,'/TADs/dist_species_specific_common_peaks_from_genes_same_tad.pdf',sep=''),width = 7, height = 7)
ggplot(common_sp_specific,aes(x=peak_species,y=log10(abs(distTSS)+1),fill=peak_species))+
geom_boxplot(notch=T)+
facet_wrap(file~.,ncol=2)+
xlab('')+ylab('log10 abs distance peak from gene')+
stat_compare_means(
  method = "wilcox.test",
  label = "p.signif",
  label.y = max(log10(abs(peak_w_genes_same_tad$distTSS)+1)+0.5),
  size=7,
  ref.group='common'
  )+
theme(
  legend.key = element_rect(fill = "white", colour = "black"),
  axis.line = element_blank()
)
dev.off()

##---------------------
## Add expression info
##---------------------
## Look at the expression level (avg normalised cpm) of the genes regulated by these peaks
## read ensembl gene ids from archive 86
## note that becuase i am using an older ensembl version and because I am focusing 
## only on genes with orthologous i will loose ~27% genes on average
ensembl_86_ids = fread("./rna_seq/de_output/ensembl_id_hugo_symbols.txt",sep='\t',header=T)

## deprecated, take just AvgExp from DE results
# ## get normalised cpms
# get_avg_cmp_counts = function(species_dir){
#   file = dir(paste('./rna_seq/',species_dir,'/output/PostAlignment/',sep=''),pattern='ipsc',full.names=T)
#   ipsc_gene_counts = fread(file,sep='\t',header=T)
  
#   genes=ipsc_gene_counts$GeneID

#   matrix=ipsc_gene_counts[,-1]
#   matrix_dge = edgeR::DGEList(matrix)
#   rownames(matrix_dge$counts)=genes
#   keep_genes = edgeR::filterByExpr(matrix_dge$counts,group=matrix_dge$samples$group)
#   matrix_dge = matrix_dge[keep_genes,, keep.lib.sizes=FALSE]
  
#   matrix_dge$counts = edgeR::cpm(matrix_dge$counts)
#   matrix_dge = edgeR::calcNormFactors(matrix_dge)

#   avg_cpm=copy(matrix_dge$counts)%>%as.data.table()
#   avg_cpm =avg_cpm[,gene:=rownames(matrix_dge$counts)][, .(avg_cpm = rowMeans(.SD)), by = gene]%>%setnames(old=1,new='EnsemblID')
#   return(avg_cpm)
# }

# chimp_avg_cpm = get_avg_cmp_counts('pantro5')
# human_avg_cpm = get_avg_cmp_counts('hg38')

# avg_cpms =list(chimp_avg_cpm,human_avg_cpm)

de_genes = fread("./rna_seq/de_output/topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",sep=' ',header=F,col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B'))

avg_gene_expr = copy(peak_w_genes_same_tad)%>%split(by='file')%>%
lapply(function(x)x=x[
  ensembl_86_ids,on='gene',nomatch=0
  ][
    de_genes,on='EnsemblID',nomatch=0
    ]
)
numb_peaks_w_gene_same_tad_w_expr =lapply(avg_gene_expr,function(x)x=x[,c(..range_keys)]%>%unique()%>%nrow())
# $chimp
# [1] 13569
# $human
# [1] 15690

purrr::map2(numb_peaks_w_gene_same_tad_w_expr,numb_peaks_w_genes_same_tad,`/`)

avg_gene_expr =rbindlist(avg_gene_expr)[,c('regulation','EnsemblID','file','AveExpr')]%>%unique()
# avg_gene_expr = purrr::map2(avg_gene_expr,avg_cpms,function(x,y)x[y,on='EnsemblID',nomatch=0])%>%rbindlist()

## add all genes to include another baseline comparison
all_gene_expr_chimp = copy(de_genes)[,regulation:='all_genes'][,file:='chimp'][,c('EnsemblID','regulation','file','AveExpr')]
all_gene_expr_human = copy(de_genes)[,regulation:='all_genes'][,file:='human'][,c('EnsemblID','regulation','file','AveExpr')]

all_genes_expr =rbind(avg_gene_expr,all_gene_expr_human,all_gene_expr_chimp)
# all_genes_expr = all_genes_expr[,log_avg_cpm:=log(avg_cpm+1)]
comparisons = list(
  c('common','common_and_specific'),
  c('common','specific'),
  c('common_and_specific','specific')
)

pdf(paste0(plot_dir,'/TADs/avg_expression_genes_peaks_same_tad.pdf',sep=''),width = 7, height = 8)
p <- ggboxplot(
  all_genes_expr, x = "regulation", y = "AveExpr",
  fill = "file",notch = T,order=c('all_genes','common','common_and_specific','specific'),
  facet.by = "file", short.panel.labs = FALSE
)+xlab(' ')+ylab('Avg gene expression')
p + stat_compare_means(
  comparisons = comparisons,
  show.legend=FALSE,
  method='wilcox.test', 
  p.adjust.method = "fdr",
  color = "black", size = 7
  ) +
theme(
  legend.position = "none",
  axis.line = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)
dev.off()

## check LOEUF score for these genes and see if there are differences between the regulation types
loeuf_scores = fread('/data/gpfs/projects/punim0586/dvespasiani/Annotation_and_other_files/human_genome/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz',sep='\t',header=T)%>%
setnames(old='gene_id',new='EnsemblID')
all_genes_expr_loeuf = copy(all_genes_expr)[
  loeuf_scores,on='EnsemblID',nomatch=0
  ][
    ,c('regulation','EnsemblID','file','avg_cpm','oe_lof_upper')
]

length(unique(all_genes_expr_loeuf$EnsemblID))
all_genes_expr_loeuf=na.omit(all_genes_expr_loeuf)
length(unique(all_genes_expr_loeuf$EnsemblID))

pdf(paste0(plot_dir,'/TADs/loeuf_scores_genes_peaks_same_tad.pdf',sep=''),width = 7, height = 8)
p <- ggboxplot(
  all_genes_expr_loeuf, x = "regulation", y = "oe_lof_upper",
  fill = "file",notch = T,order=c('all_genes','common','common_and_specific','specific'),
  facet.by = "file", short.panel.labs = FALSE
)+xlab(' ')+ylab('Upper bound of 90% CI for o/e ratio for pLoF variants')
p + stat_compare_means(
  comparisons = comparisons,
  show.legend=FALSE,
  method='wilcox.test', 
  p.adjust.method = "fdr",
  color = "black", size = 7
  ) +
theme(
  legend.position = "none",
  axis.line = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
)

dev.off()

## write down list of peaks + genes in the same TAD for the DE/DA analysis 
## also include the field with numb peaks per gene so that 
## later u can filter for those genes w multiple peaks associated for Gviz 
list_peaks_w_genes_same_tad = copy(peak_w_genes_same_tad)[,c('numb_peaks','numb_peaks_per_regulation','prop_regulation'):=NULL]
write.table(list_peaks_w_genes_same_tad,paste(tmp_files_dir,'list_peaks_w_genes_same_tad.txt',sep=''),sep='\t',col.names=T, quote=F,row.names=F)

