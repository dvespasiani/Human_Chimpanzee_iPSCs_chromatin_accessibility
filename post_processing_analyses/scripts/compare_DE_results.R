## use this script to compare my DE results with Irene's 2015

library(dplyr)
library(data.table)
library(magrittr)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))

de_dir = './rna_seq/de_output/'
##-----------------
## read DE genes 
##-----------------
## my DE results 
my_results = fread(paste(de_dir,"topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",sep=''),sep=' ',header=F,
col.names=c('genes','EnsemblID','logFC','AveExpr','t','P.Value','adj.P.Val','B'))[,significant:=ifelse(adj.P.Val<=0.05,'s','ns')]
## Irene's
igr_results= fread(paste(de_dir,"gallego_romero_etal_2015_de_genes.txt",sep=''),sep='\t',header=T)[,significant:=ifelse(adj.P.Val<=0.05,'s','ns')]

length(unique(my_results$EnsemblID))
# [1] 13094
length(unique(igr_results$EnsemblID))
# [1] 12171

## check how many genes in common and then how many DE genes in common
common_results=copy(my_results)[igr_results,on='EnsemblID',nomatch=0]
length(unique(common_results$EnsemblID))
# [1] 10970

length(unique(common_results$EnsemblID))/length(unique(igr_results$EnsemblID))
# 90.1%

## for DE genes
my_de=copy(my_results)[significant=='s']
igr_de=copy(igr_results)[significant=='s']

length(unique(my_de$EnsemblID))
# [1] 6275
length(unique(igr_de$EnsemblID))
# [1] 6204

## check how many genes in common and then how many DE genes in common
common_de_results=copy(my_de)[igr_de,on='EnsemblID',nomatch=0]
length(unique(common_de_results$EnsemblID))
# [1] 4561

length(unique(common_de_results$EnsemblID))/length(unique(igr_de$EnsemblID))
# 73.5%

## look at correlation between results
cor.test(common_results$logFC,common_results$i.logFC,method='spearman')
# data:  common_results$logFC and common_results$i.logFC
# S = 1.8217e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.9172045


cor.test(common_de_results$logFC,common_de_results$i.logFC,method='spearman')
# data:  common_de_results$logFC and common_de_results$i.logFC
# S = 693628486, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#      rho 
# 0.956137 

