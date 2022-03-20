##--------------------------------------------------------------------------------------------------
## script to create random set of genomic regions using regioneR R package
## to run this script do:
## Rscript ./make_random_genomic_regions.R new_da_results.txt random_genomic_regions.txt 3
##--------------------------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(regioneR)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))
set.seed(2022)

wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/'
outdir = 'output/random_regions/'

dir.create(file.path(wd, outdir), showWarnings = FALSE)

args = commandArgs(trailingOnly=TRUE)

input=args[1]
output=paste(wd,outdir,args[2],sep='')
numbsets = args[3]

## get DA peaks
da_results <- read_da_results(input)
da_results <- da_results[,c(..range_keys)]%>%makeGRangesFromDataFrame()

blacklist <- fread(paste(wd,'../data/ENCODE_blacklisted/hg38_blacklist_v2.bed',sep=''),sep='\t',header=F,col.names=range_keys)%>%makeGRangesFromDataFrame()

random_gr <-list()

for (i in 1:numbsets){
    print(paste('Creating random region numb:',i,sep=' '))
    random_gr[[i]] <- randomizeRegions(da_results,mask=blacklist,genome='hg38',allow.overlaps=F)%>%as.data.table()
}
random_regions <- rbindlist(random_gr)
random_regions <- random_regions[,c('width','strand'):=NULL][seqnames %in% standard_chr]%>%setorderv(c(range_keys),1)

print(paste('Saving file to ',output,sep=''))
write.table(random_regions,output,sep='\t',col.names=T,row.names=F,quote=F)

