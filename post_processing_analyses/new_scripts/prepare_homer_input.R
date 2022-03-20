library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

output_dir <- create_dir(outdir,'homer/homer_input')

## get DA peaks
da_results <- read_da_results('new_da_results.txt')
da_results <- da_results[,c(..range_keys,'DA','peakID','FDR','peak_species','logFC')]

# all_peaks <- split(da_results,by='DA')[,c(..range_keys,'peakID')]

## using peaks after DA analysis (both da and non DA)
# arrange_columns_for_homer = function(x,othercols){
#     x=x[
#     ,strand:= '+'][
#             ,c(..range_keys,'peakID','strand',..othercols)]
# }

homer_input_format <- copy(da_results)[
            ,strand:= '+'
            ][
                ,c('peakID',..range_keys,'strand')
]

filenames <-  paste0(output_dir,'allpeaks_',genome,'_homer_input.txt',sep='')
write.table(homer_input_format, file = filenames,sep='\t',col.names=T,row.names=F,quote=F)



