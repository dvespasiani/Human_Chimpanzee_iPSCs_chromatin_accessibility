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
da_file <- paste(da_dir,genome,'/da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))

all_peaks <- split(da_results,by='peak_species')%>%lapply(function(x)x=x[,c(..range_keys,'peakID')])

## using peaks after DA analysis (both da and non DA)
arrange_columns_for_homer = function(x,othercols){
    x=x[
    ,strand:= '+'][
            ,c(..range_keys,'peakID','strand',..othercols)]
}

homer_input_format <- copy(all_peaks)%>%lapply(
    function(x){
        x<-x[
            ,strand:= '+'
            ][
                ,c(..range_keys,'peakID','strand')
                ]
    }
)

filenames <-  paste0(output_dir,names(all_peaks),'_',genome,'_homer_input.txt',sep='')
mapply(write.table, all_peaks, file = filenames,sep='\t',col.names=T,row.names=F,quote=F)



