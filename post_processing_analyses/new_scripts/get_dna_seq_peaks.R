## use this script to generate a fasta file containing the sequence
## in the respective species coordinates for each the peak 

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

if (genome =='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    species_bsgenome <- Hsapiens
    }else{
    library(BSgenome.Ptroglodytes.UCSC.panTro5)
    species_bsgenome <- Ptroglodytes
}

output_dir <- create_dir(outdir,'sequences')

## get DA peaks
da_file <- paste(da_dir,genome,'/da_results.txt',sep='')
da_results <- fread(da_file,sep='\t',header=T,select=c(range_keys,'DA','peakID','peak_species'))%>%split(by='peak_species')

## get DNA sequences 
dna_sequences <- lapply(da_results,function(x){ 
    sequence <- copy(x)[
        ,sequence:=as.character(getSeq(species_bsgenome, seqnames,start, end))
        ]
    dna_seq <- DNAStringSet(sequence$sequence)
    names(dna_seq) = sequence$peakID
    return(dna_seq)
})

filenames <- paste(output_dir,names(dna_sequences),'_',genome,'_peak_sequence','.fa',sep='')
mapply(writeXStringSet, dna_sequences, file = filenames,append=F,compress=F, compression_level=NA, format="fasta")

