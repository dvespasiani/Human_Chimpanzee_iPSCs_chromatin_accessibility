## use this script to generate a fasta file containing the sequence
## in the respective species coordinates for each the peak 
##-------------------------------------------------------------------------------------------------------------
## run script from command line as 
## Rscript ./get_dna_seq_peaks.R input output
## e.g., Rscript ./scripts/get_dna_seq_peaks.R da_results.txt output/sequences/all_peaks_seq.fa
## e.g., Rscript ./scripts/get_dna_seq_peaks.R output/random_regions/random_genomic_regions.txt output/sequences/random_regions_seq.fa
##-------------------------------------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)

options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir = './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

args = commandArgs(trailingOnly=TRUE)

wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/'
input=args[1]
output=paste(wd,args[2],sep='')

if (genome =='hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    species_bsgenome <- Hsapiens
    }else{
    library(BSgenome.Ptroglodytes.UCSC.panTro5)
    species_bsgenome <- Ptroglodytes
}

## get input peaks
# peaks <- read_da_results(input)

if (input %like% 'da_results'){
    print(paste('reading in ', input,sep=''))
    peaks <- read_da_results(input)
    }else{
    peaks <- fread(input,sep='\t',header=T)
}

if ('peakID' %in% names(peaks)){
    peaks <- peaks
}else{
    peaks<-peaks[
        ,peakID:=paste('peak_',1:nrow(peaks),sep='')
    ]
}

## get DNA sequences 
dna_sequences <- copy(peaks)[,sequence:=as.character(getSeq(species_bsgenome, seqnames,start, end))]
dna_seq <- DNAStringSet(dna_sequences$sequence)
names(dna_seq) = dna_sequences$peakID

writeXStringSet(dna_seq, output, append=FALSE,compress=FALSE, compression_level=NA, format="fasta")
     
