library(ATACseqQC)
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(csaw)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'atac_seq_qc')

encode_ideal = 7
encode_acceptable = 5

## list bam
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes
param <- readParam(pe = "both",restrict=standard_chr,minq=20, dedup=TRUE)

get_bams <- function(species){
    bams <- list.files(paste0('../',genome,bamDir,sep=''), recursive = T,full.names = T,pattern="^H.*_tn5_shifted_sorted.bam$|C.*_tn5_shifted_sorted.bam$")
    return(bams)
}

bams <- get_bams(genome)

alignments <- lapply(bams, function(x)readGAlignments(x))
names(alignments) = samples_names

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

get_tsse <- function(alignment){
    tsse <- TSSEscore(alignment, txs)
    range_bp <- 100*(-9:10-.5)
    tsse_df <- data.table(range=range_bp,tsse_enrichment=tsse$values,tsse_score=tsse$TSSEscore)
    return(tsse_df)
}

tss_enrichment <- lapply(alignments,function(x)get_tsse(x))
tss_enrichment <- Map(mutate,tss_enrichment,samples=samples_names)%>%rbindlist()

pdf(paste(outplot_dir,'tss_enrichment.pdf',sep=''),width=7,height=7)
ggplot(tss_enrichment,aes(x=range,y=tsse_enrichment,col=samples))+
    geom_line(size=1)+
    geom_hline(yintercept=encode_acceptable)+
    geom_text(aes(-900,encode_acceptable,label = 'acceptable', vjust = -1),col='black')+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(-900,encode_ideal,label = 'ideal', vjust = -1),col='black')+
    ylab('Enrichment')+xlab('Distance from TSS')+
    # scale_color_manual(values=samples_palette)+
    theme(
        legend.position='bottom',
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()




