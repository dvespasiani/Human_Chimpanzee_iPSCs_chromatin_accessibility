## use this script to get and plot a summary of alignment results
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd('/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses')

scripts_dir <- './scripts/'
source(paste(scripts_dir,'utils.R',sep=''))

outplot_dir <- create_dir(plot_dir,'atac_seq_qc')

get_alignment_qc = function(genome){
    files <- list.files(paste('../',genome,'/output/Alignment/logs',sep=''),full.names=T,recursive=F,pattern='bowtie_Align_qc.log')
    qc_results <- lapply(files,function(x) {
        x <- fread(x,sep='\t',col.names='bowtie_output')
        x<-x[
        c(3,14),
        ][
            ,oar:=ifelse(bowtie_output %like% 'overall',
            gsub("\\%.*","",bowtie_output), 
            gsub("\\%.*","",gsub("[\\(\\)]", "", regmatches(bowtie_output, gregexpr("\\(.*?\\)", bowtie_output)))))
            ][
                ,fraction_reads:=as.numeric(oar)/100
                ][
                    ,class:=ifelse(bowtie_output %like% 'overall','overall_alignment_rate','uniquely_mappable')
                    ][
                        ,c('fraction_reads','class')
                        ][
                            ,metric:=ifelse(class %like% 'overall','alignment_rate','mappability')
                        ]
                        }
    )
    sample_names <- gsub(".*/","",gsub('\\_.*','',files))
    qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()
    alignment_rate <- copy(qc_results)[metric %like% 'alignment']

    return(alignment_rate)
}

hg38_alignment_qcs <- get_alignment_qc('hg38')[sample%like%'H']
pantro5_alignment_qcs <- get_alignment_qc('panTro5')[sample%like%'C']

alignment_qc <- rbind(pantro5_alignment_qcs,hg38_alignment_qcs)

plot_alignment_rate = function(x){
    encode_acceptable = 0.8
    encode_ideal = 0.95
    p = ggplot(x,aes(x=sample,y=fraction_reads,fill=sample))+ 
    geom_bar(stat="identity")+xlab('')+ylab('Bowtie2 Alignment rate')+
    geom_hline(yintercept=encode_acceptable)+
    scale_fill_manual(values=samples_palette)+
    geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
    geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    theme_classic()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        )
    return(p)
}

pdf(paste(outplot_dir,'bowtie2_alignment_qc.pdf',sep=''),width=7,height=7)
plot_alignment_rate(alignment_qc)
dev.off()

# pdf(paste(create_dir('panTro5','Alignment'),'panTro5_bowtie_alignment_rate.pdf',sep=''),width=10,height=7)
# plot_alignment_rate(pantro5_alignment_qcs,'overall alignment rate')
# dev.off()

## FRiP
get_frip <- function(genome){
    files <- list.files(paste('../',genome,'/output/PeakCalling/qc',sep=''),full.names=T,recursive=F,pattern='.txt')
    files <- stringr::str_subset(files,pattern="merged",negate = T)
    qc_results <- lapply(files,function(x){
        x<-fread(x,sep='\t',header=T,col.names='results')
        x<-x[2,][,frip:=round(as.numeric(gsub('.*=','',results)),2)]
        }
    )
    sample_names <- gsub(".*/","",gsub('\\_.*','',files))
    qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()
    qc_results <- qc_results[,c('frip','sample')][,type:='reads_in_peak']
    return(qc_results)
}

hg38_frip <- get_frip('hg38')[sample%like%'H']
pantro5_frip <- get_frip('panTro5')[sample%like%'C']

frip_qc <- rbind(pantro5_frip,hg38_frip)


plot_frip = function(x){
    encode_acceptable = 0.2
    encode_ideal = 0.3
    p = ggplot(x,aes(x=sample,y=frip,fill=sample))+ 
    geom_bar(stat='identity')+xlab('')+ylab('FRiP score')+
    geom_hline(yintercept=encode_acceptable)+
    geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))+
    scale_fill_manual(values=samples_palette)+
    theme_classic()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        )
    return(p)
}

pdf(paste(outplot_dir,'frip_qc.pdf',sep=''),width=7,height=7)
plot_frip(frip_qc)
dev.off()
