## use this script to get and plot a summary of alignment results
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd('/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')


create_dir <- function(genome,step){
    outdir <- paste(genome,'/output/',step,'/qc/',sep='')
    dir.create(outdir,showWarnings = FALSE)
    return(outdir)
}

get_alignment_qc = function(genome){
    files <- list.files(paste(genome,'/output/Alignment/logs',sep=''),full.names=T,recursive=F,pattern='bowtie_Align_qc.log')
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

hg38_alignment_qcs <- get_alignment_qc('hg38')
pantro5_alignment_qcs <- get_alignment_qc('panTro5')

plot_alignment_rate = function(x,label_y){
    encode_acceptable = 0.8
    encode_ideal = 0.95
    p = ggplot(x,aes(x=sample,y=fraction_reads,fill=class))+ 
    geom_bar(stat="identity")+xlab('')+ylab(label_y)+
    geom_hline(yintercept=encode_acceptable)+
    geom_text(aes(1,encode_acceptable,label = 'acceptable', vjust = -1))+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(1,encode_ideal,label = 'ideal', vjust = -1))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    return(p)
}

pdf(paste(create_dir('hg38','Alignment'),'hg38_bowtie_alignment_rate.pdf',sep=''),width=10,height=7)
plot_alignment_rate(hg38_alignment_qcs,'overall alignment rate')
dev.off()

pdf(paste(create_dir('panTro5','Alignment'),'panTro5_bowtie_alignment_rate.pdf',sep=''),width=10,height=7)
plot_alignment_rate(pantro5_alignment_qcs,'overall alignment rate')
dev.off()

## FRiP
get_frip <- function(genome){
    files <- list.files(paste(genome,'/output/PeakCalling/qc',sep=''),full.names=T,recursive=F,pattern='.txt')
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

hg38_frip <- get_frip('hg38')
pantro5_frip <- get_frip('panTro5')


plot_frip = function(x){
    encode_acceptable = 0.2
    encode_ideal = 0.3
    p = ggplot(x,aes(x=sample,y=frip,fill=type))+ 
    geom_bar(stat='identity')+xlab('')+
    geom_hline(yintercept=encode_acceptable)+
    geom_text(aes(1,encode_acceptable,label = 'acceptable', vjust = -1))+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(1,encode_ideal,label = 'ideal', vjust = -1))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    return(p)
}

pdf(paste(create_dir('hg38','PeakCalling'),'hg38_frip.pdf',sep=''),width=10,height=7)
plot_frip(hg38_frip)
dev.off()

pdf(paste(create_dir('panTro5','PeakCalling'),'panTro5_frip.pdf',sep=''),width=10,height=7)
plot_frip(pantro5_frip)
dev.off()