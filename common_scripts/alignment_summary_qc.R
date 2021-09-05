## use this script to get and plot a summary of alignment results
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd('/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

plot_dir = 'post_processing_analyses/output/plots/Snakemake_QCs/'

get_output=function(species_dir){
    qc = list.files(paste(species_dir,'/output/Alignment/QCs',sep=''),full.names=T,recursive=F,pattern='bowtie_Align_qc.log')
    qc = lapply(qc,function(x)fread(x,sep='\t',col.names='bowtie_output'))
    qc = lapply(qc,function(x)x=x[
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
                        )
    sra_names=list.files(paste(species_dir,'/output/Alignment/QCs',sep=''),full.names=F,recursive=F,pattern='bowtie_Align_qc.log')
    sra_names =gsub('\\_.*','',sra_names)
    qc = Map(mutate,qc,SRA=sra_names)%>%rbindlist()
    # chimp_sra = paste('SRR81764',31:41,sep='')
    # qc=qc[,species:=ifelse(SRA %in% chimp_sra,'chimpanzee','human')]
    mappability=copy(qc)[metric %like% 'mappability' ]
    one_minus_mappability = copy(mappability)[,fraction_reads:=ifelse(class %like% 'overall',fraction_reads,1-fraction_reads)][,class:='non_uniquely_mappable']
    mappability_combined=rbind(mappability,one_minus_mappability)

    alignment_rate=copy(qc)[metric %like% 'alignment']

    return(list(alignment_rate,mappability_combined))
}

hg38_alignment_qcs=get_output('./hg38')
pantro5_alignment_qcs=get_output('pantro5')

plot_alignment_rate = function(x,label_y){
    p = ggplot(x,aes(x=SRA,y=fraction_reads,fill=class))+ 
    geom_bar(stat="identity")+xlab('')+ylab(label_y)+
    # geom_hline(yintercept=0.95)+
    # geom_hline(yintercept=0.8,linetype="dashed")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    return(p)
}


pdf(paste(plot_dir,'hg38_bowtie_alignment_rate.pdf',sep=''),width=10,height=7)
plot_alignment_rate(hg38_alignment_qcs[[1]],'overall alignment rate')+
    geom_hline(yintercept=0.95)+
    geom_hline(yintercept=0.8,linetype="dashed")
dev.off()

pdf(paste(plot_dir,'pantro5_bowtie_alignment_rate.pdf',sep=''),width=10,height=7)
plot_alignment_rate(pantro5_alignment_qcs[[1]],'overall alignment rate')+
    geom_hline(yintercept=0.95)+
    geom_hline(yintercept=0.8,linetype="dashed")
dev.off()

pdf(paste(plot_dir,'hg38_bowtie_mappability.pdf',sep=''),width=10,height=7)
plot_alignment_rate(hg38_alignment_qcs[[2]],'fraction reads')
dev.off()

pdf(paste(plot_dir,'pantro5_bowtie_mappability.pdf',sep=''),width=10,height=7)
plot_alignment_rate(pantro5_alignment_qcs[[2]],'fraction reads')
dev.off()


## FRiP
get_frip = function(species_dir){
    qc = list.files(paste(species_dir,'/output/PeakCalling/QC',sep=''),full.names=T,recursive=F)
    qc = lapply(qc,function(x)
    fread(x,sep='\t',header=T,col.names='results')[
        2,
        ][
            ,frip:=round(as.numeric(gsub('.*=','',results)),2)
            ])
    sample_names=list.files(paste(species_dir,'/output/PeakCalling/QC',sep=''),full.names=F,recursive=F)
    sample_names =gsub('\\_.*','',sample_names)
    qc = Map(mutate,qc,sample=sample_names)%>%rbindlist()
    qc = qc[,c('frip','sample')][,type:='reads_in_peak']
    qc_outpeak=copy(qc)[,frip:=1-frip][,type:='reads_non_in_peak']
    
    qc_final=rbind(qc,qc_outpeak)
    return(qc_final)
}

hg38_frip=get_frip('hg38')
pantro5_frip=get_frip('pantro5')

pdf(paste(plot_dir,'hg38_frip.pdf',sep=''),width=10,height=7)
ggplot(hg38_frip,aes(x=sample,y=frip,fill=type))+ 
geom_bar(position="fill", stat="identity")+xlab('')+
geom_hline(yintercept=0.7)+
geom_hline(yintercept=0.8,linetype="dashed")+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf(paste(plot_dir,'pantro5_frip.pdf',sep=''),width=10,height=7)
ggplot(pantro5_frip,aes(x=sample,y=frip,fill=type))+ 
geom_bar(position="fill", stat="identity")+xlab('')+
geom_hline(yintercept=0.7)+
geom_hline(yintercept=0.8,linetype="dashed")+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()