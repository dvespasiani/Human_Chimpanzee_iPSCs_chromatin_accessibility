## script used to get annotate:
## orth peaks, for QC that the peaks are actually meaningful

library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ChIPpeakAnno)
library(org.Pt.eg.db)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Ptroglodytes.UCSC.panTro5)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Ptroglodytes.UCSC.panTro5.refGene)
library(biomaRt)
library(liftOver)
library(ggthemes)
library(ggplot2)
library(ggpubr)

options(width=200)

setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/')

result_dir = './post_processing_analyses/output/'
orth_peakDir = "output/PeakCalling/orthologous_peaks/"
plot_dir  = './post_processing_analyses/output/plots/peak_annotation/'

range_cols = c('seqnames','start','end')

scripts_dir = './post_processing_analyses/scripts/'
source(paste(scripts_dir,'reusable_functions.R',sep=''))


## read orth peaks
read_peaks = function(file){
  peak = fread(file,sep='\t',header=T)
  return(peak)
}

human_ort_peaks = read_peaks(paste('./hg38/',orth_peakDir,'orthologous_peaks.bed',sep=''))%>%makeGRangesFromDataFrame(keep.extra.columns=T)
chimp_ort_peaks = read_peaks(paste('./pantro5/',orth_peakDir,'orthologous_peaks.bed',sep=''))%>%makeGRangesFromDataFrame(keep.extra.columns=T)
 
##-----------------------------------------------------
## 1) Annotation orthologous peaks
## This is done to QC the peaks called by macs2
## check the distribution around TSSs in particular
##-----------------------------------------------------

orth_peaks = list(chimp_ort_peaks,human_ort_peaks)
names(orth_peaks) = c('chimp','human') 

## get list of pantro3 genes (beware that pantro3=patro5)
## keep only genes on standard chromosomes
pantro_mart <- useMart(biomart="ensembl", dataset="ptroglodytes_gene_ensembl")
human_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## Annotations
pantro3_tss_annotation  <- getAnnotation(pantro_mart, featureType="TSS")
pantro3_tss_annotation <- pantro3_tss_annotation[seqnames(pantro3_tss_annotation) %in% c(1:23, "X", "Y")] ## keep only standard chrs (also for human)
TSS.chimp.panTro5 <- pantro3_tss_annotation ## yes, i know....

data(TSS.human.GRCh38)
TSS.human.GRCh38 <- TSS.human.GRCh38[seqnames(TSS.human.GRCh38) %in% c(1:23, "X", "Y")]

## Annotate all peaks with respect to TSS
hg38_anno <- annotatePeakInBatch(
  orth_peaks[[2]], 
  AnnotationData=TSS.human.GRCh38,maxgap = 5000L)%>%
  addGeneIDs(orgAnn="org.Hs.eg.db",feature_id_type='ensembl_gene_id',
  IDs2Add='symbol')

pantro5_anno <- annotatePeakInBatch(
  orth_peaks[[1]], 
  AnnotationData=TSS.chimp.panTro5,maxgap = 5000L)%>%
  addGeneIDs(orgAnn="org.Pt.eg.db",feature_id_type='ensembl_gene_id',
  IDs2Add='symbol')


## Plot some QC plots for all orth peaks

## 1) pie chart with % of overlap between peaks and TSS
pdf(paste0(plot_dir,'human_peak_ovelap_features.pdf',sep=''),width = 7, height = 7)
pie1(table(hg38_anno$insideFeature)) # change it to hg38_anno
dev.off()

pdf(paste0(plot_dir,'chimp_peak_ovelap_features.pdf',sep=''),width = 7, height = 7)
pie1(table(pantro5_anno$insideFeature))
dev.off()

## 2) Distribution of peaks around TSS
plot_distributon_around_TSS = function(peaks,annot_data){
  p <- binOverFeature(peaks, annotationData=annot_data,FUN=length)+
   ylab("number peaks")+xlab('distance from nearest TSS')
   }
pdf(paste0(plot_dir,'human_distribution_peaks_near_TSS.pdf',sep=''),width = 7, height = 7)
plot_distributon_around_TSS(hg38_anno, TSS.human.GRCh38)
dev.off()

pdf(paste0(plot_dir,'chimp_distribution_peaks_near_TSS.pdf',sep=''),width = 7, height = 7)
plot_distributon_around_TSS(pantro5_anno,TSS.chimp.panTro5)
dev.off()

##----------------------------------------
## Now get orth tss and separate peaks
## into those near orth tss and those not
##----------------------------------------
orth_genes = getBM( attributes = c('ptroglodytes_homolog_ensembl_gene','ensembl_gene_id'), ## get all orth genes
  mart = human_mart)%>%as.data.table()

orth_genes = orth_genes[,orth:=ifelse(ptroglodytes_homolog_ensembl_gene=='','no_ort','orth')]

human_genes_with_orth = orth_genes[orth=='orth']$ensembl_gene_id
chimp_genes_with_orth = orth_genes[orth=='orth']$ptroglodytes_homolog_ensembl_gene

add_orth_info = function(x,orth_genes){
  df =copy(x)%>%as.data.table()
  df = df[,has_orth_gene:= ifelse(feature%in%orth_genes,'yes_orth','no_orth')]
  df_orth = copy(df)[has_orth_gene=='yes_orth']
  df_no_orth = copy(df)[has_orth_gene=='no_orth']
  all_peaks = list(df_orth,df_no_orth)%>%lapply(function(x)x=makeGRangesFromDataFrame(x,keep.extra.columns=T))
  return(all_peaks)

}

pantro5_anno_orth = add_orth_info(pantro5_anno,chimp_genes_with_orth)
lapply(pantro5_anno_orth,function(y)count_peaks(y)) 

hg38_anno_orth = add_orth_info(hg38_anno,human_genes_with_orth)
lapply(hg38_anno_orth,function(y)count_peaks(y))

##--------------------------------------------------
## Convert all chimp peaks into hg38 coord
##--------------------------------------------------
chains_path <- './data/LiftOver_chains/'
pantro5_anno_hg38coord = lapply(pantro5_anno_orth,function(x)liftPeaks(x,'panTro5ToHg38.over.chain'))
names(pantro5_anno_hg38coord)=c('yes_orth','no_orth')

## Check the distribution of both peaks near orth tss
## and the remaining peaks around genomic features 
## use hg38 genomic annotations as they are better than the pantro5 ones
annotate_peaks_genome = function(x,organism){
  df =copy(x) %>% lapply(function(y)y=makeGRangesFromDataFrame(y,keep.extra.columns=F))
  df=lapply(df,function(y)y=assignChromosomeRegion(
    y, nucleotideLevel=FALSE,
    precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"), 
    TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
    )
  df = lapply(df,function(y)y=as.data.table(y$percentage))
}

human_peaks_distribution <- annotate_peaks_genome(hg38_anno_orth,'human')
chimp_peaks_hg38_distribution <- annotate_peaks_genome(pantro5_anno_hg38coord,'chimp')
  
plot_genome_peak_anno = function(df1,df2){
  df = rbind(df1[,species:='human'],df2[,species:='chimp'])
  colnames(df)[1:3]=c('Genomic_element','Fraction_peaks','x_axis')
  plot<-ggplot(df,aes(x=x_axis,y=Fraction_peaks,fill=Genomic_element))+
    geom_bar(position="fill", stat="identity")+
     scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9",'darkorchid2','chocolate2','linen','firebrick1'))+
     xlab(' ')+
    theme(
      panel.background =element_rect(fill = 'white', colour = 'black',size=1),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
     axis.line = element_blank())
    return(plot)
}

pdf(paste0(plot_dir,'Distribution_peaks_near_orth_TSS_across_genomic_elements.pdf',sep=''),width = 7, height = 7)
plot_genome_peak_anno(human_peaks_distribution[[1]],chimp_peaks_hg38_distribution[[1]])
dev.off()

pdf(paste0(plot_dir,'Distribution_peaks_NOT_near_orth_TSS_across_genomic_elements.pdf',sep=''),width = 7, height = 7)
plot_genome_peak_anno(human_peaks_distribution[[1]],chimp_peaks_hg38_distribution[[2]])
dev.off()

## GO on orth genes only to confirm results are consistent
GO_enrichment = function(ensemblID,orgAnn){
  df <-getEnrichedGO(
  ensemblID, 
  orgAnn=orgAnn, 
  multiAdjMethod="BH",
  condense=FALSE)%>%
  lapply(function(x)as.data.table(x))
  df_bp = df[[1]]
  
  return(df_bp)
}

orth_genes_GO <- GO_enrichment(orth_genes[orth=='orth']$ensembl_gene_id,"org.Hs.eg.db")

human_non_orth_GO <- GO_enrichment(hg38_anno_orth[[2]]$feature, "org.Hs.eg.db")
chimp_non_orth_GO <- GO_enrichment(pantro5_anno_orth[[2]]$feature, "org.Pt.eg.db")

## make plot of top 30 enriched terms
go_enrichment_plot=function(x,species){
  df<-copy(x)
  df<-df[,c('go.term','BH.adjusted.p.value')]%>%unique()
  
  if(nrow(df)>30){
    df<-df[1:30,]
    df<-df[,species:=species]
  }else{
    df = df[,species:=species]
  }

  plot <- ggplot(df, aes(x=reorder(go.term,-log10(BH.adjusted.p.value)), y=-log10(BH.adjusted.p.value),fill=species)) +
    geom_bar(stat = 'identity',position = 'dodge')+
    xlab(" ") +
    ylab("\n -Log10 BH adj (P) \n ") +
    scale_y_continuous(breaks = round(seq(0, max(-log10(df$BH.adjusted.p.value))), 1)) +
    # scale_fill_manual(name= " ",values=color)+
   theme(
      legend.position='none',
      legend.background=element_rect(),
      axis.text.x=element_text(angle=0, hjust=1.10),
      axis.text.y=element_text(angle=0, vjust=0.8),
      axis.title=element_text(),
      axis.line = element_line(color = "black",size = 0, linetype = "solid"),
     panel.background =element_rect(fill = 'white', size = 0.5,colour = 'black'),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      title=element_text()) +
    coord_flip()
    
    return(plot)
  
}

pdf(paste0(plot_dir,'human_non_orth_GO.pdf',sep=''),width = 7, height = 7)
go_enrichment_plot(human_non_orth_GO,'human')
dev.off()

pdf(paste0(plot_dir,'chimp_non_orth_GO.pdf',sep=''),width = 7, height = 7)
go_enrichment_plot(chimp_non_orth_GO,'chimp')
dev.off()

pdf(paste0(plot_dir,'orth_genes_GO.pdf',sep=''),width = 7, height = 7)
go_enrichment_plot(orth_genes_GO,'chimp')
dev.off()

## count and plot the number/proportion of peaks nearby orth/non orth genes
count_annotations = function(x,species){
  count = function(y,species){
    counts = y[,numb:=.N,by=has_orth_gene][,totnumb:=.N][,proportion:=numb/totnumb][,species:=species]
    return(counts)
    }

  df = copy(x)%>%lapply(function(y)y=as.data.table(y))%>%rbindlist()
  df_peaks = copy(df)[,c('seqnames','start','end','has_orth_gene')]%>%unique()%>%count(species)
  return(df_peaks)
}

pantro5_anno_hg38coord = Map(mutate,pantro5_anno_hg38coord,'has_orth_gene'=names(pantro5_anno_hg38coord))

chimp_counts = count_annotations(pantro5_anno_hg38coord,'chimp')
human_counts = count_annotations(hg38_anno_orth,'human')


plot_counts = function(x,y){
  all_species =rbind(x,y)
  all_species = all_species[,c('has_orth_gene','numb','totnumb','proportion','species')]%>%unique()

  raw_counts_plot = ggplot(all_species,aes(x=species, y=numb, fill=has_orth_gene)) + 
    geom_bar(position="stack", stat="identity")
  
  proportions_plot = ggplot(all_species,aes(x=species, y=proportion, fill=has_orth_gene)) + 
    geom_bar(position="stack", stat="identity")
    
  arranged_plots = ggarrange(raw_counts_plot, proportions_plot, ncol = 2,common.legend = TRUE)

  return(arranged_plots)
}

pdf(paste0(plot_dir,'all_species_numb_peaks_near_orth_non_orth_genes.pdf',sep=''),width = 7, height = 7)
plot_counts(human_counts,chimp_counts)
dev.off()


# ##------------------------
# ## 2) Annotate DA peaks
# ##------------------------
# da_peaks_dir = './post_processing_analyses/output/DA/peaks/'

# read_da_peaks = function(file){
#     df = fread(paste(da_peaks_dir,file,sep=''),sep='\t',header=T,select=c(range_cols,'FDR','logFC','direction'))
#     df = df[,peak_type:=ifelse(FDR < 0.01,'da','non_da')]%>%unique()
#     setkeyv(df,range_cols)
#     return(df)
# }

# human_all_tmm_peaks = read_da_peaks('all_tmm_human_hg38coord_peaks.txt')
# human_all_tmm_peaks = human_all_tmm_peaks[peak_type=='da'][,c(..range_cols,'direction')][,species:='human']

# human_all_tmm_peaks%>%count_peaks()

# chimp_all_tmm_peaks = read_da_peaks('all_tmm_chimp_pantro5_coord_peaks.txt')
# chimp_all_tmm_peaks = chimp_all_tmm_peaks[peak_type=='da'][,c(..range_cols,'direction')][,species:='chimp']

# chimp_all_tmm_peaks%>%count_peaks()

# ## annotate DA peaks with respect to TSSs and add orthologous info
# human_da_peaks <- annotatePeakInBatch(
#   makeGRangesFromDataFrame(human_all_tmm_peaks,keep.extra.columns=T), 
#   AnnotationData=TSS.human.GRCh38,maxgap = 5000L)%>%
#   addGeneIDs(orgAnn="org.Hs.eg.db",feature_id_type='ensembl_gene_id',
#   IDs2Add='symbol')%>%add_orth_info(human_genes_with_orth)

# as.data.table(human_da_peaks[[1]])%>%count_peaks()

# chimp_da_peaks <- annotatePeakInBatch(
#   makeGRangesFromDataFrame(chimp_all_tmm_peaks,keep.extra.columns=T), 
#   AnnotationData=TSS.chimp.panTro5,maxgap = 5000L)%>%
#   addGeneIDs(orgAnn="org.Pt.eg.db",feature_id_type='ensembl_gene_id',
#   IDs2Add='symbol')%>%add_orth_info(chimp_genes_with_orth)

# as.data.table(chimp_da_peaks[[1]])%>%count_peaks()

# ## convert DA chimp to hg38coords
# chimp_da_peaks_hg38 = lapply(chimp_da_peaks,function(x)x=liftPeaks(x,'panTro5ToHg38.over.chain'))

# ## annotate peaks across genomic elements
# human_da_peaks_annot = annotate_peaks_genome(human_da_peaks)
# human_da_peaks_annot  = lapply(human_da_peaks_annot,function(x)x=x[,species:='human'])

# chimp_da_peaks_hg38_annot = annotate_peaks_genome(chimp_da_peaks_hg38)
# chimp_da_peaks_hg38_annot  = lapply(chimp_da_peaks_hg38_annot,function(x)x=x[,species:='chimp'])


# ## Plot distribution DA peaks across genomic elements

# ## near ort Tss
# pdf(paste0(plot_dir,'Distribution_DA_peaks_near_orth_genes_across_genomic_elements.pdf',sep=''),width = 7, height = 7)
# plot_genome_peak_anno(human_da_peaks_annot[[1]],chimp_da_peaks_hg38_annot[[1]])
# dev.off()

# ## and not near orth Tss
# pdf(paste0(plot_dir,'Distribution_DA_peaks_not_near_orth_genes_across_genomic_elements.pdf',sep=''),width = 7, height = 7)
# plot_genome_peak_anno(human_da_peaks_annot[[2]],chimp_da_peaks_hg38_annot[[2]])
# dev.off()


# ## GO enrichment for orth genes near DA peaks
# chimp_da_orth_genes_GO <- GO_enrichment(chimp_da_peaks[[1]]$feature, "org.Pt.eg.db") ## NO ENRICHMENT
# hg38_da_orth_genes_GO <- GO_enrichment(human_da_peaks[[1]]$feature, "org.Hs.eg.db")


# pdf(paste0(plot_dir,'human_DA_orth_GO.pdf',sep=''),width = 7, height = 7)
# go_enrichment_plot(hg38_da_orth_genes_GO,'DA_peaks')
# dev.off()

# # pdf(paste0(plot_dir,'chimp_DA_orth_GO.pdf',sep=''),width = 7, height = 7)
# # go_enrichment_plot(chimp_da_orth_genes_GO,'DA_peaks')
# # dev.off()

# ## and for non ort genes near DA peaks
# chimp_da_non_orth_genes_GO <- GO_enrichment(chimp_da_peaks[[2]]$feature, "org.Pt.eg.db") 
# hg38_da_non_orth_genes_GO <- GO_enrichment(human_da_peaks[[2]]$feature, "org.Hs.eg.db")

# pdf(paste0(plot_dir,'human_DA_non_orth_GO.pdf',sep=''),width = 7, height = 7)
# go_enrichment_plot(hg38_da_non_orth_genes_GO,'human')
# dev.off()

# pdf(paste0(plot_dir,'chimp_DA_non_orth_GO.pdf',sep=''),width = 7, height = 7)
# go_enrichment_plot(chimp_da_non_orth_genes_GO,'chimp')
# dev.off()




