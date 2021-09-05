## script to test if there is any enrichment for specific tf/cluster of tf
## between da peaks and non da peaks between the two species
##  and between peaks near orth genes and not


library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggthemes);library(ggplot2);library(ggpubr)
library(viridis);library(viridisLite)
library(ggrepel)

# hg38_fimo = snakemake@input[[]]
# pantro5_fimo = snakemake@input[[]]
# motif_cluster_file = snakemake@input[[]]
# peaks = snakemake@input[[]]

options(width=150)

orth_peaks_dir = './Orthologous_peaks/'
da_peaks_dir = './output/DA/peaks/'
plot_dir = '/data/gpfs/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/plots/'
merging_keys =  c("seqnames", "start",'end')

read_fimo = function(x){
    df= fread(x,sep='\t',header=T)
    df=df[,motif_alt_id:=NULL]%>%
    na.omit()%>%
    setnames(old=c('start','stop','p-value','q-value'),new=c('motif_start','motif_end','p_value','q_value'))
    df=df[q_value<=0.01
    ][
      ,c(merging_keys,'species') := tstrsplit(sequence_name, "_", fixed=TRUE)
      ][
        ,start:=as.numeric(start)
        ][
          ,end:=as.numeric(end)][
            ,c('sequence_name','score','p_value','matched_sequence','strand'):=NULL
          ][ 
            ,.SD[which.min(q_value)], by=.(seqnames,start,end,motif_start,motif_end)] ## this keeps only the moost signif motif per peak (or >1 but that have diff motif ranges)
    return(df)
    }

hg38_fimo_results = read_fimo('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/fimo_meme/human/fimo.tsv')
pantro5_fimo_results = read_fimo('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/output/fimo_meme/chimp/fimo.tsv')

all_fimo_results = rbind(hg38_fimo_results,pantro5_fimo_results)

## now get peaks with all the other info
read_peaks = function(x,organism){
    df = fread(paste(da_peaks_dir,x,sep=''),
    sep='\t',header=T,select=c('seqnames','start','end','width','FDR','direction','has_orth_gene'))
    df = df[,peak_type:=ifelse(FDR<=0.01,'da','non_da')][,species:=organism]
    return(df)
}
all_human_tmm_peaks = read_peaks('all_tmm_human_peaks.txt','human')
all_chimp_hg38coord_tmm_peaks = read_peaks('all_tmm_chimp_peaks_hg38coord.txt','chimp')

all_tmm_peaks = rbind(all_human_tmm_peaks,all_chimp_hg38coord_tmm_peaks)

motif_cluster = fread('/data/projects/punim0586/dvespasiani/Files/Annotation_and_other_files/Motifs_clusters/motifs_clusters',sep=' ',header = T)%>%setnames(old='Motif',new='motif_id')

## add cluster info to tfs
all_fimo_results = all_fimo_results[motif_cluster,on='motif_id',nomatch=0]
all_fimo_results = all_tmm_peaks[all_fimo_results,on=c(merging_keys,'species'),nomatch = 0]
all_fimo_results = all_fimo_results[
  ,c('seqnames','start','end','peak_type','species','Cluster','Name','motif_id','has_orth_gene') ## add orth column
  ]%>%unique()
 
##--------------------------------------------------------
## first split by orth non orth peaks
## and test if there is an enrichment for any given motif
## in da vs non-da peaks 
##--------------------------------------------------------
all_fimo_results_split = all_fimo_results%>%split(as.factor(all_fimo_results$has_orth_gene))
names(all_fimo_results_split)=c('no_orth','yes_orth')

make_matrix = function(x,organism){
  motif_matrix = copy(x)
  motif_matrix=motif_matrix[species==organism][,motif_id:=NULL]%>%unique()

  motif_matrix=motif_matrix[
        ,number_motifs_peaktype:=.N,by=.(Name,peak_type)
      ][
        ,tot_numb_peaks_per_type:=.N,by=.(peak_type)
      ][
        ,number_motifs_notin_peaktype:=tot_numb_peaks_per_type-number_motifs_peaktype
      ][
        ,tot_numb_peaks_per_type:=NULL
      ][,c(merging_keys,'has_orth_gene'):=NULL]%>%unique()

  motif_matrix_da = copy(motif_matrix)[peak_type=='da'][,peak_type:=NULL]%>%setnames(old='number_motifs_peaktype',new='number_motifs_da_peaks')
  motif_matrix_nonda = copy(motif_matrix)[peak_type=='non_da'][,peak_type:=NULL]%>%setnames(old='number_motifs_peaktype',new='number_motifs_nonda_peaks')
  
  motif_matrix_final = motif_matrix_da[motif_matrix_nonda,on=c('Name','Cluster','species')]%>%na.omit() ## remove those with na values in da peaks (i.e. no entries for da peaks)
  return(motif_matrix_final)
}

fisher_pvalues=function(a){
  pvals=copy(a)[,-c(1:3)]
  pvals=apply(pvals, 1, 
              function(x) {
                tbl <- matrix(as.numeric(x), ncol=2, byrow=T)
                fisher.test(tbl)$p.value
              })
  
  pvalues_table=data.table(pval=pvals)
  pvalues_table=pvalues_table[
    ,adj_p:=p.adjust(pval,method = 'fdr')][
      ,log10_p_adjust:=-log10(adj_p)
      ][
        ,significant_score:=ifelse(`adj_p`<=0.0001,'****',
                                   ifelse(`adj_p`>0.0001 &`adj_p`<=0.001,'***',
                                          ifelse(`adj_p`>0.001 & `adj_p`<=0.01,'**',
                                                 ifelse(`adj_p`>0.01 & `adj_p`<=0.05,'*',' '))))
        ]
  
  return(pvalues_table)
}

motif_enrichment = function(x,organism){
  motif_matrix  = make_matrix(x,organism)
  pvalues = fisher_pvalues(motif_matrix)
  enrichment_matrix = cbind(motif_matrix,pvalues)
  enrichment_matrix=enrichment_matrix[
      ,ratio_da_nonda_motifs_in_peaktype := number_motifs_da_peaks/number_motifs_nonda_peaks
      ][
      ,ratio_da_nonda_motifs_notin_peaktype := number_motifs_notin_peaktype/i.number_motifs_notin_peaktype
      ][
        ,odds_ratio:=ratio_da_nonda_motifs_in_peaktype/ratio_da_nonda_motifs_notin_peaktype
        ][
          ,log10_number_motifs_da_peaks:=log10(number_motifs_da_peaks)
        ][
            ,c('Name','Cluster','odds_ratio','pval','adj_p',
              'significant_score','log10_p_adjust','log10_number_motifs_da_peaks')
            ]
  return(enrichment_matrix)

}

human_orth_motif_enrichmnet = motif_enrichment(all_fimo_results_split[[2]],'human')
chimp_orth_motif_enrichmnet= motif_enrichment(all_fimo_results_split[[2]],'chimp')

human_nonorth_motif_enrichmnet = motif_enrichment(all_fimo_results_split[[1]],'human')
chimp_nonorth_motif_enrichmnet= motif_enrichment(all_fimo_results_split[[1]],'chimp')


## plot enrichment as if only u knew how do it
tf_enrich_plot=function(x){
  df=copy(x)[,'Log10 total number DA peaks per cluster':=log10_number_motifs_da_peaks]
  range_x_axis = range(log2(df$odds_ratio))
  gradient=scale_colour_viridis(aes(`Log10 total number aSNPs per cluster`),option="viridis",discrete = F)
  text=ifelse(!df$significant_score%in%' ',df$Name,'')

  ggplot(df,aes(x=log2(odds_ratio),y=log10_p_adjust,label = text,col=log10_number_motifs_da_peaks))+
    geom_point(size=2)+
    geom_vline(xintercept=0, linetype="dashed", color = "black",size=0.5)+
    geom_hline(yintercept=1.5, linetype="dashed", color = "black",size=0.5)+
    geom_text_repel(size = 5,color='black',
                    box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.5, "lines")
    )+
    gradient+
    xlab('\n log2 odds ratio \n')+
    ylab('-Log10 (P-adj)')+
    theme(
      strip.background.y = element_blank(),
      strip.background.x =element_blank(),
      panel.spacing=unit(1, "lines"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background =element_rect(fill = 'white', colour = 'black',size=1),
      legend.position = "bottom",
      legend.key = element_rect(fill = "white", colour = "black"),
      axis.line = element_blank()
      )
}


pdf(paste(plot_dir,'human_hocomoco_enrich_da_vs_nonda_peaks_orth_genes.pdf',sep=''),width = 6,height=6)
tf_enrich_plot(human_orth_motif_enrichmnet)
dev.off()

pdf(paste(plot_dir,'chimp_hocomoco_enrich_da_vs_nonda_peaks_orth_genes.pdf',sep=''),width = 6,height=6)
tf_enrich_plot(chimp_orth_motif_enrichmnet)
dev.off()

pdf(paste(plot_dir,'human_hocomoco_enrich_da_vs_nonda_peaks_NON_orth_genes.pdf',sep=''),width = 6,height=6)
tf_enrich_plot(human_nonorth_motif_enrichmnet)
dev.off()

pdf(paste(plot_dir,'chimp_hocomoco_enrich_da_vs_nonda_peaks_NON_orth_genes.pdf',sep=''),width = 6,height=6)
tf_enrich_plot(chimp_nonorth_motif_enrichmnet)
dev.off()

##----------------------------------------------------------------------------
## given there is enrichment for some tfbs in da over non-da peaks
## now do it again but by chromstate 
## See if here there are difference between orth/non orth peaks 
## otherwise this is another evidence in support of a shared chromatin state
## between the species at the iPSCs state
##----------------------------------------------------------------------------

## read chrom state
chrom_state_dir = '/data/projects/punim0586/dvespasiani/Files/Annotation_and_other_files/Roadmap_data/ChromHMM_15states_cell_lines_combined/Continuous_states/'
ipsc_chromstate = fread(paste(chrom_state_dir,'IPSC.gz',sep=''),sep=' ',header=T)

setkeyv(ipsc_chromstate,merging_keys)
setkeyv(all_fimo_results,merging_keys)

tf_peak_annotation = foverlaps(all_fimo_results,ipsc_chromstate,type='within')[
  ,c('start','end','cell_type','motif_id'):=NULL]%>%
  setnames(old=c('i.start','i.end'),new=c('start','end'))%>%
  na.omit()%>% ## removes peaks in sex chromosomes[,c('start','end'):=NULL]%>%
  unique()

tf_peak_annotation_split = tf_peak_annotation%>%split(as.factor(tf_peak_annotation$has_orth_gene))


zscore_chromstate = function(x,organism){
  matrix = copy(x)[species==organism][
          ,number_motifs_peaktype_chromstate:=.N,by=.(Name,peak_type,chrom_state)
        ][,c(merging_keys,'has_orth_gene','species'):=NULL]%>%unique()
      
  matrix_da = copy(matrix)[peak_type=='da'][,peak_type:=NULL]%>%setnames(old='number_motifs_peaktype_chromstate',new='number_motifs_da_peaks_chromstate')
  matrix_nonda = copy(matrix)[peak_type=='non_da'][,peak_type:=NULL]%>%setnames(old='number_motifs_peaktype_chromstate',new='number_motifs_nonda_peaks_chromstate')
    
  matrix_final = matrix_da[matrix_nonda,on=c('Name','Cluster','chrom_state')]%>%na.omit() ## remove those with na values in da peaks (i.e. no entries for da peaks)
  matrix_final = matrix_final[
    ,ratio_numb_daVSnonda_peaks_per_chromstate_motif:=number_motifs_da_peaks_chromstate/number_motifs_nonda_peaks_chromstate
  ][ 
    , chromstate_wise_zscore:=.SD[ , scale(ratio_numb_daVSnonda_peaks_per_chromstate_motif)], by=.(chrom_state)
    ][
      ,c('chrom_state','Name','chromstate_wise_zscore')][,species:=organism]%>%unique()
  return(matrix_final)

}



human_erichm_chromstate_nonorth_peaks = zscore_chromstate(tf_peak_annotation_split[[1]],'human')
chimp_erichm_chromstate_nonorth_peaks = zscore_chromstate(tf_peak_annotation_split[[1]],'chimp')

all_species_enrich_chromstate_nonorth_peaks = rbind(human_erichm_chromstate_nonorth_peaks,chimp_erichm_chromstate_nonorth_peaks)

human_erichm_chromstate_orth_peaks = zscore_chromstate(tf_peak_annotation_split[[2]],'human')
chimp_erichm_chromstate_orth_peaks = zscore_chromstate(tf_peak_annotation_split[[2]],'chimp')

all_species_enrich_chromstate_orth_peaks = rbind(human_erichm_chromstate_orth_peaks,chimp_erichm_chromstate_orth_peaks)


heatmap_matrix =function(x){
  df = copy(x)[
    ,rowname:=paste(species,chrom_state,sep='.')]
  df=df[,c('species','chrom_state'):=NULL]%>%dcast(rowname~Name,value.var='chromstate_wise_zscore')
  rowname=df$rowname
  df = df[,c('rowname'):=NULL]
  df[is.na(df)]=0
  df_matrix = as.matrix(df)
  rownames(df_matrix)=rowname
  return(df_matrix)
}

matrix_enrich_chromstate_orth_peaks = heatmap_matrix(all_species_enrich_chromstate_orth_peaks)
matrix_enrich_chromstate_nonorth_peaks = heatmap_matrix(all_species_enrich_chromstate_nonorth_peaks)

## plot heatmap

## Legends 
species_colors=structure(c('#C99E10','#9B4F0F'),names=c('chimp','human')) 


nihroadmap_colors=as.data.table(ipsc_chromstate)[
  ,col:=plyr::revalue(`chrom_state`,c('1_TssA'='#FF0000','2_TssAFlnk'='#FF6E00','3_TxFlnk'='#32CD32','4_Tx'='#008000',
                                      '5_TxWk'='#006400','6_EnhG'='#C2E105','7_Enh'='#FFFF00',
                                      '8_ZNF/Rpts'='#66CDAA','9_Het'='#8A91D0','10_TssBiv'='#CD5C5C','11_BivFlnk'='#E9967A',
                                      '12_EnhBiv'='#BDB76B','13_ReprPC'='#3A3838','14_ReprPCWk'='#808080',
                                      '15_Quies'='#DCDCDC'))
  ][
    ,c('chrom_state','col')
    ] %>% unique()


chromHMM_colors=nihroadmap_colors$col
names(chromHMM_colors)=nihroadmap_colors$chrom_state

library(ComplexHeatmap)
library(circlize)
library(viridis)

enrichment_heatmap=function(enrichmatrix){
  x = Heatmap(
    enrichmatrix, 
      border = T, 
      col= viridis(1000),
      row_dend_reorder = T,
      row_order=order(as.character(gsub("^.*\\.", "", rownames(enrichmatrix)))),
      show_heatmap_legend = F,
      show_row_names = F,
      row_title =" ",
      show_column_names = T,
      # column_order =order(as.factor(readr::parse_number(gsub("^.*\\.", "",colnames(enrichmatrix))))),
      # column_split = as.factor(readr::parse_number(gsub("^.*\\.", "", colnames(enrichmatrix)))),
      column_title =' ',
      right_annotation = HeatmapAnnotation(
      which='row',

      species=anno_simple(
        sort(gsub("\\..*", "",rownames(enrichmatrix))),
        border = T,
        height=unit(0.5,'cm'),
        col=species_colors),
        
        chrom_state=anno_simple(
          gsub("^.*\\.", "", rownames(enrichmatrix)),
          border=T,
          height = unit(0.5,'cm'),
          col=chromHMM_colors),
          
          show_annotation_name = F)
        )
  return(x)
}


pdf(paste(plot_dir,'all_species_hocomoco_chromstate_zscore_davsnonda_peaks_orth_genes.pdf',sep=''),width = 6,height=6)
enrichment_heatmap(matrix_enrich_chromstate_orth_peaks)
dev.off()

pdf(paste(plot_dir,'all_species_hocomoco_chromstate_zscore_davsnonda_peaks_NON_orth_genes.pdf',sep=''),width = 6,height=6)
enrichment_heatmap(matrix_enrich_chromstate_nonorth_peaks)
dev.off()
