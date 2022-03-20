## directories 
tads_dir <-  './output/TADs/'
da_dir <- './output/DA/'
plot_dir <- './output/plots/'
outdir <- './output/'
bamDir  <- "/output/Alignment/Files"
genome <- 'hg38'

## vectors
samples_names <- c(
  "C3647",
  "C3649",
  "C3651",
  "C40280",
  "C4955",
  "C8861",
  "H19101",
  "H19114",
  "H20961",
  "H28834",
  "HUtt45",
  "HUtt60"
)

range_keys <- c('seqnames','start','end')
species_names <- c('chimp','human')
peak_type <- c('da','non_da')
chain_path <- '/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/data/LiftOver_chains/'
standard_chr <- paste0("chr", c(1:23,'2A','2B', "X", "Y")) # only use standard chromosomes

chrom_states <- c(
  "1_TssA","2_TssAFlnk","3_TxFlnk",
  "4_Tx","5_TxWk","6_EnhG",
  "7_Enh","8_ZNF/Rpts","9_Het",
  "10_TssBiv","11_BivFlnk","12_EnhBiv",
  "13_ReprPC","14_ReprPCWk","15_Quies"
)

##----------------
## color palette
##----------------
samples_palette <- c(
  "#219ebc","#468189","#f3722c",
  "#f8961e","#f9c74f","#f94144",
  "#90be6d","#f4e9cd","#585123",
  "#43aa8b","#577590","#f58549"
)
names(samples_palette) <- samples_names

species_palette <- c('#ff7d00','#15616d')
names(species_palette) = species_names

da_palette <- c('#E9D8A6','#0A9396')
names(da_palette) <- peak_type

chrom_state_colors <- c(
  '#FF0000','#FF6E00','#32CD32',
  '#008000','#006400','#C2E105',
  '#FFFF00','#66CDAA','#8A91D0',
  '#CD5C5C','#E9967A','#BDB76B',
  '#3A3838','#808080','#DCDCDC'
)
names(chrom_state_colors) <- chrom_states

##---------------------------
## functions
##---------------------------
# export_file <- function(directory,filename){
#     genomes <- c('hg38','panTro5')
#     dirs <- list()
#     for (g in genomes){
#         dirs[[g]] <- paste(directory,g,sep='')
#         dir.create(dirs[[g]], showWarnings = FALSE)
#     }
#     file <- purrr::map2(dirs,genomes,function(x,y){
#         f <- paste(x,'/',y,'_',filename,sep='')
#         return(f)
#     })

#     return(file)
# }

create_dir <- function(base_dir,path){
  dir <- paste(base_dir,path,'/',sep='')
  dir.create(dir,showWarnings=F,recursive = TRUE)
  return(dir)
}

export_file <- function(directory,genome,filename){
    dir <- paste(directory,genome,sep='/')
    dir.create(dir, showWarnings = FALSE)
    file <- paste(dir,'/',genome,'_',filename,sep='')
    return(file)
}

## annotate peaks in chromatin states
annotate_peaks <- function(peaks){
  annotated <- foverlaps(copy(peaks),ipsc_chromstate,type='any')%>%na.omit()
  annotated <- annotated[
    ,overlap:=ifelse(i.start<start,i.end-start,end-i.start),by=.(cell_type)
    ][
        ,.SD[which.max(overlap)], by=.(peakID,cell_type)
        ][
            ,c(range_keys[-1],'overlap'):=NULL
            ]%>%setnames(old=c('i.start','i.end'),new=c(range_keys[-1]))
  return(annotated)
}

## counts peaks in chromatin states
count_peaks_chromstate <- function(peaks){
  counts <- copy(peaks)
  counts<-counts[
      ,numb_peaks_chromstate:=.N,by=.(chrom_state,cell_type)
      ][
        ,numb_peaks:=.N,by=.(cell_type)
        ][
            ,c('cell_type','chrom_state','numb_peaks_chromstate','numb_peaks')
            ]%>%unique()
  return(counts)
}
## calculate OR
calculate_or <- function(peaks_oi,peaks_noi,merging_keys){
  peaksoi_vs_peaksnoi_or <- merge(peaks_oi,peaks_noi,by=c(merging_keys))%>%
  dplyr::select(c(contains('numb'),contains(all_of(merging_keys))))

  fisher_test <- copy(peaksoi_vs_peaksnoi_or)%>%split(by=c(merging_keys))%>%
  lapply(
    function(x){
    x <- x[,c(merging_keys):=NULL]%>%as.numeric()%>%matrix(nrow=2,byrow=T)%>%fisher.test()
    x <- data.table(
        'p'=x$p.value,
        'odds_ratio'=x$estimate,
        'lower_ci'=x$conf.int[[1]],
        'upper_ci'=x$conf.int[[2]]
        )
    }
  )
  or_results <- Map(mutate,fisher_test,elements=names(fisher_test))%>%rbindlist()
  or_results <- p_significance(or_results,or_results$p)
  return(or_results)
}

## plot OR 
plot_or <- function(peaks){
  p <- ggplot(peaks,aes(x=factor(elements,levels=chrom_states),y=log(odds_ratio),fill=elements))+
    geom_violin(trim=T,scale = "width")+
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1),binwidth=0.2)+
    scale_fill_manual(values = chrom_state_colors)+
    geom_hline(yintercept=0,linetype='dashed')+
    xlab('chromatin state') + ylab('log OR')+
    theme_classic()+
    theme(
      legend.position = "bottom",
      axis.text.x =element_blank(),
      axis.ticks.x =element_blank()
    )
  return(p)
}

## empirical permutations
permute_data <- function(metric, peak_type, n=10000){
  score_distribution = c()
  observed_value = diff(by(metric,peak_type,mean))
  for(i in 1:n){
      score_distribution[i]=diff(by(metric, sample(peak_type, length(peak_type), replace = F), mean))
  }
  df_score_distribution = data.table(permuted_scores=score_distribution)
  zscore = (observed_value-mean(score_distribution))/sd(score_distribution)
  pvalue = 2*pnorm(q=abs(zscore), lower.tail=FALSE)
  df_stat =  data.table(zscore = zscore,pval=pvalue)
  return = list(observed_value,df_score_distribution,df_stat)
  names(return) = c('observed_value','permuted_scores','stat_results')
  return(return)
}

## read files 
read_da_results <-  function(file,significance){
    df <- fread(paste(da_dir,file,sep=''),sep='\t',header=T)%>%setnames(
      old=c('human_seqnames','human_start','human_end','human_peakID'),new=c(range_keys,'peakID'))
    setkeyv(df,range_keys)
    return(df)
}

## TADs
read_tads = function(file){
  tad = fread(paste(tads_dir,file,sep=''),sep='\t',header=T)%>%setnames(old='disc_species',new='species')
  return(tad)

}
## iPSC chrom state hg38 (these files contain also info for sex chromosomes)
read_chromstate  = function(dir){
    state = list.files(dir,full.names=T,recursive=F)%>%
    lapply(
        function(x)fread(x,sep='\t',header=F,col.names=c(range_keys,'chrom_state'))
        )
    names(state) = gsub('\\_.*','', list.files(dir,full.names=F,recursive=F))

    state = Map(mutate,state,cell_type=names(state))%>%rbindlist()
    state = state[seqnames!='chrM']
    return(state)

}

##-------------------
## count numb peaks
##-------------------
count_peaks = function(x){
    x=copy(x)%>%as.data.table()
    x=x[,c(..range_keys)]%>%unique()%>%nrow()

    print(paste('total number of peaks:', x,sep= ' '))
}

##---------------------------
## bin peak sizes/distances
##---------------------------
bin_distance = function(x,column){
    x=x[
    ,binned_column:= ifelse(..column < 50, '0-49',
                    ifelse(..column >= 50 & ..column < 151,'50-150',
                    ifelse(..column >= 151 & ..column < 301,'151-300',
                    ifelse(..column >= 301 & ..column < 451,'301-450',
                    ifelse(..column >= 451 & ..column < 601,'451-600',
                    ifelse(..column >= 601 & ..column < 751,'601-750',
                    ifelse(..column >= 751 & ..column < 901,'751-900',
                    ifelse(..column >= 901 & ..column < 1001,'901-1000',
                    ifelse(..column >= 1001 & ..column < 1151,'1001-1150',
                    ifelse(..column >= 1151 & ..column < 1301,'1151-1301',
                    ifelse(..column >= 1301 & ..column < 1501,'1301-1500',
                    ifelse(..column >= 1501 & ..column < 2001,'1501-2000',
                    ifelse(..column >= 2001 & ..column < 3001,'2001-3000',
                    ifelse(..column >= 3001 & ..column < 4001,'3001-4000',
                    ifelse(..column >= 4001 & ..column < 5001,'4001-5000',
                    ifelse(..column >= 5001 & ..column < 6001,'5001-6000','>6000'
                    ))))))))))))))))
                    ]
    return(x)
}

##-------------
## LiftOver
##-------------
liftPeaks <-  function(peaks,chain_file){
    chain <- rtracklayer::import.chain(chain_file)
    peaks_gr = makeGRangesFromDataFrame(peaks,keep.extra.columns=T)
    seqlevelsStyle(peaks_gr) = "UCSC" 
    lifted_coord <- liftOver(peaks_gr, chain)%>%reduce(min.gapwidth=50L)%>%unlist()%>%as.data.table()
    return(lifted_coord)
}

convert_coord <- function(peaks,chain_file){
    chain <- rtracklayer::import.chain(paste(chain_path,chain_file,sep=''))
    original_peak_df <- copy(peaks)[,width:=end-start]
    original_peak_gr <- makeGRangesFromDataFrame(original_peak_df,keep.extra.columns=T)
    seqlevelsStyle(original_peak_gr) = "UCSC" 
    names(original_peak_gr)=original_peak_gr$peakID
    
    lifted_peaks = liftOver(original_peak_gr, chain)

    lifted_peaks  = unlist(lifted_peaks)%>%as.data.table()
    keep = copy(as.data.table(lifted_peaks))[,c('peakID','seqnames')]%>%unique()
    keep = keep[,diff_chr:=.N,by=.(peakID)][diff_chr==1][seqnames %in% standard_chr]

    peaks_to_keep = copy(lifted_peaks)[
        peakID %in% keep$peakID
        ][
            original_peak_df[,c('width','peakID')],on=c('peakID'),nomatch=0
            ]
    merge_peaks = copy(peaks_to_keep)[
        ,start:=min(start),by=.(peakID)
        ][
            ,end:=max(end),by=.(peakID)
            ][
                ,width:=end-start
                ][
                    ,size_change:=round(width/i.width,2)
                    ][
                        size_change>=0.8 & size_change<=1.2
                        ][
                            ,c('width','i.width','strand','size_change'):=NULL
                            ]%>%unique()
    setkeyv(merge_peaks,range_keys)

    return(merge_peaks)
}

##--------------
## pvalues
##--------------
get_pval =function(x){
  df=copy(x)
  test_stat = t.test(df$enrichment,mu= 1)$statistic
  pval = t.test(df$enrichment,mu= 1)$p.value
  degrees_freedom = t.test(df$enrichment,mu= 1)$parameter
  df=df[,test_stat:= test_stat][,p := pval][,df := degrees_freedom]
  return(df)
  }

##
adjust_pvalues=function(x){
  df=copy(x)
  pvals_df=copy(df)
  pvals_df=pvals_df$p
  pvals_df_adjusted = p.adjust(pvals_df,'fdr')%>%as.data.table() %>% setnames('p.adj')
  pvals_df_adjusted = p_significance(pvals_df_adjusted,pvals_df_adjusted$p.adj)
  df_final=cbind(df,pvals_df_adjusted)
  return(df_final)
}

##
p_significance = function(x,p){
  df=copy(x)
  df=df[
    ,p.signif:= ifelse(p<=0.0001,'****',
                       ifelse(p>0.0001 & p <=0.001,'***',
                              ifelse(p>0.001 & p<=0.01,'**',
                                     ifelse(p>0.01 & p<=0.05,'*',' '))))
    ]
    return(df)
}

##
pval_vector=function(df,columns){
  x = copy(df)[,c(..columns)] %>% unique()
  x=x[, c("custom_column", "chrom_state") := tstrsplit(x[[1]], ".", fixed=TRUE)]
  x=x[order(as.factor(readr::parse_number(gsub("^.*\\.", "",x$chrom_state)))),]
  x=split(x,by='chrom_state')
  x=lapply(x,function(y)y=y[order(y$custom_column)][,p.signif:=ifelse(p.signif==' ',' ','*')])%>%rbindlist()
  named_vector=structure(as.character(x$p.signif),names=as.character(x[[1]]))
  return(named_vector)  
}

##-------------------------------
## color legends for heatmaps
##-------------------------------

chromstatus_color=function(df,column){
  x=copy(df)[,c(..column)]
    x=x[
        ,chrom_state:=gsub("^.*\\.","", x[[1]])
        ][
    ,col :=plyr::revalue(`chrom_state`,c('1_TssA'='khaki3','2_TssAFlnk'='khaki3','3_TxFlnk'='khaki3',
                                         '4_Tx'='khaki3','5_TxWk'='khaki3','6_EnhG'='khaki3',
                                         '7_Enh'='khaki3','8_ZNF/Rpts'='khaki3','9_Het'='ivory3',
                                         '10_TssBiv'='ivory3','11_BivFlnk'='ivory3','12_EnhBiv'='ivory3',
                                         '13_ReprPC'='ivory3','14_ReprPCWk'='ivory3','15_Quies'='ivory3'))
    ]%>%unique()
    chromcols = x$col
    names(chromcols)=x$chrom_state

    return(chromcols)
}
## 
nihroadmap_colors = function(df,column){
    x=copy(df)[,c(..column)]
     x=x[
        ,chrom_state:=gsub("^.*\\.","", x[[1]])
        ][
        ,col:=plyr::revalue(
            `chrom_state`,c(
                '1_TssA'='#FF0000',
                '2_TssAFlnk'='#FF6E00',
                '3_TxFlnk'='#32CD32',
                '4_Tx'='#008000',
                '5_TxWk'='#006400',
                '6_EnhG'='#C2E105',
                '7_Enh'='#FFFF00',
                '8_ZNF/Rpts'='#66CDAA',
                '9_Het'='#8A91D0',
                '10_TssBiv'='#CD5C5C',
                '11_BivFlnk'='#E9967A',
                '12_EnhBiv'='#BDB76B',
                '13_ReprPC'='#3A3838',
                '14_ReprPCWk'='#808080',
                '15_Quies'='#DCDCDC'))
                ][,c('chrom_state','col')] %>% unique()
                                    
    chromcols = x$col
    names(chromcols)=x$chrom_state

    return(chromcols)
}


## 
colors_custom_column = function(df,column,custom_colors){
  column_values = copy(df)[
    ,c(..column)
    ]
column_values=column_values[
      ,custom_col:=gsub("\\..*","",column_values[[1]])
      ][
        ,'custom_col']%>%dplyr::pull()%>%unique()
  names(custom_colors) = column_values
  return(custom_colors)
}