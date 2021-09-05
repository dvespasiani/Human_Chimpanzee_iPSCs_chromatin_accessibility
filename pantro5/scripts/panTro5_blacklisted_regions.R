## Use this script to convert genoimic coordinates of the ENCODE blacklisted regions from hg38 to panTro5
library(data.table);library(dplyr);library(magrittr)
library(rtracklayer);library(GenomicRanges)
library(R.utils)

hg38_blacklisted=fread('../data/ENCODE_blacklisted/hg38_blacklist_v2.bed',sep='\t',header=F,
col.names = c('seqnames','start','end','feature'))
liftover_chain=import.chain('../data/LiftOver_chains/hg38ToPanTro5.over.chain')

liftover=function(file,chain){
  gr=copy(file)%>%makeGRangesFromDataFrame(keep.extra.columns = T,ignore.strand=T)

  gr=liftOver(gr,chain) 
  gr=unlist(gr) %>% as.data.table()
  return(gr)
}
panTro5_blacklisted=liftover(hg38_blacklisted,liftover_chain)
panTro5_blacklisted=panTro5_blacklisted[,c('width','strand','feature'):=NULL]

## remove those where start==end (this might be structural variants in humans)
panTro5_blacklisted=panTro5_blacklisted[!start==end]%>%unique()%>%setorderv(c('seqnames','start','end'),1)

## if start/end is the same across different rows use min/max start/end
panTro5_blacklisted=panTro5_blacklisted[, .SD[which.min(start)], by=.(seqnames,end)]
panTro5_blacklisted=panTro5_blacklisted[, .SD[which.max(end)], by=.(seqnames,start)]

# panTro5_gz=gzfile("../data/ENCODE_blacklisted/panTro5_blacklist_hg38lifted.bed.gz", "w")

write.table(test,"../data/ENCODE_blacklisted/panTro5_blacklist_hg38lifted.bed",sep='\t',col.names = F,quote=F,row.names = F)

