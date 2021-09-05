library(data.table)
library(magrittr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# # install devtools if not previously installed
# if (!is.element('devtools', installed.packages()[,"Package"])) install.packages('devtools')

# # # install dependencies
# devtools::install_github("demuellae/muLogR")
# devtools::install_github("demuellae/muRtools")
# devtools::install_github("demuellae/muReportR")
# # install ChrAccR
# devtools::install_github("GreenleafLab/ChrAccR", dependencies=TRUE)

# # hg38 annotation package
# install.packages("https://muellerf.s3.amazonaws.com/data/ChrAccR/data/annotation/ChrAccRAnnotationHg38_0.0.1.tar.gz")

suppressWarnings(library(ChrAccR))

## annotation file
## for this analysis it must contain:
## sample identifiers 
## bamFilename
## full path of bamFilename

# sampleAnnotFn <- file.path("tcells", "samples.tsv")
# bamDir <- file.path("tcells", "bam")
# sampleAnnot <- fread(sampleAnnotFn, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# # Download an ATAC peak set from a Pan-cancer dataset [doi:10.1126/science.aav1898]
# cancerPeaks <- read.table("https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec", sep="\t", header=TRUE, comment.char="", quote="", stringsAsFactors=FALSE)
# cancerPeaks <- GenomicRanges::GRanges(cancerPeaks[,"seqnames"], IRanges::IRanges(start=cancerPeaks[,"start"], end=cancerPeaks[,"end"], names=cancerPeaks[,"name"]))
# cancerPeaks <- sort(cancerPeaks) # sort
# # a list of custom region sets
# regionSetList <- list(
#   tiling500bp = muRtools::getTilingRegions("hg38", width=500L, onlyMainChrs=TRUE),
#   cancerPeaks = cancerPeaks
# )
## import data from bam file
# dsa_fromBam <- DsATAC.bam(sampleAnnot, "bamFilenameFull", "hg38", regionSets=regionSetList, sampleIdCol="sampleId")

## path to Tn5 shifted files
peakQC_dir='./output/PeakCalling/QC'
sampleAnnot_file <- list.files(path=peakQC_dir,pattern = "\\SampleAnnotation.txt$",full.names = T)
sampleAnnot <- read.table(sampleAnnot_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)


## Get MACS2 pooled filtered peaks
filtered_peaks<- fread("./output/PeakCalling/Files/human_macs2_default_peaks_filtered.narrowPeak.gz", sep="\t", header=F, stringsAsFactors=F,select=c(1:3))
colnames(filtered_peaks)<-c('seqnames','start','end')

filtered_peaks_gr <- makeGRangesFromDataFrame(filtered_peaks)%>%sort()

## make regionSetList
regionSetList <- list(
  tiling500bp = muRtools::getTilingRegions("hg38", width=500L, onlyMainChrs=TRUE),
  filtered_peaks = filtered_peaks_gr
)

## import data from bam file
dsa<- DsATAC.bam(sampleAnnot, bamFiles="bamfilesFull", genome="hg38", regionSets=regionSetList)

## make GRanges object containing TSS coordinates
tssGr <- muRtools::getAnnotGrl.gencode("gencode.v27")[["gene"]]
tssGr <- tssGr[elementMetadata(tssGr)[,"gene_type"]=="protein_coding"]
tssGr <- promoters(tssGr, upstream=0, downstream=1)
tssGr

## compute TSS enrichment per sampleId
samples_id=getSamples(dsa)

tsse_all=list()
for (i in samples_id){
    tsse_all[[i]]=getTssEnrichment(dsa, i, tssGr)
}
names(tsse_all)=samples_id

## calculate and plot the TSS enrichment score,
## i.e. number of insertions at the TSS/ number of insertion in background regions
tsse_enrich = do.call("rbind", lapply(tsse_all, "[[", 2))
tsse_enrich=as.data.frame(tsse_enrich)
tsse_enrich=mutate(tsse_enrich,'Sample'=rownames(tsse_enrich))
colnames(tsse_enrich)[1]<-'TSS_Enrichment'

## write file
write.table(tsse_enrich,"./output/PeakCalling/QC/TSS_enrichment.txt",sep='\t',col.names=T,quote=F,row.names=F)