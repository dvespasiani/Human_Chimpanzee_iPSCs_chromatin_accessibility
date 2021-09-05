## Original author IGR, 
## I've simply adjusted it to suit my repo organisation and sometimes cleaned

### THIS SCRIPT IS GENERALLY IDENTICAL TO THAT USED TO ANALYSE THE RNA-SEQ DATA FROM OUR IPSC PAPER IN 2015, BUT THE ENSEMBL BUILD HAS BEEN UPDATED TO HG89. 

### For analysis of DE in iPSCs from human and chimpanzee origin
### Combines what used to be main_analysis_preprocessing.R and main_analysis.R into a single, cleaner script, ommited steps and more detailed descriptions of choices can be found there. 
## Seriously. Read the other one to rationalise some choices and look at some explanations!
### IGR 08.18.14

### Last update: 18.11.01
### Removed samples H18489 and C40210 in keeping with the rest of the analyses.

### load libraries
library(data.table)
library(magrittr)
library(dplyr)
library(ggplots)
library(plyr)
library(RColorBrewer)
library(edgeR) 
library(biomaRt)


options(width=150)
setwd('/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/rna_seq/')

input_dir = './output/PostAlignment/'
plotdir = './plots/'
output_dir = './de_output/'

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

#######################################################
### 1. READ IN GENES AND GENERATE TABLES AND META FILES 
#######################################################

chimp.counts <- fread(paste('./pantro5/',input_dir,"ipsc_chimp_counts.txt",sep=''), header=T)
human.counts <- fread(paste('./hg38/',input_dir,"ipsc_human_counts.txt",sep=''), header=T)

ipsc.counts <- merge(chimp.counts, human.counts,by="GeneID", all=T)
names(ipsc.counts)[1] <- "EnsemblID"
dim(ipsc.counts)
# [1] 42172    13

exon.lengths = fread('orthoexon_files/metaOrthoExon_hg38_panTro5_ensembl86_0.97_ortho_FALSE_pc.txt',header=T)

gene.lengths <- ddply(exon.lengths, "gene", function(x) colSums(x[,c(6,11)]))
gene.lengths$diff <- abs(gene.lengths$widthHs - gene.lengths$widthpanTro5)

#remove unwanted individiuals (already removed).
# ipsc.counts <- ipsc.counts[,colnames(ipsc.counts) != "C8887D" & colnames(ipsc.counts) != "Hutt40"]

#drop the annoying Y chromosome genes
ychr.drop <- c("ENSG00000183878", "ENSG00000198692", "ENSG00000234652", "ENSG00000243576", "ENSG00000252472")
gene.lengths <- gene.lengths[!gene.lengths$gene %in% ychr.drop,]

ipsc.counts <- ipsc.counts[!ipsc.counts$EnsemblID %in% ychr.drop,] 
rownames(ipsc.counts) <- ipsc.counts$EnsemblID ### this may seem dumb, but it makes downstream comparisons a lot easier - certainly much more so than pulling things out of ipsc.genes.dge$genes....
dim(ipsc.counts)
# [1] 42168    13 # What are these genes?

##########################################################################################
### 1.5 WEEDING OUT ALL THE UNDESIRABLE RIBOSOMAL PROTEIN GENES, WHICH MESS UP EVERYTHING. 
##########################################################################################
# listMarts(archive = TRUE)
ensembl.mart <- useMart(
    host='https://dec2016.archive.ensembl.org',
    biomart='ENSEMBL_MART_ENSEMBL', 
    dataset='hsapiens_gene_ensembl'
)

ensembl.biotype <- getBM(
    attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'), 
    filters = 'ensembl_gene_id', 
    values = gene.lengths$gene,
    mart = ensembl.mart
)

ensembl.go <- getBM(
    attributes = c('ensembl_gene_id', 'go_id'), 
    filters = 'ensembl_gene_id', 
    values = gene.lengths$gene,
    mart = ensembl.mart
)

ensembl.biotype$slim <- ensembl.biotype$gene_biotype
ensembl.biotype[grepl("pseudogene", ensembl.biotype$gene_biotype),]$slim <- "pseudogene"

# check the ribosomal genes
ribosomal <- ensembl.go[ensembl.go$go_id %in% "GO:0005840",]
dim(ribosomal)

pseudogene <- ensembl.biotype[ensembl.biotype$slim %in% "pseudogene",]
dim(pseudogene)

genes.to.discard <- unique(c(ribosomal$ensembl_gene_id, pseudogene$ensembl_gene_id))
length(genes.to.discard)
# [1] 7703

    # # try to do this only once, and generate a new gene list which can then be used to filter the genes downstream:
    # ensembl.mart.75 <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
    # ensembl.biotype <- getBM(attributes = c('ensembl_gene_id', 'gene_biotype', 'external_gene_name'), filters = 'ensembl_gene_id', values = gene.lengths$gene, mart = ensembl.mart.75)
    # ensembl.go <- getBM(attributes = c('ensembl_gene_id', 'go_id'), filters = 'ensembl_gene_id', values = gene.lengths$gene, mart = ensembl.mart.75)
    
    # ensembl.biotype$slim <- ensembl.biotype$gene_biotype
    # ensembl.biotype[grepl("pseudogene", ensembl.biotype$gene_biotype),]$slim <- "pseudogene"
    
    # # check the ribosomal genes
    # ribosomal <- ensembl.go[ensembl.go$go_id %in% "GO:0005840",]
    # dim(ribosomal)

    # pseudogene <- ensembl.biotype[ensembl.biotype$slim %in% "pseudogene",]
    # dim(pseudogene)
    
    # genes.to.discard <- unique(c(ribosomal$ensembl_gene_id, pseudogene$ensembl_gene_id))
    # length(genes.to.discard)
    # # [1] 7178
    
    # write.table(genes.to.discard, file="ribosomal_and_pseudo_genes_to_discard_ensembl89.txt")  # This file is correct
 
# genes.to.discard <- read.table(file="ribosomal_and_pseudo_genes_to_discard_ensembl89.txt")

ipsc.counts <- ipsc.counts[!ipsc.counts$EnsemblID %in% genes.to.discard]
dim(ipsc.counts)
# [1] 34465    15

#prepare meta information and plot parameters:
ipsc.counts.names <- colnames(ipsc.counts[,-1])
rpkm.col <- c(rep("darkorchid4", 6), rep("black", 6))
samples.meta <- data.frame(ipsc.counts.names, rpkm.col)
names(samples.meta) <- c("line", "col")
samples.meta$cex <- 1.5

    #alternative colours etc:
    samples.meta$species <- c(rep("chimp", 6), rep("human", 6))
    samples.meta$pch.species <- c(rep(16,6), rep(15,6))
    samples.meta$pch.method <- c(rep(16,5), 15, 14, 15, 15, 17, 17, 14)
    samples.meta$sex <- c('F','M','F','M','M','F','M','M','F','F','M','F') ## i've added this
    samples.meta$col.sex <- ifelse(samples.meta$sex=='M', "yellow3", "plum4")
    samples.meta$col.method <- ifelse(samples.meta$pch.method == 16, "plum4", ifelse(samples.meta$pch.method == 15, "yellow3", ifelse(samples.meta$pch.method == 14, "lightpink", "firebrick")))
    samples.meta$pch.method.strict <- ifelse((samples.meta$pch.method == 16 | samples.meta$pch.method == 14), 16, 15)    
    samples.meta$sex <- ifelse(samples.meta$col.sex == "plum4", "m", "f")
	# and this is for EDAseq:
	rownames(samples.meta) <- samples.meta$line

# Now drop all the samples you don't want, because it's a real pain to drop them above since he code is kind of suspect:
ipsc.counts <- ipsc.counts[,colnames(ipsc.counts) != "C40210A" & colnames(ipsc.counts) != "H18489"]
samples.meta <- samples.meta[rownames(samples.meta) != "C40210A" & rownames(samples.meta) != "H18489",]

save(samples.meta, file="samples_meta_expression_final.Rda")

#########################################
### 2. CALCULATE NORMALISED LIBRARY SIZES
#########################################

# reads into edgeR, calculate TMM and then CPM, writing out the intermediate steps:
ipsc.genes.dge <- DGEList(counts=as.matrix(ipsc.counts[,2:13]), genes=ipsc.counts$EnsemblID)
ipsc.genes.dge <- calcNormFactors(ipsc.genes.dge)
save(ipsc.genes.dge, file="ipsc_genes_dge_ipsc_final_no_ribo_no_ribo.Rda")
write.table(ipsc.genes.dge$samples, file="TMM.normFactors_ipsc_genes_ipsc_final_no_ribo.txt", sep="\t", quote=F, row.names=T)

## some barplots about mapping stats
## "Raw" library sizes:
pdf(paste(plotdir,"mapped_reads_raw_ipsc_final_no_ribo.pdf",sep=''))
mp <- barplot(sort(colSums(ipsc.genes.dge$counts)), ylab="Number of reads mapped to orthologous exons", xlab="", col="darkgrey", xaxt="n")
text(mp, -200000, srt = 45, adj = 1, labels = names(sort(colSums(ipsc.genes.dge$counts))), xpd = TRUE, cex=0.8)
dev.off()

# normalised library sizes
pdf(paste(plotdir,"mapped_reads_normalised_ipsc_final_no_ribo.pdf",sep=''))
mp <- barplot(sort(ipsc.genes.dge$samples$lib.size * ipsc.genes.dge$samples$norm.factor), ylab="Normalized library sizes", xlab="", xaxt="n", col="darkgrey")
text(mp, -200000, srt = 45, adj = 1, labels = row.names(ipsc.genes.dge$samples[order(ipsc.genes.dge$samples$lib.size * ipsc.genes.dge$samples$norm.factor), ]), xpd = TRUE, cex=0.8)
dev.off()

## number of genes expressed
pdf(paste(plotdir,"number_of_genes_expressed_ipsc_final_no_ribo.pdf",sep=''))
some.counts <- apply(ipsc.genes.dge$counts, 2, function(x) { sum(x > 0) })
mp <- barplot(sort(some.counts), ylab="Number of genes with at least 1 read", xlab="", xaxt="n", col="darkgrey")
text(mp, -500, srt = 45, adj = 1, labels = names(sort(some.counts)), xpd = TRUE, cex=0.8)
dev.off()

## Perform rarefaction curves for number of expressed genes vs. proportion of pool mRNA
## As in Ramskold D, Wang ET, Burge CB, Sandberg R. 2009. An abundance of ubiquitously expressed genes revealed by tissue transcriptome sequence data. PLoS Comput Biol 5:e1000598.
## This gives an idea of the complexity of transcriptome in different tissues

## using ortho Exons genes aggregates
pdf(paste(plotdir,"plots/rarefaction_curves_ipsc_final_no_ribo.pdf",sep=''))
plot(1:length(ipsc.genes.dge$counts[,1]), cumsum(sort(ipsc.genes.dge$counts[,1], decreasing=T)/sum(ipsc.genes.dge$counts[,1])), log="x", type="n", xlab="Number of genes", ylab="Fraction of reads pool", ylim=c(0,1)) ## initialize the plot area
for (sample in colnames(ipsc.genes.dge)){
  lines(1:length(ipsc.genes.dge$counts[,sample]), cumsum(sort(ipsc.genes.dge$counts[,sample], decreasing=T)/sum(ipsc.genes.dge$counts[,sample])), col=as.character(samples.meta[samples.meta$line %in% sample,]$col))
}
dev.off()

######################################
### 3. DO NOT NORMALISE FOR GC CONTENT
######################################

## Calculate and export normalized CPMs
cpm.norm <- cpm(ipsc.genes.dge) ## values are not logged (contrary to voom output)
# write.table(cpm.norm, file="new.CPM.TMM.norm_ipsc_final_no_ribo.txt", sep="\t", quote=F, row.names=T)

## Voom requires a design matrix as input
## To make contrasts easier to formulate, we rename factors species and tissue in a single factor
design <- model.matrix(~ 0 + samples.meta$species)
colnames(design) <- c("chimp", "human")

### this design FORCES THE INTERCEPT TO 0 (doing model.matrix ~ etc would NOT), which means we have to use the constrast function down below - need to read more about that, however. 

###################################################################################
### 4. RETAIN ONLY GENES WITH LOG2 CPM > 1 IN AT LEAST 4 INDIVIDUALS IN ONE SPECIES
###################################################################################

par(mfrow=c(1,1))

### Calculate log2 CPM
### Non-normalised data
cpm.norm <- cpm(ipsc.genes.dge, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25) ### I may come back and try this with the Irizarry prior of 0.5, because it makes more intuitive sense to me... but this one is simply -2, so really, no difference
pdf(paste(plotdir,"cpm.density_ipsc_final_no_ribo.pdf",sep=''))
plotDensities(cpm.norm, group=samples.meta$species) 
abline(v=1)
dev.off()

#filter on 4 or more observations of one read per species, not throughout. 
ipsc.genes.dge.filtered <- ipsc.genes.dge[rowSums(cpm.norm[,1:6] > 1) >= 3 | rowSums(cpm.norm[,7:12] > 1) >= 3 , ] 
dim(ipsc.genes.dge.filtered)
# [1] 12674    12
pdf(paste(plotdir,"cpm.filtered.density_ipsc_final_no_ribo.pdf",sep=''))
plotDensities(ipsc.genes.dge.filtered, group=samples.meta$species)
dev.off()

# Recalculate TMM
ipsc.genes.dge.filtered <- calcNormFactors(ipsc.genes.dge.filtered) ## recalculate norm factors
save(ipsc.genes.dge.filtered, file="ipsc.genes.dge.TMM.filtered_final_no_ribo.Rda")


#########################################################################
### 5. LOESS-NORMALISE CPM WITHOUT RANDOM EFFECT USING VOOM - NOT VOOMMOD 
#########################################################################

## Voom on filtered nonGC normalized data, with cyclic loess normalization, without random effect
cpm.loess.norm.norandom.voom <- voom(ipsc.genes.dge.filtered, design, normalize.method="cyclicloess", plot=F) 
pdf(paste(plotdir,"voom.loess.norm.norandom.density_ipsc_final_no_ribo.pdf",sep=''))
plotDensities(cpm.loess.norm.norandom.voom, group=samples.meta$species) # doesn't look much different from the GC-normalised data, which makes sense, given the very high correlation coefficients above... 
dev.off()
# save(cpm.loess.norm.norandom.voom, file="voom.loess.norm.norandom_ipsc_final_no_ribo.Rda")


### For the sake of completeness, also generate the quantile normalised files, even though there's no use for them downstream. 
## Voom on filtered nonGC normalized data, with quantile normalization, without random effect
cpm.quantile.norm.norandom.voom <- voom(ipsc.genes.dge.filtered, design, normalize.method="quantile", plot=F) 
pdf(paste(plotdir,"voom.quantile.norm.norandom.density_ipsc_final_no_ribo.pdf",sep=''))
plotDensities(cpm.quantile.norm.norandom.voom, group=samples.meta$species) #And again, another good example of quantile normalisation working like it is supposed to. 
dev.off()
# save(cpm.quantile.norm.norandom.voom, file="voom.quantile.norm.norandom_ipsc_final_no_ribo.Rda")


########################################################################
### 6. CALCULATE RPKM BY SPECIES AND TEST FOR THE CONSEQUENCES OF EDASEQ
########################################################################

## we could fit the model here, but we transform the CPMs to RPKMs to be sure to remove any difference of read counts due to difference in orthologous exon length between species
## It was really at the RPKM level that things became clearly skewed, so let's take loess and quantile (both GC) forward and test that

## no GC normalized, cyclic loess normalized data
rpkm.loess.norm.norandom.voom <- cpm.loess.norm.norandom.voom
rpkm.loess.norm.norandom.voom$E[,1:6] <- rpkm.loess.norm.norandom.voom$E[,1:6] - log2(gene.lengths[gene.lengths$gene %in% rpkm.loess.norm.norandom.voom$genes[,1],3]/1000)
rpkm.loess.norm.norandom.voom$E[,7:12] <- rpkm.loess.norm.norandom.voom$E[,7:12] - log2(gene.lengths[gene.lengths$gene %in% rpkm.loess.norm.norandom.voom$genes[,1],2]/1000)

pdf(paste(plotdir,"rpkm.voom.loess.norm.density_ipsc_final_no_ribo.pdf",sep=''))
plotDensities(rpkm.loess.norm.norandom.voom, group=samples.meta$species) # not as good as the GC-normalised version, for some reason!
dev.off()
# save(rpkm.loess.norm.norandom.voom, file="voom.RPKM.loess.norm_ipsc_final_no_ribo.Rda")


##############################################################
### 7. PRINCIPAL COMPONENT ANALYSES OF THE NORMALISED DATA SET
##############################################################

# Non GC-corrected data:
ipsc.norm.rpkm.pca <- prcomp(t(rpkm.loess.norm.norandom.voom$E), scale=T, center=T)
pca.var <- ipsc.norm.rpkm.pca$sdev^2/sum(ipsc.norm.rpkm.pca$sdev^2)

pdf(paste(plotdir,"pc1_vs_pc2_loess_norandom_voom_ipsc_final_no_ribo.pdf",sep=''))
plot(
    ipsc.norm.rpkm.pca$x[,1], 
    ipsc.norm.rpkm.pca$x[,2], 
    pch=samples.meta$pch.species, 
    col=as.character(samples.meta$col.method), 
    cex=1.5,
    xlab=paste("PC1 (", round(pca.var[1]*100, digits=2), "% of variance)", sep=""), 
    ylab=paste("PC2 (", round(pca.var[2]*100, digits=2), "% of variance)", sep=""), 
    main=paste("Loess GC norm RPKM at gene level (n=", dim(rpkm.loess.norm.norandom.voom$E)[1], ")", sep="")
)
legend(
    x="topleft", 
    col=c("plum4", "yellow3", "firebrick4", "lightpink"), 
    pch=c(16,15,15,15), 
    pt.cex=1.5, 
    legend=c("Chimpanzee (fib)", "Yoruba (LCL)", "Caucasian (LCL)", "Caucasian (fib)"),
    bty="n", cex=0.8
)
dev.off()    

### 07.30.14 testing the correlation between the pcs and a couple of variables:

# Generate additional samples.meta columns:
samples.meta$provenance <- ifelse(samples.meta$col.method == "plum4", "fib", ifelse(samples.meta$col.method == "lightpink", "fib", "lcl"))
samples.meta$species.provenance <- paste(samples.meta$species, samples.meta$provenance, sep=".")
samples.meta$ethnicity <- c(rep("chimp", 6), "yri", "caucasian", "yri", rep("caucasian", 3))
samples.meta$population <- paste(samples.meta$provenance, samples.meta$ethnicity, sep=".")

# Define the testing function
ipsc_final.assoc <- function(pca.data){
    all.pcs <- data.frame()
    for (i in 1:12){
        species <- anova(lm(pca.data$x[,i] ~ samples.meta$species))$Pr[1]
        provenance <- anova(lm(pca.data$x[,i] ~ as.factor(samples.meta$provenance)))$Pr[1]
        species.provenance <- anova(lm(pca.data$x[,i] ~ as.factor(samples.meta$species.provenance)))$Pr[1]
        ethnicity <- anova(lm(pca.data$x[,i] ~ as.factor(samples.meta$ethnicity)))$Pr[1]
        single.pc <- c(i, species, provenance, species.provenance, ethnicity)
        all.pcs <- rbind(all.pcs, single.pc)
    }
    names(all.pcs) <- c("PC", "species", "provenance", "species.provenance", "ethnicity")
    return(all.pcs)
}

ipsc.assoc <- ipsc_final.assoc(ipsc.norm.rpkm.pca)
ipsc.assoc

#    PC      species provenance species.provenance    ethnicity
# 1   1 6.088398e-12 0.06749507       8.221944e-10 4.842926e-11
# 2   2 9.878898e-01 0.13690313       4.223080e-01 5.948352e-01
# 3   3 9.622052e-01 0.71185172       9.794129e-01 7.495273e-01
# 4   4 8.293710e-01 0.79333455       1.705073e-01 4.411588e-01
# 5   5 9.734219e-01 0.16125381       4.345390e-01 8.021754e-01
# 6   6 9.395908e-01 0.51515351       5.098223e-01 9.200204e-01
# 7   7 9.974407e-01 0.91831901       8.529425e-01 8.796126e-01
# 8   8 9.262865e-01 0.11541429       8.093314e-02 3.798712e-02
# 9   9 9.429005e-01 0.92598070       9.317351e-01 9.170798e-01
# 10 10 9.895060e-01 0.75262700       9.743682e-01 8.418393e-01
# 11 11 9.611264e-01 0.84499518       9.785548e-01 9.984041e-01
# 12 12 5.473281e-04 0.45222916       1.245736e-03 3.137560e-03

# The above table suggests that the humans are indeed segregating by ethnicity and not cell type, but the presence of chimps confounds things, since they have a third ethnic label. So I have two choices:

# 1. don't touch the PCA, but subset to the columns I want. 
    ipsc_final.assoc.cols <- function(pca.data){
        all.pcs <- data.frame()
        for (i in 1:12){
            provenance <- anova(lm(pca.data$x[7:12,i] ~ as.factor(samples.meta$provenance[7:12])))$Pr[1]
            ethnicity <- anova(lm(pca.data$x[7:12,i] ~ as.factor(samples.meta$ethnicity[7:12])))$Pr[1]
            population <- anova(lm(pca.data$x[7:12,i] ~ as.factor(samples.meta$population[7:12])))$Pr[1]
            single.pc <- c(i, provenance, ethnicity, population)
            all.pcs <- rbind(all.pcs, single.pc)
        }
        names(all.pcs) <- c("PC", "provenance", "ethnicity", "population")
        return(all.pcs)
    }

ipsc.human.cols.assoc <- ipsc_final.assoc.cols(ipsc.norm.rpkm.pca)
ipsc.human.cols.assoc
#    PC provenance  ethnicity population
# 1   1 0.65778183 0.13241483 0.36750958
# 2   2 0.81891163 0.01385217 0.00824398
# 3   3 0.05206393 0.54831256 0.19534615
# 4   4 0.95474192 0.47928217 0.76658353
# 5   5 0.27374333 0.60909826 0.60262564
# 6   6 0.43349791 0.90304844 0.64683071
# 7   7 0.99391497 0.50883370 0.77461012
# 8   8 0.74101913 0.83906557 0.95388686
# 9   9 0.37115166 0.51019035 0.68935227
# 10 10 0.85890848 0.97284461 0.98499792
# 11 11 0.30083319 0.83758821 0.59929763
# 12 12 0.95372543 0.24354244 0.40917443

# 2. Deleted after conversations with Rachel, John and Yoav. I agree with all of them. 


##################################################################################################
### 8. PREPARE FOR DE TESTING, CHECK DISPERSIONS AND OTHER INDICATORS OF SUCCESSFUL NORMALISATION. 
##################################################################################################

# Define the model matrices:
# Note that the second one is problematic because it has more contrasts than voom had, but I can't really work out how to nest something right now. 
design <- model.matrix (~ 0 + samples.meta$species)
colnames(design) <- c("chimp", "human")

design.sex.interact <- model.matrix(~ 0 + samples.meta$species:samples.meta$sex)
colnames(design.sex.interact) <- c("chimp.female", "human.female", "chimp.male", "human.male")

## No GC normalisation + loess norm, without random effect
fit.loess.norm.norandom <- lmFit(rpkm.loess.norm.norandom.voom, design)
fit.loess.norm.norandom <- eBayes(fit.loess.norm.norandom)

fit.loess.norm.norandom.sex.interact <- lmFit(rpkm.loess.norm.norandom.voom, design.sex.interact)
fit.loess.norm.norandom.sex.interact <- eBayes(fit.loess.norm.norandom.sex.interact)

# Non-GC normalised MA plots:
## If 'MA' is an 'MArrayLM' object, then the plot is a fitted model MA-plot in which the estimated coefficient is on the y-axis and the average A-value is on the x-axis.
pdf(paste(plotdir,"model_MA_plots_design_loess.norm.norandom.voom_ipsc_final_no_ribo.pdf",sep=''))
    par(mfrow=c(1,2))
    limma::plotMA(fit.loess.norm.norandom, coef=1, xlab="average coefficient", ylab="estimated coefficient") ## for a bunch of exons with medium intensity (-2<cpm<2), in chimp heart they are either highly expressed or not expressed. Highly expressed genes (on average) seem to be less variable in individual samples
    abline(a=0, b=1,col=pal[1])
    limma::plotMA(fit.loess.norm.norandom, coef=2, xlab="average coefficient", ylab="estimated coefficient")
    abline(a=0, b=1,col=pal[1])
dev.off()

pdf(paste(plotdir,"model_MA_plots_design.sex.interact_loess.norm.norandom.voom_ipsc_final_no_ribo.pdf",sep=''))
    par(mfrow=c(2,2))
    limma::plotMA(fit.loess.norm.norandom.sex.interact, coef=1, xlab="average coefficient", ylab="estimated coefficient") ## for a bunch of exons with medium intensity (-2<cpm<2), in chimp heart they are either highly expressed or not expressed. Highly expressed genes (on average) seem to be less variable in individual samples
    abline(a=0, b=1,col=pal[1])
    limma::plotMA(fit.loess.norm.norandom.sex.interact, coef=2, xlab="average coefficient", ylab="estimated coefficient")
    limma::plotMA(fit.loess.norm.norandom.sex.interact, coef=3, xlab="average coefficient", ylab="estimated coefficient")
    limma::plotMA(fit.loess.norm.norandom.sex.interact, coef=4, xlab="average coefficient", ylab="estimated coefficient")
dev.off()

## Sigma vs A plot. After a linear model is fitted, this checks constancy of the variance with respect to intensity level.
pdf(paste(plotdir,"model_SA_plots_loess.norm.norandom.voom_ipsc_final_no_ribo.pdf",sep=''))
    par(mfrow=c(1,2))
    plotSA(fit.loess.norm.norandom, main="species")
    plotSA(fit.loess.norm.norandom.sex.interact, main="species.sex")
dev.off()

## There is a non flat trend across all datasets (quite similar): big problem?
## The flatness seems the best for microarrays, but there is often high variability for low intensity probes
## If there is a strong trend, eBayes(trend=T) should be used, but this is not valid with RNA-seq/voom.
## see: https://stat.ethz.ch/pipermail/bioconductor/2013-July/053835.html
## and: http://permalink.gmane.org/gmane.science.biology.informatics.conductor/54146
## Final thing not clear: is the red curve the trend accounted for, or the remaining trend?

## lmFit output object (see ?gls.series, called by lmFit)
## - coefficients: numeric matrix containing the estimated coefficients for each linear model. Same number of rows as 'M', same number of columns as 'design'.
## - stdev.unscaled: numeric matrix conformal with 'coef' containing the unscaled standard deviations for the coefficient estimators. The standard errors are given by 'stdev.unscaled * sigma'.
## - sigma: numeric vector containing the residual standard deviation for each gene.
## - df.residual: numeric vector giving the degrees of freedom corresponding to 'sigma'
## - correlation: inter-duplicate or inter-block correlation
## - qr: QR decomposition of the generalized linear squares problem, i.e., the decomposition of 'design' standardized by the Choleski-root of the correlation matrix defined by 'correlation'
## See also: help("MArrayLM-class")

## Boxplot of the residuals (difference between the true value and fitted value) for the 5 datasets tested
pdf(paste(plotdir,"model_residual_boxplots_loess.norm.norandom.voom_ipsc_final_no_ribo.pdf",sep=''))
boxplot(log2(fit.loess.norm.norandom$sigma), log2(fit.loess.norm.norandom.sex.interact$sigma), names=c("species", "species:sex"))
dev.off()
### The distribution of the residuals across all three looks very very similar... and very, very normal. Yay!

### This is clearly the bit where I am most lacking in awareness, and should read the limma manual from beginning to end, to understand what these plots show and what I want them to show. 

## To compare no norm, and cyclic loess, no random: find which scheme produces the largest standard deviation on average (either measured by the median value or the mean value) 
## fit.GC.loess.norm$sigma has the lowest median of residuals (very slight trend)


##################################################
### 9. IDENTIFY DE GENES ACCORDING TO BOTH SCHEMES
##################################################

ipsc.cm <- makeContrasts(
                    ChimpvsHuman = (chimp-human),
                    levels=design)

ipsc.cm.sex.interact <- makeContrasts(
                    ChimpvsHuman = ((chimp.female+chimp.male)-(human.female+human.male)),
                    FemalevsMale = ((chimp.female+human.female)-(chimp.male+human.male)),
                    levels=design.sex.interact)
 
### Non GC data
## fit the contrasts
fit2.loess.norm.norandom <- contrasts.fit(fit.loess.norm.norandom, ipsc.cm)
fit2.loess.norm.norandom <- eBayes(fit2.loess.norm.norandom)

fit2.loess.norm.norandom.sex.interact <- contrasts.fit(fit.loess.norm.norandom.sex.interact, ipsc.cm.sex.interact)
fit2.loess.norm.norandom.sex.interact <- eBayes(fit2.loess.norm.norandom.sex.interact)

## topTables
topSpecies.loess.norm.norandom <- topTable(fit2.loess.norm.norandom, coef=1, adjust="BH", number=Inf, sort.by="p")
tests.loess.norm.norandom <- decideTests(fit2.loess.norm.norandom, adjust.method="BH", method="separate", p.value=0.01)
table(tests.loess.norm.norandom)
# tests.loess.norm.norandom
#   -1    0    1 
# 2098 8450 2126   

topSpecies.loess.norm.norandom.sex.interact <- topTable(fit2.loess.norm.norandom.sex.interact, coef=1, adjust="BH", number=Inf, sort.by="p")
topSpecies.loess.norm.norandom.sex.fm <- topTable(fit2.loess.norm.norandom.sex.interact, coef=2, adjust="BH", number=Inf, sort.by="p")

### Write out all output:

# Non GC data:
write.table(topSpecies.loess.norm.norandom, paste(output_dir,"topSpecies.loess.norm.norandom_ipsc_final_no_ribo.out",sep=''), row.names=T, col.names=T, quote=F) ## use this file, containing list of DE genes between chimp and humans
write.table(topSpecies.loess.norm.norandom.sex.interact, paste(output_dir,"topSpecies.loess.norm.norandom.sex.interact_ipsc_final_no_ribo.out",sep=''), row.names=T, col.names=T, quote=F)
write.table(topSpecies.loess.norm.norandom.sex.fm, paste(output_dir,"topSpecies.loess.norm.norandom.sex.fm_ipsc_final_no_ribo.out",sep=''), row.names=T, col.names=T, quote=F)

##########################################################################################################################################
### 10. PLOT SOME RESULTS - THE PROTEIN CODING THRESHOLDING MAKES THE HEATMAPS A LOT NICER, BECAUSE THERE'S NO SECONDARY 0 EXPRESSION PEAK
##########################################################################################################################################

#Define the function
check_results <- function(count_table, top, n) {
  cols <- colorRampPalette(pal)(n)
  plot(1,1, xlim=c(1, length(count_table[1,])), ylim=c(min(count_table[rownames(count_table) %in% top$genes[1:n],]$E), max(count_table[rownames(count_table) %in% top$genes[1:n],]$E)), type="n", xaxt="n", ylab="Normalized expression", xlab="")
  text(1:length(count_table[1,]), min(count_table[rownames(count_table) %in% top$genes[1:n],]$E)-0.07*(max(count_table[rownames(count_table) %in% top$genes[1:n],]$E)-min(count_table[rownames(count_table) %in% top$genes[1:n],]$E)), srt = 45, adj = 1,labels = colnames(count_table), xpd = TRUE)
  axis(side=1, at=c(1:length(count_table[1,])), labels=rep("", 12))
  abline(v=1:length(count_table[1,]), lty=3, col="lightgray")
 
  if (n > length(top[,1])){
    n <- length(top[,1])
  }
  for (i in 1:n){
    lines(1:length(count_table[1,]), count_table[rownames(count_table) %in% top$genes[i],]$E, col=cols[i])
   }
}

# Define the number of genes to plot:
gene.plots <- c(10,50,100,200,500,1000)
colors <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)) ## We want red to be higher expression and blue lower

# Non GC norm:
for (i in 1:6){
pdf(file=paste("check_results_topSpecies.loess.norm.norandom_", gene.plots[i], "_ipsc_final_no_ribo.pdf", sep=""))
check_results(rpkm.loess.norm.norandom.voom, topSpecies.loess.norm.norandom, gene.plots[i])
dev.off()

pdf(file=paste("heatmap_topSpecies.loess.norm.norandom_", gene.plots[i], "_ipsc_final_no_ribo.pdf", sep=""))
heatmap.2(rpkm.loess.norm.norandom.voom$E[rownames(topSpecies.loess.norm.norandom[1:gene.plots[i],]),], col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=samples.meta$species, labRow=NA, ColSideColors=pal[as.integer(as.factor(samples.meta$species))], cexCol = 1.5)
dev.off()
}

pdf(file="heatmap_topSpecies.loess.norm.norandom_all_genes_ipsc_final_no_ribo.pdf")
heatmap.2( rpkm.loess.norm.norandom.voom$E[rownames(topSpecies.loess.norm.norandom[1:dim(topSpecies.loess.norm.norandom)[1],]),], scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=samples.meta$species, labRow=NA, ColSideColors=pal[as.integer(as.factor(samples.meta$species))], cexCol = 1.5)
dev.off()

pdf(file="heatmap_topSpecies.loess.norm.norandom_all_DE_ipsc_final_no_ribo.pdf")
heatmap.2( rpkm.loess.norm.norandom.voom$E[topSpecies.loess.norm.norandom$adj.P.Val < 0.01,], scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=samples.meta$species, labRow=NA, ColSideColors=pal[as.integer(as.factor(samples.meta$species))], cexCol = 1.5)
dev.off()
