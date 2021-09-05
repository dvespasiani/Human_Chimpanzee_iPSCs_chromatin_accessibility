* change MDS into PCA for fig 1
* convert all annotations to hg38 
* keep go as supp to show that u are using peaks near orth genes.
* use fimo for peak in chrom states and see if there is any enrichment (especially for bivalent peaks)
* check if there is sign enrichment of da peaks over bivalent states (do this for human and chimp separately)
* keep peaks that have orthologous regions in the other 2 genomes for separate analyses
* make table with intial number of peaks and then with the filtered ones

From Irene's paper (on liftOver):
Chain regions used to generate non-overlapping 100 bp  windows. Liftover then used to test whether each of these windows
could be lifted over to a single site in the chimp genome and then back again to original location in human genome.
They allowed for max 20% changes in size during each lifover process, if changes were >  reads discarded.

* calculate uniqueness of each 50 k-mers in both genomes using GEM (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030377)
From Irene's paper (on this k-mer uniqueness):
they removed windows where > 20% of the 50-mers in either species were not uniquely mappable in their respective genomes.

## tips from laura
get the hg38 coordinates of the pantro5 peaks (and vice versa) then keep those that pass the 2 filters u mentioned: 1) a liftover back in which i will retain only those peaks that are included in the original consensus set; 2) those that have the same nearest gene (if any). I will then merge the resulting peaks (in the respective species-coordinates) into a single peak file and use this to count the reads overlapping each region

## Things to discuss with IGR
1. Chimp files are odd, evidences:
    * snakemake pipeline takes ages for chimp whereas for human doenst (e.g. atac-seq QCs) even though commands are the same
    * some chimp peaks are insanely wide whereas human arent (this shouldnt happen as I've filtered for proper paired reads during post-aligment QCs)

2. also some human peaks are longish albeit not as much as the chimp ones (human longest is < 5kb)

I'd keep all peaks <=  1kb for both species and then separately investigate the longer ones (if they are in the consensus peak set)

weqy3. How to perform DA analysis, specifically which peak to keep after liftover
