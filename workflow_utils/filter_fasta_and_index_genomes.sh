## script to set up genomes 

module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load samtools/1.9
module load bowtie2/2.3.5.1 
module load ucsc/21072020

homedir='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/'
ucsc_links=('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz','https://hgdownload.soe.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz')

## genomes
genomes=("hg38" "panTro5")


## download fasta files from UCSC
for g in "${genomes[@]}"; do
    for f in ${homedir}${g}'/data/genome/'*.fa.gz; do
    if grep -q ${g} "$f"; then
    echo "$(echo $(basename $f)| cut -f 1 -d '/') found"
    else 
    echo  "$(echo $(basename $f)| cut -f 1 -d '/') NOT found, download it"
fi
done
done


## filter fasta files retaining only standard chromosomes

gunzip panTro5.fa.gz 
grep '>' panTro5.fa | grep -v '_' | sed 's/^>//' | xargs samtools faidx panTro5.fa > panTro5_filtered.fa
rm panTro5.fa
rm panTro5.fa.fai
mv panTro5_filtered.fa panTro5.fa 
gzip panTro5.fa 

gunzip hg38.fa.gz 
grep '>' hg38.fa | grep -v '_' | sed 's/^>//' | xargs samtools faidx hg38.fa > hg38_filtered.fa
rm hg38.fa
rm hg38.fa.fai
mv hg38_filtered.fa hg38.fa 
gzip hg38.fa


## index genome
bowtie2-build hg38.fa.gz ../Index_genome/hg38


## get chromosome sizes
faToTwoBit hg38.fa.gz hg38.2bit
twoBitInfo hg38.2bit stdout | sort -k2rn > hg38.chrom.sizes
