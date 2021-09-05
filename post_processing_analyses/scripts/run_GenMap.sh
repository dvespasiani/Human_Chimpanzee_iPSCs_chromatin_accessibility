
## load atac conda env as this is where program is installed
conda activate atac

## index genomes
genmap index -FD /data/projects/punim0586/dvespasiani/tmp_hg38_pantro5/ -I /data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/Genome_mappability/indexed_genome

## calculate mappability
genmap map -K 50 -E 2 -I Genome_mappability/indexed_genome/ -O Genome_mappability/mappability/ --bed --exclude-pseudo -T 50

## polish file by removing info on unwanted chromosomes
filename=$(find *bed)
awk '$1!="/random/ || /Un/ || /./"' $filename > cleaned

rm $filename
awk '!/Un/' cleaned > $filename
rm cleaned

## run this on /var/local/tmp is a fast NVMe drive local to the compute node.
mkdir -p /var/local/tmp/davide/indexed_genome
cp /var/local/tmp/davide/indexed_genome/

genmap map -K 50 -E 2 -I /var/local/tmp/davide/indexed_genome/ -O Genome_mappability/mappability/ --bed --exclude-pseudo -T 50 --memory-mapping

rm -Rf /var/local/tmp/davide/indexed_genome

