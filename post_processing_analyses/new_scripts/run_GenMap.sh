
# conda activate atac
# source /usr/local/module/spartan_new.sh
# module load web_proxy
# module load gcc/8.3.0 openmpi/3.1.4
# module load python/3.7.4

# helpFunction(){
#    echo ""
#    echo "Usage: $0 -g genome"
#    echo -e "\t-g genome to index"
#    exit 1 # Exit script after printing help
# }

# while getopts "g:" flag; do
#     case "${flag}" in
#         g) genome="$OPTARG";;
#         ?) helpFunction ;; # Print helpFunction
#     esac
# done

## tmp directories
mkdir /data/gpfs/projects/punim0595/dvespasiani/tmp
export TMPDIR=/data/gpfs/projects/punim0595/dvespasiani/tmp

wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility'
g='hg38'
genome_dir="$wd/$g/data/genome"
index_dir="$wd/$g/genome_mappability/indexed_genome"
mappability_dir="$wd/$g/genome_mappability/mappability"

## create dirs
if [ ! -d "$mappability_dir" ]; then
  mkdir -p "$mappability_dir";
fi

## index genome
genmap index -FD $genome_dir -I  $index_dir

## calculate k-mer mappability
genmap map -K 50 -E 2 -I $index_dir -O $mappability_dir -bg --exclude-pseudo -T 6 -v





# ## load atac conda env as this is where program is installed
# conda activate atac

# declare -a genomes=("hg38" "panTro5")

# fasta_dir='/data/projects/punim0586/dvespasiani/tmp_hg38_pantro5'
# wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility'

# # for i in "${genomes[@]}"; do
# #     fasta_file="${fasta_dir}/${i}.fa"
# #     index_dir="${wd}/${i}/genome_mappability/index"
# #     mappability_dir="${wd}/${i}/genome_mappability/mappability"
# #     kmer=50
    
# #     ## first index genomes
# #     echo "Indexing $fasta_file"
# #     genmap index -F $fasta_file -I $index_dir
    
# #     ## then calculate mappability
# #     echo "Calculating $kmer mappability for $fasta_file"
# #     genmap map -K $kmer -E 2 -I $index_dir -O $mappability_dir --bed --exclude-pseudo -T 50

# #     echo "Polishing files"
# #     cd $mappability_dir
# #     filename=$(find *bed)
# #     awk '$1!="/random/ || /Un/ || /./"' $filename > cleaned
# #     # rm $filename
# #     awk '!/Un/' cleaned > cleaned_${filename}
# #     rm cleaned

# # done
# wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses'
# index_dir="${wd}/Genome_mappability/indexed_genome"
# mappability_dir="${wd}/Genome_mappability/mappability"

# ## index genomes
# genmap index -FD /data/projects/punim0586/dvespasiani/tmp_hg38_pantro5/ -I /data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses/Genome_mappability/indexed_genome

# ## create all mappability dir
# if [ ! -d "$mappability_dir" ]; then
#   mkdir -p "$mappability_dir";
# fi

# ## calculate mappability
# genmap map -K 50 -E 2 -I $index_dir -O $mappability_dir --bed --exclude-pseudo -T 50

# ## polish file by removing info on unwanted chromosomes
# filename=$(find *bed)
# awk '$1!="/random/ || /Un/ || /./"' $filename > cleaned

# rm $filename
# awk '!/Un/' cleaned > $filename
# rm cleaned

# # ## run this on /var/local/tmp is a fast NVMe drive local to the compute node.
# # mkdir -p /var/local/tmp/davide/indexed_genome
# # cp /var/local/tmp/davide/indexed_genome/

# # genmap map -K 50 -E 2 -I /var/local/tmp/davide/indexed_genome/ -O Genome_mappability/mappability/ --bed --exclude-pseudo -T 50 --memory-mapping

# # rm -Rf /var/local/tmp/davide/indexed_genome

