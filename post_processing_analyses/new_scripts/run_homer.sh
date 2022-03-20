# ## load modules
# source /usr/local/module/spartan_new.sh
# module load web_proxy
# module load gcc/8.3.0 openmpi/3.1.4
# module load python/3.7.4 

# ## actual script
# wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses"
# input_dir="$wd/output/homer/homer_input"
# output_dir="$wd/output/homer/homer_output"

# for f in "$input_dir"/* ; do

# homer_output="$output_dir/$(echo $(basename $f)| cut -f 1 -d '_')/"

# echo "Running homer on: $(basename $f)"

#     if [[ "$f" == "pantro"* ]];then
#         findMotifsGenome.pl "$f" panTro5 "$homer_output" -size given -preparse -p 10
#     else
#         findMotifsGenome.pl "$f" hg38 "$homer_output" -size given -preparse -p 10 
#     fi

# done

## load modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load python/3.7.4 

## actual script
wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses"
input_dir="$wd/output/sequences"
output_dir="$wd/output/homer/homer_output"

findMotifs.pl ${input_dir}/all_peaks_seq.fa \
human $output_dir -fasta ${input_dir}/random_regions_seq.fa \
-mask

