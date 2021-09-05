## load modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load python/3.7.4 

## actual script
wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses"
input_dir="$wd/output/homer/homer_input"
output_dir="$wd/output/homer/homer_output"

for f in "$input_dir"/*{specific,common}* ; do

homer_output="$output_dir/$(echo $(basename $f)| cut -f 1 -d '.')/"
echo "Running homer on: $(basename $f)"

    if [[ "$f" == *"pantro"* ]];then
        findMotifsGenome.pl "$f" panTro5 "$homer_output" -size given -preparse -p 10
    else
        findMotifsGenome.pl "$f" hg38 "$homer_output" -size given -preparse -p 10 
    fi

done

# for d in "$output_dir"/*; do
# cd  "$d/homerResults/"
#     for mf in *.motif ; do
#         head -n 1  "$mf" >> denovo_file.txt
#         awk -v OFS='\t' '{$1=$1}' denovo_file.txt
#         awk -F'\t' -v OFS='\t' '{ $7 = "'$mf'"}1'  < denovo_file.txt  >> new_denovo_file.txt  ## adding new field with value = filename 
#         rm "$d/homerResults/"denovo_file.txt
#     done
#     ## add colnames at the end
#     ## check out field names for .motif file here: http://homer.ucsd.edu/homer/motif/index.html
#     ## field names: 1) consensus_sequence 2) motif_name 3)log_odds 4)log_pval 5)0 6)motif_occurrences
#     echo -e "consensus_sequence\tmotif_name\tlog_odds\tlog_pval\t0\tmotif_occurrences\tfile_name" | cat - new_denovo_file.txt  > homer_denovo_motifs.txt
#     rm "$d/homerResults/"new_denovo_file.txt
# done
