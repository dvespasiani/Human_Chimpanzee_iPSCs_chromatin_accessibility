## load modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load meme/5.1.1-python-3.7.4


wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses"
input_dir="$wd/output/sequences"
output_dir="$wd/output/fimo_meme"

hocomoco="/data/projects/punim0586/dvespasiani/Annotation_and_other_files/MEME_files/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"


for f in "$input_dir"/* ; do
    fimo_output="$output_dir/$(echo $(basename $f)| cut -f 1 -d '.')/"
    echo "Running fimo on $(basename $f) using $(basename $hocomoco) motifs"
    fimo --bfile --motif-- --thresh 0.01 --o "$fimo_output" "$hocomoco" "$f"
done
