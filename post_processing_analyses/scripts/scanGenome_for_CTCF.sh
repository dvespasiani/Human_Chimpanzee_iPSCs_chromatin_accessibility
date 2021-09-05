## use this script to find the peaks containing CTCF binding sites
wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses'
motif_file="$wd/data/homer_ctcf_ctcfl_pwm.motif"
genome_wide_ctcf="$wd/output/homer/homer_output/ctcf/genome_wide_ctcf"
target_seq_ctcf="$wd/output/homer/homer_output/ctcf/target_sequences"

## create all dirs
if [ ! -d "$ctcf" ]; then
  mkdir -p "$genome_wide_ctcf" && mkdir -p "$target_seq_ctcf";
fi

## declare the array first
declare -a genomes=('hg38' 'panTro5')

## then loop through it. remember the [@]!!
for genome in "${genomes[@]}"; do
    scanMotifGenomeWide.pl "$motif_file" "$genome" -bed -keepAll -mask > "$genome_wide_ctcf/${genome}_ctcf.bed"
done
 
## Merge overlapping ranges
cd "$genome_wide_ctcf"
module load bedtools/2.27.1 

bedtools merge -i hg38_ctcf.bed > hg38_ctcf_merged.bed
bedtools merge -i panTro5_ctcf.bed > panTro5_ctcf_merged.bed

rm *_ctcf.bed

## Find CTCFs motif only in the target fasta files
fasta_dir="$wd/output/sequences"

for fa in "$fasta_dir"/*.fa; do
    basename_fa="$(echo $(basename $fa)| cut -f 1 -d '.')"
    findMotifs.pl "$fa" fasta  "$target_seq_ctcf" -find "$motif_file" > "$target_seq_ctcf/${basename_fa}_discovered_ctcf.txt"
done

## Finally remove tmp files if still there
cd "$target_seq_ctcf"
rm *.tmp 
rm motifFindingParameters.txt
