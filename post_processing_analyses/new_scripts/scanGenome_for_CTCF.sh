## use this script to find the peaks containing CTCF binding sites
wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses'
motif_file="$wd/data/homer_ctcf_ctcfl_pwm.motif"
genome_wide_ctcf="$wd/output/homer/homer_output/ctcf/genome_wide_ctcf"
target_seq_ctcf="$wd/output/homer/homer_output/ctcf/target_sequences"

## create all dirs
if [ ! -d "$ctcf" ]; then
  mkdir -p "$genome_wide_ctcf" && mkdir -p "$target_seq_ctcf";
fi

## loop through array of genome assemblies
declare -a genomes=('hg38' )

for g in "${genomes[@]}"; do
    scanMotifGenomeWide.pl "$motif_file" "$g" -bed -keepAll -mask > "$genome_wide_ctcf/${g}_ctcf.bed"
done
 
## Merge overlapping ranges
cd "$genome_wide_ctcf"
module load bedtools/2.27.1 

for g in "${genomes[@]}"; do
  bedtools merge -i ${g}_ctcf.bed > ${g}_tmp_ctcf.bed
done

echo -e "seqnames\tstart\tend" | cat - ${g}_tmp_ctcf.bed > ${g}_ctcf_merged.bed

## get the DNA seq for all hg38 CTCF peaks and merge the 2 files (i.e., bed + fasta)
Rscript get_dna_seq_peaks.R output/homer/homer_output/ctcf/genome_wide_ctcf/hg38_ctcf_merged.bed output/sequences/hg38_all_ctcf.fa

cd ${wd}/output/sequences/
awk 'BEGIN{RS=">"}{gsub("\n","\t",$0); print ""$0}' hg38_all_ctcf.fa > tmp.bed
tail -n +2  tmp.bed >  tmp_out.bed
echo -e "peakID\tsequence" | cat - tmp_out.bed > hg38_all_ctcf_seq.bed
rm tmp* && rm hg38_all_ctcf.fa
mv hg38_all_ctcf_seq.bed $genome_wide_ctcf

cd $genome_wide_ctcf

paste --delimiters='\t' hg38_ctcf_merged.bed hg38_all_ctcf_seq.bed > hg38_all_ctcf_final.bed

rm *_ctcf.bed && rm hg38_ctcf_merged.bed && rm hg38_all_ctcf_seq.bed

## Find CTCFs motif only in the target fasta files
fasta_dir="$wd/output/sequences"

for fa in "$fasta_dir"/*.fa; do
    basename_fa="$(echo $(basename $fa)| cut -f 1 -d '.')"
    findMotifs.pl "$fa" fasta  "$target_seq_ctcf" -find "$motif_file" > "$target_seq_ctcf/${basename_fa}_discovered_ctcf.txt" -norevopp
done

## Finally remove tmp files if still there
cd "$target_seq_ctcf"
rm *.tmp 
rm motifFindingParameters.txt
