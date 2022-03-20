
## use this script to find the peaks containing CTCF binding sites
base_dir='/data/projects/punim0595/dvespasiani'
working_dir="$base_dir/Human_Chimpanzee_iPSCs_chromatin_accessibility/post_processing_analyses"
homer_dir="$base_dir/homer/motifs"

input_fa="$working_dir/output/sequences/all_peaks_seq.fa"
output_file="$working_dir/output/homer/homer_output/"


# other pluripotency motifs: Rex1, Dax1, Tcl1, Esrrb, Sall4, Klf2, Klf4, Klf5, Stat3 and Tcf3
## declare array of all motifs you want to find
declare -a motifs=("oct4.motif" "sox2.motif" "nanog.motif" "ctcf.motif")

for m in "${motifs[@]}"; do
    basename_motif="$(echo $(basename $m)| cut -f 1 -d '.')"
    
    outdir="$output_file/$basename_motif"
    if [ ! -d "$outdir" ]; then
      mkdir -p "$outdir";
    fi
    
    findMotifs.pl $input_fa fasta "$outdir" -find "$homer_dir/$m" > "$outdir/${basename_motif}_discovered.txt" -norevopp -fdr
    rm $outdir/*.tmp 
    rm $outdir/motifFindingParameters.txt
done
