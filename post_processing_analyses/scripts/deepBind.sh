## this code is what I used to predicted the affinity of CTCFs for the underlying binding sites

## unfortunately I cannot yet run deepbind on the cluster so i have to do it locally (no big deal is v quick)
## another thing is that I modified the example.ids file in order to contain only the  D00328.003 # CTCF (SELEX) model

deepbind_dir='/Users/dvespasiani/deepbind'
ctcf_seq_dir='/Users/dvespasiani/Desktop/atac_plots/ctcf_peak_sequences'
outdir='/Users/dvespasiani/Desktop/atac_plots/deepbind_predictions'

## create dirs
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir";
fi

cd "$deepbind_dir"

for seq in "$ctcf_seq_dir"/*sequences.txt ; do
deepbind_output="$outdir/$(echo $(basename $seq)| cut -f 1 -d '_')"

    ./deepbind  --echo  example.ids < $seq > "${deepbind_output}_ctcf_affinity_predictions.txt"

done
