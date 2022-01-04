## this code is what I used to predicted the affinity of CTCFs for the underlying binding sites

## unfortunately I cannot yet run deepbind on the cluster so i have to do it locally (no big deal is v quick)
## another thing is that I modified the example.ids file in order to contain only the  D00328.003 # CTCF (SELEX) model

## loop extracts sequence and peakID field names from bed file 

wd='/Users/davidevespasiani'
deepbind_dir="${wd}/deepbind"
input_dir="${wd}/Desktop/Projects/human_chimp_atac/fimo_ctcf" 
output_dir="${wd}/Desktop/Projects/human_chimp_atac/deepbind_predictions"

if [ ! -d "$output_dir" ]; then
   mkdir -p "$output_dir";
fi


cd "$deepbind_dir"

for f in "$input_dir"/*ctcf_peaks.txt ; do
  file_basename=$(echo $(basename $f)| cut -f 1 -d '_')
  awk '{print $8}' $f | tail -n +2 > ${file_basename}_sequences.txt 
  awk '{print $4}' $f >  ${file_basename}_peakIDs.txt
  deepbind_output="$output_dir/${file_basename}"
  ./deepbind  --echo  example.ids < ${file_basename}_sequences.txt > "${deepbind_output}_tmp.txt"
  paste -d '\t' ${file_basename}_peakIDs.txt "${deepbind_output}_tmp.txt" > "${deepbind_output}_ctcf_affinity_predictions.txt"
  rm "${deepbind_output}_tmp.txt"
  rm ${file_basename}_sequences.txt
  rm ${file_basename}_peakIDs.txt
done




