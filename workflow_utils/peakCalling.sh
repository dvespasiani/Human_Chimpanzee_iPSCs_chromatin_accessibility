for f in output/Post_alignment/Files/*tn5_shifted_sorted.bam; do
file_basename="$(echo $(basename $f)| cut -f 1 -d '_')"
macs2 callpeak --format BAMPE --treatment $f \
    --keep-dup all --outdir output/PeakCalling/Files/ \
    --name ${file_basename}_macs2_default --shift 100 --nomodel -B --SPMR \
    --extsize 200 --pvalue 0.01 --call-summits -g 2792339170 
done

for p in  output/PeakCalling/Files/*_macs2_default_peaks.narrowPeak; do
file_basename="$(echo $(basename $p)| cut -f 1 -d '_')"
sort -k 8gr,8gr $p | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | \
head -n 300000 | gzip -nc > output/PeakCalling/Files/${file_basename}_macs2_default_peaks_sorted.narrowPeak.gz
done


