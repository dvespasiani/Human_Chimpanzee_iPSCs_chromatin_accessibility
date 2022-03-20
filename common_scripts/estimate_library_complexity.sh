## script from encode https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
## input filtered bams (dup removed)
module add bedtools/2.29.2

basedir='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/panTro5'
wd="${basedir}/output/Alignment"
indir="${wd}/Files"
outdir="${wd}/qc"

tmp_file="${outdir}/tmp_file.qc"
tmp_filename="${outdir}/tmp_filename.txt"
final_outfile="${outdir}/alignment_qc_pbc.txt"

for f in "${indir}"/*_nochrM_encodefiltered_fixmate_rmorphanread_nodup.bam ; do
    echo "$(echo $(basename $f)| cut -f 1 -d '_')" > ${tmp_filename}

    bedtools bamtobed -i ${f} | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
    END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > results.txt

    paste -d "\t" ${tmp_filename} results.txt >> ${tmp_file}
done 

echo -e "sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2" | cat - ${tmp_file} > ${final_outfile}

rm ${tmp_file} && rm results.txt && rm ${tmp_filename}
