## script used to liftover the consensus peaks coordinates using reciprocal best hits chain files from ucsc
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load ucsc

wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility"
chainsDir="$wd/data/LiftOver_chains"
consensusPeaksDir="$wd/post_processing_analyses/output/files/consensusPeaks"

hg38TopanTro5="$chainsDir/hg38.panTro5.rbest.chain"
panTro5Tohg38="$chainsDir/panTro5.hg38.rbest.chain"

hg38Peaks="$consensusPeaksDir/human_consensus_peaks.bed"
panTro5Peaks="$consensusPeaksDir/chimp_consensus_peaks.bed"

lift_peaks() {

    fbname="$(basename "$1" .bed)"
    chain="$2"

    tmp_file="$consensusPeaksDir/tmp"
    touch $tmp_file

    tail -n +2 $1 > $tmp_file

    echo
    echo "liftovering $1"
    echo
    liftOver -bedPlus=3 -tab $tmp_file $chain "${tmp_file}_mapped.bed" "${tmp_file}_unmapped.bed"

    colnames="seqnames\tstart\tend\tpeakID"

    echo -e $colnames | cat - "${tmp_file}_mapped.bed" > "$consensusPeaksDir/${fbname}_mapped.bed"
    echo -e $colnames | cat - "${tmp_file}_unmapped.bed" | sed '/^#/d' >  "$consensusPeaksDir/${fbname}_unmapped.bed" | sed '/^#/d' 

    echo
    echo "done"
    echo
}


lift_peaks $hg38Peaks $hg38TopanTro5
lift_peaks $panTro5Peaks $panTro5Tohg38

rm $consensusPeaksDir/tmp*
