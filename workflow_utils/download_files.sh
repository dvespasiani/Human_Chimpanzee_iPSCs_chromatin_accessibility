## script used to download SRA files (if they dont exist) and encode blacklist regions
## this is done for humans and chimp separately

function check_fastq_files () {
  dir=$1
  expected_no_fastqc=$2
  
  count_files=`ls -1 $dir/*.fastq.gz 2>/dev/null | wc -l`

if [ $count_files == expected_no_fastqc ] 
then
  echo "All files are here"
elif [ $count_files != expected_no_fastqc ]
then
  echo "Downloading files"
  `fastq-dump --split-files --skip-technical --gzip $(<$dir/SraAccList.txt) --outdir $dir `
fi
  
}

echo "Cheking presence human files..."
check_fastq_files "hg38/data/samples" 12

echo "Cheking presence chimp files..."
check_fastq_files "panTro5/data/samples" 10
 
#ENCODE blacklist
echo "Getting hg38 ENCODE blacklisted regions"
encode_dir="data/ENCODE_blacklisted"

if [ ! -d $encode_dir ]; then
  mkdir -p $encode_dir
else
  echo "Directory exists, downloading files"
fi

cd $encode_dir/
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip ./hg38.blacklist.bed.gz

## liftOver
echo "Liftover to create panTro5 blacklisted regions"

liftOver hg38.blacklist.bed  \
../LiftOver_chains/hg38ToPanTro5.over.chain  \
panTro5_blacklist_v2.bed   panTro5_unmapped.bed

sed '/^#/d' panTro5_unmapped.bed  > panTro5_unmapped.bed

echo 'Now you are all set' 
