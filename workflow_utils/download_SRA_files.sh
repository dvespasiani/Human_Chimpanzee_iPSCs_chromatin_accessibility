## script used to download SRA files (if they dont exist) and encode blacklist regions
## this is done for humans and chimp separately

helpFunction()
{
   echo ""
   echo "Usage: $0 -o outdir -f file"
   echo -e "\t-o Directory where you want to save the downloaded files"
   echo -e "\t-f Relative path to the accession txt file containing all the SRA files to download"
   echo -e "\t base dir: $wd"
   exit 1 # Exit script after printing help
}

while getopts "o:f:" flag; do
    case "${flag}" in
        o) outdir="$OPTARG";;
        f) file="$OPTARG";;
        ?) helpFunction ;; # Print helpFunction 
    esac
done

## Begin script
wd="/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility"

# load modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load sra-toolkit/2.10.5-centos_linux64  

## run the fastq-dump command like this otherwise it doesnt work
command="fastq-dump --split-files --skip-technical --gzip "$(<$wd/$file)" --outdir "$wd/$outdir""
$command



