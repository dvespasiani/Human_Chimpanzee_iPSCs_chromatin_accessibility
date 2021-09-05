
cd /data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/

## activate conda env
# conda activate atac
sinteractive \
 --account=punim0586  \
 --ntasks=1 \
 --cpus-per-task=1 \
 --mem=20000 \
 --time=12:00:00 \
 --partition=mig
 
# list of modules
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load python/3.7.4 
module load r/4.0.0  
# module load meme/5.1.1-python-3.7.4
# module load bedtools/2.27.1 ## untill conda env is not created
# module load sra-toolkit/2.10.5-centos_linux64  #for downloading files, otherwise unnecessary
# module load ucsc/21072020 ## for converting wig to bigwig



