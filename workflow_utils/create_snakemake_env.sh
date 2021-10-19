## load all modules necessary for snakemake when running outside conda env
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load r/4.0.0  
module load snakemake/5.18.1
module load pybedtools/0.8.1-python-3.7.4
module load samtools/1.9
module load samstat/1.5.1
module load openssl/1.1.1f
module load bowtie2/2.3.5.1 
module load ucsc/21072020
module load trimmomatic/0.39-java-11.0.2
# module load --ignore-cache fastqc/0.11.9-java-11
module load macs2/2.2.7.1-python-3.7.4
module load deeptools/3.3.1-python-3.7.4
module load picard/2.6.0-java-11.0.2
module load star/2.7.3a
module load subread/2.0.0