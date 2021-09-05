## This script is used to liftover the bam files between the 2 assemblies 

helpFunction()
{
   echo ""
   echo "Usage: $0 -d input_dir -c chain_file"
   echo -e "\t-d Directory where the bams are stored"
   echo -e "\t-c name of the chain file"
#    echo -e "\t base dir: $wd"
   exit 1 # Exit script after printing help
}

while getopts "d:c:" flag; do
    case "${flag}" in
        d) input_dir="$OPTARG";;
        c) chain_file="$OPTARG";;
        ?) helpFunction ;; # Print helpFunction 
    esac
done


wd='/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility'
chain_dir="$wd/data/LiftOver_chains"

## run command
cd "$input_dir"
shopt -s extglob
# for bam in !(*merged*)_blacklist_removed.bam; do 
#     CrossMap.py bam -a "$chain_dir/$chain_file" $bam "$(echo $(basename $bam) | cut -f 1 -d '.')_lifted" 
# done

## lift this bam locally only
CrossMap.py bam -a "$chain_dir/$chain_file" merged_blacklist_removed.bam  merged_blacklist_removed_lifted

## then run these commands on the hg38 dir
mv  merged_blacklist_removed_lifted.sorted.bam chimp_hg38_merged_blacklist_removed_lifted.sorted.bam
cp merged_blacklist_removed.bam human_hg38_merged_blacklist_removed.bam
samtools merge -c -p chimp_human_hg38_blacklist_removed.bam  human_hg38_merged_blacklist_removed.bam chimp_hg38_merged_blacklist_removed_lifted.sorted.bam



# for SRA in SRR8176431 SRR8176432 SRR8176433 SRR8176434 SRR8176435 SRR8176436 SRR8176437 SRR8176438 SRR8176439 SRR8176440 SRR8176441;
# do
#    samtools view -h merged_blacklist_removed_lifted.sorted.bam | awk -v S=$SRA '($0 ~ /^@/ || substr($1,1,length(v))==v)' | samtools view -S  -b  -o ${SRA}_blacklist_removed_lifted.sorted.bam - 
# done

# samtools merge -c -p C3647_blacklist_removed_lifted.sorted.bam SRR8176431_blacklist_removed_lifted.sorted.bam SRR8176432_blacklist_removed_lifted.sorted.bam SRR8176433_blacklist_removed_lifted.sorted.bam
# mv SRR8176434_blacklist_removed_lifted.sorted.bam C3649_blacklist_removed_lifted.sorted.bam
# samtools merge -c -p C3651_blacklist_removed_lifted.sorted.bam SRR8176435_blacklist_removed_lifted.sorted.bam SRR8176436_blacklist_removed_lifted.sorted.bam
# mv SRR8176437_blacklist_removed_lifted.sorted.bam C4955_blacklist_removed_lifted.sorted.bam
# mv SRR8176438_blacklist_removed_lifted.sorted.bam C8861_blacklist_removed_lifted.sorted.bam
# samtools merge -c -p C40280_blacklist_removed_lifted.sorted.bam SRR8176439_blacklist_removed_lifted.sorted.bam SRR8176440_blacklist_removed_lifted.sorted.bam SRR8176441_blacklist_removed_lifted.sorted.bam
# rm SRR*


## liftover bams after alignment
for SRA in  *.nochrM.encodefiltered.fixmate.rmorphanread.nodup.bam;do
    CrossMap.py bam -a ../../../../data/LiftOver_chains/panTro5ToHg38.over.chain $SRA  ${SRA}.lifted
done

