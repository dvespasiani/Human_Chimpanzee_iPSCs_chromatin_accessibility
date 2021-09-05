#!/bin/bash

#SBATCH --partition=mig
#SBATCH --nodes=1
#SBATCH --job-name="name"
#SBATCH --account="punim0586"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50000
#SBATCH --time=1-23:00:00
#SBATCH --mail-user=dvespasiani@student.unimelb.edu.au
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# The modules to load:
source /usr/local/module/spartan_new.sh
module load web_proxy
module load gcc/8.3.0 openmpi/3.1.4
module load python/3.7.4 
module load r/4.0.0  
module load deeptools/3.3.1-python-3.7.4
module load samtools/1.9

## script (needs to be re-written)
cd /data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/pantro5/output/Post_alignment/Files/combined/Lifted/

for SRA in  *.lifted.sorted.bam;do
    output_file=$(echo $(basename $SRA)| cut -f 1 -d '.')_lifted_tn5_shifted.bam
    alignmentSieve -b $SRA -o $output_file --ATACshift
done

for SRA in  *_lifted_tn5_shifted.bam;do
    output_file=$(echo $(basename $SRA)| cut -f 1 -d '_')
    samtools sort $SRA -O BAM -o ${output_file}_lifted_tn5_shifted_sorted.bam
done

# for SRA in *_lifted_tn5_shifted_sorted.bam;do
#     samtools index  $SRA ${SRA}.bai
# done

samtools merge -c -p C3647_lifted_tn5_shifted_sorted.bam SRR8176431_lifted_tn5_shifted_sorted.bam SRR8176432_lifted_tn5_shifted_sorted.bam SRR8176433_lifted_tn5_shifted_sorted.bam
mv SRR8176434_lifted_tn5_shifted_sorted.bam C3649_lifted_tn5_shifted_sorted.bam
samtools merge -c -p C3651_lifted_tn5_shifted_sorted.bam SRR8176435_lifted_tn5_shifted_sorted.bam SRR8176436_lifted_tn5_shifted_sorted.bam
mv SRR8176437_lifted_tn5_shifted_sorted.bam C4955_lifted_tn5_shifted_sorted.bam
mv SRR8176438_lifted_tn5_shifted_sorted.bam C8861_lifted_tn5_shifted_sorted.bam
samtools merge -c -p C40280_lifted_tn5_shifted_sorted.bam SRR8176439_lifted_tn5_shifted_sorted.bam SRR8176440_lifted_tn5_shifted_sorted.bam SRR8176441_lifted_tn5_shifted_sorted.bam

for file in  C*;do
    samtools index -c $file ${file}.bai
done

# samtools merge -c -p H19101_tn5_shifted_sorted.bam SRR8176442_tn5_shifted_sorted.bam  SRR8176443_tn5_shifted_sorted.bam
# samtools merge -c -p H20961_tn5_shifted_sorted.bam SRR8176444_tn5_shifted_sorted.bam  SRR8176445_tn5_shifted_sorted.bam
# samtools merge -c -p H28834_tn5_shifted_sorted.bam SRR8176446_tn5_shifted_sorted.bam  SRR8176447_tn5_shifted_sorted.bam
# samtools merge -c -p HUtt45_tn5_shifted_sorted.bam SRR8176448_tn5_shifted_sorted.bam  SRR8176449_tn5_shifted_sorted.bam
# samtools merge -c -p HUtt60_tn5_shifted_sorted.bam SRR8176450_tn5_shifted_sorted.bam  SRR8176451_tn5_shifted_sorted.bam
# samtools merge -c -p H19114_tn5_shifted_sorted.bam SRR8176452_tn5_shifted_sorted.bam  SRR8176453_tn5_shifted_sorted.bam


chimp_dir=/data/projects/punim0595/dvespasiani/Human_Chimpanzee_iPSCs_chromatin_accessibility/pantro5/output/Post_alignment/Files/combined/Lifted

samtools merge -c -p chimp_human_hg38_tn5_shifted_merged.bam  \
H19101_tn5_shifted_sorted.bam  H20961_tn5_shifted_sorted.bam H28834_tn5_shifted_sorted.bam \
HUtt45_tn5_shifted_sorted.bam HUtt60_tn5_shifted_sorted.bam H19114_tn5_shifted_sorted.bam \
${chimp_dir}/C3647_lifted_tn5_shifted_sorted.bam ${chimp_dir}/C3649_lifted_tn5_shifted_sorted.bam \
${chimp_dir}/C3651_lifted_tn5_shifted_sorted.bam ${chimp_dir}/C4955_lifted_tn5_shifted_sorted.bam \
${chimp_dir}/C40280_lifted_tn5_shifted_sorted.bam ${chimp_dir}/C8861_lifted_tn5_shifted_sorted.bam



