## use this script to create a spreadsheet with file name and counts (either numb peaks or numb reads)
## run this script from QC directory
#!/bin/sh

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

## run command
echo "Getting filenames and count"
echo $input_dir
for file in "$wd/$input_dir/$@"; do
    if [[ "$file" = *.gz ]]; then
        printf '%s\n' "$file"  >> names.out
        gunzip -nc "$file" | wc -l >> counts.out
    elif [["$file" = *.bam]]; then
        printf '%s\n' "$file"  >> names.out
        samtools view "$file" | wc -l >> counts.out 
    fi
done

# for file in "$@"; do
#     if [[ "$file" = *.gz ]]; then
#         printf '%s\n' "$file"  >> names.out
#         gunzip -nc "$file" | wc -l >> counts.out
#     elif [["$file" = *.bam]]; then
#         printf '%s\n' "$file"  >> names.out
#         samtools view "$file" | wc -l >> counts.out 
#     fi
# done

paste -d "\t" names.out counts.out > stat_qc.out
rm names.out 
rm counts.out 

## Call custom python script to make final spreadsheet

echo "Making spreadsheet"
python "$wd"/common_scripts/make_spreadsheet.py stat_qc.out $chain_file

rm stat_qc.out 