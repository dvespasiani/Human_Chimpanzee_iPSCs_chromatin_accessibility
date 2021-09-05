module load bedtools/2.27.1 
module load ucsc/21072020


## split the file from eres et al into one for human and one for chimp 
## based on the species the tad was originally discovered
## and get the coordinates of the orthologous tad 
awk -F'\t' '$7 == "Human" ' human_chimp_orth_tads.txt | tail -n +1 > tmp_human.bed ## remove header for human
awk -v OFS='\t' '{ print $1, $2, $3, $7}' tmp_human.bed > human_hg38_tads.bed ## this will have hg38 coord

awk -F'\t' '$7 == "Chimp" ' human_chimp_orth_tads.txt > tmp_chimp.bed
awk -v OFS='\t' '{ print $1, $2, $3, $7}' tmp_chimp.bed > chimp_hg38_tads.bed ## this will have hg38 coord

rm tmp*

## intersection to find common tads in hg38coord
bedtools intersect -c -f 0.9 -r  -a human_hg38_tads.bed -b chimp_hg38_tads.bed > shared_hg38_tads.bed


## do the same for pantro5 
awk -F'\t' '$7 == "Human" ' human_chimp_orth_tads.txt | tail -n +1 > tmp_human.bed ## remove header for human
awk -v OFS='\t' '{ print $4, $5, $6, $7}' tmp_human.bed > human_pantro5_tads.bed ## this will have pantro5 coord

awk -F'\t' '$7 == "Chimp" ' human_chimp_orth_tads.txt > tmp_chimp.bed
awk -v OFS='\t' '{ print $4, $5, $6, $7}' tmp_chimp.bed > chimp_panto5_tads.bed ## this will have pantro5 coord

rm tmp*

bedtools intersect -c -f 0.9 -r  -a chimp_panto5_tads.bed -b human_pantro5_tads.bed > shared_pantro5_tads.bed

## shared files contain the list of tads identified in the relevant species in the species-specific coordinates
## so shared_pantro5_tads will have labelled the tads identified in chimp that are either species specific or in common
# head  shared_pantro5_tads.bed
# head  shared_hg38_tads.bed

## if you then count the words in each file and calculate common tads/(common tads + human specific + chimp specific) this matches ~ 43% claimed by Ittai Eres


## so now make multiples file:
## 1)common tads in hg38 coord
## 2) human specific tads in hg38 coord
## 3) chimp specific tads in pantro5 coords
field_names="seqnames\tstart\tend\tdisc_species\tlabel"
echo -e "$field_names" | cat - shared_hg38_tads.bed > hg38_tads.bed
echo -e "$field_names" | cat - shared_pantro5_tads.bed > pantro5_tads.bed

# awk -F'\t' '$5 != "0" ' shared_hg38_tads.bed > tmp1.bed | echo -e "$field_names" | cat - tmp1.bed > common_hg38_tads.bed
# awk -F'\t' '$5 == "0" ' shared_hg38_tads.bed > tmp2.bed | echo -e "$field_names" | cat - tmp2.bed > human_specific_tads.bed
# awk -F'\t' '$5 == "0" ' shared_pantro5_tads.bed > tmp3.bed | echo -e "$field_names" | cat - tmp3.bed > chimp_specific_tads.bed

rm chimp_panto5_tads.bed 
rm chimp_hg38_tads.bed
rm human_pantro5_tads.bed 
rm human_hg38_tads.bed
rm tmp*
rm shared*







