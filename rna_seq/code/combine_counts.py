## script used to combine counts across samples to generate a single file per species
import os 
import pandas as pd
import sys
from functools import reduce

input_dir = snakemake.input[0]
output_file = snakemake.output[0]

# input_dir = os.getcwd() + '/rna_seq/pantro5/output/PostAlignment/'

count_files = [x for x in os.listdir(input_dir) if x.endswith('_final.txt')]
count_files = [input_dir + x for x in count_files]

counts = []
for filename in count_files:
    df = pd.read_csv(filename, index_col=None,sep='\t')
    counts.append(df)

combined_counts = reduce(lambda x, y: pd.merge(x, y, left_index=True, on=['GeneID'], right_index=True), counts)

sample_names = [os.path.basename(n).rsplit('_')[0] for n in os.listdir(input_dir) if n.endswith('_final.txt')]
sample_names.insert(0, "GeneID")

combined_counts.columns =  sample_names

combined_counts.to_csv(output_file,sep ='\t',index=False,header=True)

