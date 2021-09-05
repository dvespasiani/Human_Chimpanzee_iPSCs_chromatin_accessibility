import os 
import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

stat_file = pd.read_csv(input_file,sep='\t',header=None)
stat_file.columns = ['File','Number']
stat_file['File'] = stat_file['File'].str.replace(r'.bam|.gz', '')
stat_file['File']= stat_file['File'].str.split('/').str[-1].str.strip()

from pandas import ExcelWriter
writer = ExcelWriter(output_file)
stat_file.to_excel(writer,'Sheet1',index=False)
writer.save()
