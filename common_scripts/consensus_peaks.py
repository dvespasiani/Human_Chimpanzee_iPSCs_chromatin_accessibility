## script used to find consensus peak using pybedtools
## see this thread https://www.biostars.org/p/13516/

import os
import pybedtools
import pandas as pd
import matplotlib.pyplot as plt

# input_dir = snakemake.input[0]
# output_dir = snakemake.output[0]
# fraction_overlap = snakemake.params[0]
# species = snakemake.param[1]

# peak_files = [x for x in os.listdir(input_dir) if x.endswith('filtered.narrowPeak.gz')]
# peak_files = [input_dir + x for x in peak_files]

fraction_overlap = 0.5

input_dir = os.getcwd() + '/output/PeakCalling/Files/'
output_dir = os.getcwd() + '/output/PeakCalling/'

peak_files = [x for x in os.listdir(input_dir) if x.endswith('filtered.narrowPeak.gz')]
peak_files = [input_dir + x for x in peak_files]


## verbosely specify bedtools objects....
sample1 = pybedtools.BedTool(peak_files[0])
sample2 = pybedtools.BedTool(peak_files[1])
sample3 = pybedtools.BedTool(peak_files[2])
sample4 = pybedtools.BedTool(peak_files[3])
sample5 = pybedtools.BedTool(peak_files[4])
sample6 = pybedtools.BedTool(peak_files[5])


## ...and make the combinations
colnames = ['chrom','start','end']

def peak_overlap (p1,p2):
    ## using reciprocal overlap here also reduces the number of combinations you'll need to specify downstream
    overlap = p1.intersect(p2, f=fraction_overlap, r=True)
    return(overlap)


##
s12 = peak_overlap(sample1,sample2)
s13 = peak_overlap(sample1,sample3)
s14 = peak_overlap(sample1,sample4)
s15 = peak_overlap(sample1,sample5)
s16 = peak_overlap(sample1,sample6)

##
s23 = peak_overlap(sample2,sample3)
s24 = peak_overlap(sample2,sample4)
s25 = peak_overlap(sample2,sample5)
s26 = peak_overlap(sample2,sample6)

##
s34 = peak_overlap(sample3,sample4)
s35 = peak_overlap(sample3,sample5)
s36 = peak_overlap(sample3,sample6)

##
s45 = peak_overlap(sample4,sample5)
s46 = peak_overlap(sample4,sample6)

##
s56 = peak_overlap(sample5,sample6)


combined_peaks = s12.cat(s13,s14,s15,s16,s23,s24,s25,s26,s34,s35,s36,s45,s46,s56)
combined_peaks = combined_peaks.sort()


## report table with number of consensus peaks per fraction of overlap
## and a table with number of consensus peaks after merging (see below)

species = 'human'

qc_number_peaks = pd.DataFrame({
    'species':[species],
    'fraction overlap':[fraction_overlap],
    'number consensus peaks':[len(combined_peaks.to_dataframe())]
    })

qc_number_peaks.columns = [''] * len(qc_number_peaks.columns)

qc_number_peaks.to_csv(output_dir + '/QC/numb_consensus_peaks.txt',sep='\t',header=True,index=None,mode='a')


dist_10bp = combined_peaks.merge(d=10)
dist_20bp = combined_peaks.merge(d=20)
dist_50bp = combined_peaks.merge(d=50)
dist_70bp = combined_peaks.merge(d=70)
dist_100bp = combined_peaks.merge(d=100)
dist_150bp = combined_peaks.merge(d=150)


qc_number_merged_peaks = pd.DataFrame({
    'species':[species]*6,
    'max merging distance btwn peaks':[10,20,50,70,100,150],
    'number consensus peaks':[
        len(dist_10bp.to_dataframe()),
        len(dist_20bp.to_dataframe()),
        len(dist_50bp.to_dataframe()),
        len(dist_70bp.to_dataframe()),
        len(dist_100bp.to_dataframe()),
        len(dist_150bp.to_dataframe())
        ]
    })

qc_number_merged_peaks.to_csv(output_dir + '/QC/numb_merged_consensus_peaks_'+str(fraction_overlap)+'.txt',sep='\t',header=True,index=None)

combined_peaks = combined_peaks.merge(d=50)

## NB distance between peaks can be calculated as:
## (distance between consecutive peak ends) - (peak width)
# test = combined_peaks.to_dataframe()
# test['diff_end'] = test.groupby('chrom').end.diff()
# test = test[test['diff_end'].notna()]
# test['diff_end'] = test['diff_end'].astype(int)
# test['distance_between_peaks'] = test['diff_end'] - (test['end'] -test['start'])
# # test = test[test['distance_between_peaks'] == min(test['distance_between_peaks'])]
# print(len(test.index))

##------------------------
## Write consensus peaks
##------------------------
final_peaks  = combined_peaks.to_dataframe()
final_peaks.to_csv(output_dir+'/ConsensusPeaks/Third_way_consensus_peak.bed', header=True, index=None, sep='\t')