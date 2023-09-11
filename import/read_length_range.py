import os
import glob 
import pandas as pd 
import matplotlib.pyplot as plt
sample_list = glob.glob("*_read_length.txt")
#samtools view test.bam | awk '{print length($10)}'
for sample in sample_list: 
    sample_name = sample.split("_read_length.txt")[0]
    read_length_file = sample_name + "_read_length.txt"
    df = pd.read_csv(read_length_file , header = None, delim_whitespace=True)
    bin_list = list(range(0,max(df[0]) + 100, 100))
    df_count = pd.DataFrame(df[0].value_counts(bins = bin_list))
    read_length_range_file = sample_name + "_read_length_range.txt"
    df_count.to_csv(read_length_range_file)

