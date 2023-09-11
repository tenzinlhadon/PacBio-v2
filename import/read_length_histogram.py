import os
import glob 
import pandas as pd
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/' 
import matplotlib.pyplot as plt
import sys
sample_arg = sys.argv[1]
sample_list = glob.glob("*_read_length.txt")
for sample in sample_list: 
    sample_name = sample.split("_read_length.txt")[0]
    ref = sample_name.replace(sample_arg,"")
    read_length_file = sample_name + "_read_length.txt"
    data = pd.read_csv(read_length_file , header = None, delim_whitespace=True)
    x = data[0]
    binwidth = 50
    plt.hist(x, bins=range(min(x),max(x) + binwidth, binwidth),edgecolor="blue", color="white")
    plt.title(ref)
    plt.xlabel('Read length')
    plt.ylabel('No of reads') 
    histogram_plot = sample_name + "_read_length_histogram.png"
    plt.savefig(histogram_plot)
    plt.show()
    df = pd.read_csv(read_length_file , header = None, delim_whitespace=True)
    bin_list = list(range(0,max(df[0]) + 100, 100))
    df_count = pd.DataFrame(df[0].value_counts(bins = bin_list))
    read_length_range_file = sample_name + "_read_length_range.csv"
    df_count['Read length'] = df_count.index
    df_count.columns = ['Read count', 'Read length']
    df_count = df_count[['Read length', 'Read count']]
    df_count.to_csv(read_length_range_file, index = False)
