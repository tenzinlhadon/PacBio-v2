import sys
import shutil 
import pandas as pd 
import numpy as np
import os 

file_arg = sys.argv[1] + sys.argv[2] + "_counts.temp.csv"
read_count = int(sys.argv[3])
print(file_arg)
sample = sys.argv[1]
print("sample")
print(sample)
output_file = file_arg.replace(".temp","")
output_file_percentage = sample + sys.argv[2] + "_reference_percentage.csv"
with open(
    file_arg, 'r') as r, open(
        output_file, 'w') as o:
      
    for line in r:
        #strip() function
        if line.strip():
            o.write(line)
def sort_percentage(file_name):
    df = pd.read_csv(file_name, sep = ",")
    df['percent'] = (df['aligned_reads'] / read_count) * 100
    df = df.sort_values(['aligned_reads'], ascending = [False])
    df.to_csv(file_name, index=False)
    df.columns = ['reference','aligned_reads',sample]
    df_new = df[['reference',sample]]
    df_new = df_new[df_new['reference'] != 'combined_ref']
    df_new.to_csv(output_file_percentage, index = False)
    print(df)

sort_percentage(output_file)

multisample_contamination_file = "/analysis/sample_files/Contamination_percentage" + sys.argv[2] + ".csv"

if os.path.exists(multisample_contamination_file ):
    df_sample_wise = pd.read_csv(output_file_percentage )
    df_multisample = pd.read_csv(multisample_contamination_file  )
    print(df_sample_wise)
    print(df_multisample)
    merged_df = pd.merge(df_multisample,df_sample_wise, on='reference', how='inner')
    merged_df = merged_df[merged_df['reference'] != 'combined_ref']
    print(merged_df)
    temp_file = "/analysis/sample_files/Contamination_percentage" + sys.argv[2] + "_temp.csv"
    merged_df.to_csv(temp_file, index = False)
    shutil.move(temp_file,multisample_contamination_file )
else:
    df_sample_wise = pd.read_csv(output_file_percentage)
    df_sample_wise.to_csv(multisample_contamination_file  , index = False)





