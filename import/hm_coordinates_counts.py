import pandas as pd
import subprocess
import os
from os import listdir
import shutil
import numpy as np
import json

import glob
alignment_csv = '/analysis/reference_masking_output/Homologous_sequences/'
# Get a list of all CSV files in the directory
csv_files = glob.glob(alignment_csv + '*_homologous_regions.csv')
# Combine all CSV files into a single DataFrame
combined_data = pd.concat([pd.read_csv(file) for file in csv_files])
# Write the combined data to a new CSV file
combined_data.to_csv('/analysis/reference_masking_output/Homologous_sequences/combined.csv', index=False)


###########################################################################################
#Find homolgous sequences present between all the reference genomes
##########################################################################################

file_name = "/analysis/reference_masking_output/Homologous_sequences/combined.csv"
colnames = ['Query','Subject','Query start','Query end']
df = pd.read_csv(file_name, usecols = colnames)
df['Query'] = df['Query'].astype('string')
df['Subject'] = df['Subject'].astype('string')
df['Subject'] = df['Subject'].str.split("'").str[1] + df['Subject'].str.split("'").str[2]
df['Query'] = df['Query'].str.split("'").str[1] + df['Query'].str.split("'").str[2]
df = pd.DataFrame(df, columns = ['Query','Subject'])
df2 = df[['Query','Subject']].groupby('Query').apply(lambda x: x['Subject'].unique())
df2 = pd.DataFrame({'Query':df2.index, 'Subject':df2.values})
df2['Homologous_sequences'] =  df2.apply(lambda x:  np.append(x['Subject'], x['Query']), axis=1).to_numpy()
df2 = df2[['Query','Homologous_sequences']]

list_1 = df2['Homologous_sequences']
#print(list_1[:5])
list_dicts = []
for item in list_1:
    d = {}
    for b in item:
        i = b.split(':')
        if i[0] in d:
            d[i[0]].append(i[1])
        else:
    
            d[i[0]] = [i[1]]
      
       
        
    list_dicts.append(d)
df2['dict'] = list_dicts

df2['Homologous_sequences'] = [','.join(map(str, l)) for l in df2['Homologous_sequences']]

df3= df2[~df2.Homologous_sequences.str.split(",").apply(frozenset).duplicated(keep='last')].copy()
df3['Homologous_sequences_labels'] = ["Homologous sequence-"+str(i+1) for i in range(0,len(df3))]
df3.head()

df4 = df3[["Homologous_sequences_labels","dict"]]
df5 = df4["dict"].apply(pd.Series )
df5.index = df4['Homologous_sequences_labels']
ref = (df5.stack()
   .reset_index(level=1)
   .groupby(level=0, sort=False)
   ['level_1'].apply(list)
)
df5.insert(0, "Reference genomes carrying homologous sequences",ref )
df5.fillna("n/a",inplace = True)


##############################################################################
# write the data frame with homologous sequence coordinates into a .csv file 
# One of the two final result files 
############################################################################

path = "/analysis/count_reads"
path_hm_counts = "/analysis/reference_masking_output/Homologous_sequences/"
if os.path.exists(path_hm_counts) == False:
    os.mkdir(path_hm_counts)
homolgous_regions_coordinates = path_hm_counts + "Homologous_regions_coordinates.csv"
df5.to_csv(homolgous_regions_coordinates)



