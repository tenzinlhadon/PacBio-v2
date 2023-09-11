import os
import glob 
import pandas as pd 
import numpy as np
import sys 
sample_list = glob.glob("*aln.temp.txt")
read_count = int(sys.argv[1])


for sample in sample_list: 
    sample_name = sample.split("_consensus_aln.temp.txt")[0]
    df = pd.read_csv(sample, sep = "\t", header = None)
    df.columns = [ "Cluster","size","configuration","%_identity","alignment_length","mismatches","gap_openings","query_start","query_end","subject_start","subject_end","E_value","bit_score"]  
    df = df.sort_values(by=['Cluster','alignment_length', '%_identity'], ascending = [True,False, False])
    # Keep only the first occurrence of duplicates in Cluster
    df = df.drop_duplicates(subset='Cluster', keep='first')
    df['size%'] = (df['size']/read_count) * 100
    df_clus = df[df['size'] > 10] #take cluster with read support greater than 10
    df_clus = df_clus.sort_values(['Cluster','alignment_length','%_identity','size'], ascending = [True,False, False, False])
    df_clus = df_clus.sort_values(['size','alignment_length', '%_identity'], ascending = [False,False, False])
    # df_clus = df_clus.drop(['size%'], axis=1, errors='ignore')
    cluster_file= sample_name + "_consensus_aln.csv"
    df_clus.to_csv(cluster_file, index=False)


