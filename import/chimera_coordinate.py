import pandas as pd
import glob
import sys 
import os 
total_reads = int(sys.argv[1])
sample = sys.argv[2]
annotation_file = sample + "_annotation.txt"
if (os.path.isfile(annotation_file)) == True: 
    sample = (glob.glob("*annotation.txt")[0]).split("_annotation.txt")[0]
    annotation_file = sample + "_annotation.txt"
    df = pd.read_csv(annotation_file, sep = "\t", header = None)
    df.columns = ['query_chr','query_start', 'query_end','query', 'subject_chr','subject_start', 'subject_end', 'subject']
    df = df[['query', 'subject']]
    list_1 = df['query']
    df2 = df[['query','subject']].groupby('query').apply(lambda x: x['subject'].unique())
    df2 = pd.DataFrame({'query':df2.index, 'subject':df2.values})
    df2['subject'] = df2['subject'].apply(sorted)
    df2['subject'] = ['|'.join(map(str, l)) for l in df2['subject']]
    df2 = df2.sort_values(by=['subject'], ascending=True)
    print(df2.head())
    df_counts = df2.subject.value_counts().to_frame()
    df_counts.columns = ['No_of_reads']
    df_counts['Chimeric_species'] = df_counts.index
    sum = df_counts['No_of_reads'].sum()
    df_counts["chimeric_reads_%"] = (df_counts['No_of_reads']/sum)*100
    df_counts['total_reads_%'] = (df_counts['No_of_reads']/total_reads)*100
    df_counts = df_counts[['Chimeric_species', 'No_of_reads', 'chimeric_reads_%', 'total_reads_%']]
    chimeric_species_file = sample + "_chimeric_species.csv"
    df_counts.to_csv(chimeric_species_file, index = False)