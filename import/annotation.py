import pandas as pd 
annotation_df= pd.read_csv("/analysis/features_BED_files/combined_annotation.temp.bed", sep = "\t", header = None)
annotation_df[3] = annotation_df[0] + "-" + annotation_df[3]
annotation_df.to_csv("/analysis/features_BED_files/combined_annotation.bed", header = False, index = False, sep = "\t")