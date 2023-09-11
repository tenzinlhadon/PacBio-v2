import os
import pandas as pd
import sys 

sample = sys.argv[1]

# Output file to store merged values
output_file = sample + "_alignment_per_ref_metrics.csv"

# Get a list of input files
input_files = [f for f in os.listdir() if f.endswith(".stats")]

# Create an empty DataFrame to store the merged values
merged_df = pd.DataFrame()

# Loop through input files and merge values
for input_file in input_files:
    df = pd.read_csv(input_file)    
    if merged_df.empty:
        merged_df = df
    else:
        merged_df = merged_df.merge(df, on='Metrics', how='left')

# shift combined_ref to the third column 
columns = merged_df.columns.tolist()
columns.remove('combined_ref')
columns.insert(1,'combined_ref')
merged_df = merged_df.reindex(columns=columns)


# Write merged DataFrame to the output file
merged_df.to_csv(output_file, index=False)

print(f"Merged values written to {output_file}.")