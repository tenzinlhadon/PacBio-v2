import subprocess
import pandas  as pd
import numpy as np


Read the bed file into a pandas DataFrame

sample = "m64409e_230203_155352"
df = pd.read_csv(bed_file, sep='\t', header=None)
df[1] = df[1] + 1
# Define a function to calculate the clipped bases
def get_clipped_bases(cigar):
    bases = 0
    for i in range(len(cigar)):
        if cigar[i].isdigit():
            bases = bases * 10 + int(cigar[i])
        else:
            if i < len(cigar) - 1 and cigar[i] == 'M':
                return 1
            elif cigar[i] == 'H' or cigar[i] == 'S':
                return bases + 1
    return bases if bases > 0 else 1
# Apply the function to calculate clipped bases and create a new column
df['Clipped Bases'] = df[6].apply(get_clipped_bases)
df['position'] = df[0] + ":" + df[1].astype(str) + "-" + df[2].astype(str) 
df.columns = ['ref', 'start', 'end', 'read','score', 'strand', 'cigar', 'read_start_pos', 'ref_coordinate']
df = df.sort_values(by=['read', 'read_start_pos'], ascending = True)

# for annotation 
df_annotation = df[['ref','start', 'end','read', 'read_start_pos']].copy()
df_annotation['start'] = df_annotation['start'] - 1
df_annotation.to_csv("chimera_annotation_input.bed", sep ="\t", header = None, index = False)

non_overlapping_regions = "bedtools subtract -a chimera_annotation_input.bed -b /analysis/features_BED_files/combined_annotation.temp.bed  > nonoverlappingregions.bed"
result = subprocess.run(non_overlapping_regions, shell = True)

overlapping_regions = "bedtools intersect -a chimera_annotation_input.bed -b /analysis/features_BED_files/combined_annotation.temp.bed  -wa -wb > annotated_unknown.bed"
result = subprocess.run(overlapping_regions, shell = True)

df_annotated = pd.read_csv("annotated_unknown.bed", sep = "\t", header = None)

# non overlapping region, fill coordinates 
mask = (df_annotated[5] == '.') & (df_annotated[6] == -1) & (df_annotated[7] == -1)
df_annotated.loc[mask, [6, 7, 8]] = df_annotated.loc[mask, [1, 2, 0]].values
# convert to 1 coordinate
df_annotated[1] = df_annotated[1] + 1
df_annotated[6] = df_annotated[6] + 1

# Iterate over each row of the DataFrame
interval1_start = df_annotated[1] 
interval1_end = df_annotated[2]
interval2_start = df_annotated[6] 
interval2_end = df_annotated[7]
overlap_start = np.maximum(interval1_start, interval2_start)
overlap_end = np.minimum(interval1_end, interval2_end)
df_annotated[9] = overlap_start
df_annotated[10] = overlap_end

# Print the updated DataFrame
df_annotated[11] = df_annotated.apply(lambda row: f"{row[8]}:{row[9]}-{row[10]}", axis=1)
df_annotated[11] = df_annotated[11].str.strip()
df_annotated = df_annotated.iloc[:,[3,4,9,11]]
df_annotated.columns = ['read','read_start_pos','sub_region_start','annotation']
# df_annotated.to_csv("annotated_overlapped.csv", header = None, index = False)

df_non_overlapping = pd.read_csv("nonoverlappingregions.bed", sep = "\t", header = None)
# m64409e_230203_155352/1/ccs     1       4851    cassette-3'ITR                         :4851-4981
df_non_overlapping[1] = df_non_overlapping[1] + 1
df_non_overlapping[5] = df_non_overlapping.apply(lambda row: f"{row[0]}:{row[1]}-{row[2]}", axis=1)
df_non_overlapping[5] = df_non_overlapping[5].str.strip()
df_non_overlapping = df_non_overlapping.iloc[:,[3,4,1,5]]
df_non_overlapping.columns = ['read','read_start_pos','sub_region_start','annotation']
df_non_overlapping.to_csv("annotated_non_overlapping.bed", sep = "\t", header = None, index = False)

df_annotated_all = pd.concat([df_annotated,df_non_overlapping])

print(df_annotated_all.head())
df_annotated_all.columns = ['read','read_start_pos','sub_region_start','annotation']
df_annotated_all = df_annotated_all.sort_values(by=['read','read_start_pos','sub_region_start'], ascending = True)
# df_annotated_all.to_csv("df_known_unknown.csv")
df_annotated_grouped = df_annotated_all.groupby('read')['annotation'].agg('__'.join).reset_index()
df_annotated_grouped.to_csv("annotated_grouped_1.csv", index = False)
df_annotated_chimeric_species_count = df_annotated_grouped['annotation'].value_counts().reset_index()
df_annotated_chimeric_species_count.columns = ['Annotated_chimeric_species','Read_count']
df_annotated_chimeric_species_count = df_annotated_chimeric_species_count.sort_values(by = "Read_count", ascending = False)
df_annotated_chimeric_species_count = df_annotated_chimeric_species_count.applymap(lambda x: x.strip() if isinstance(x, str) else x)
total_read_count =  3190173
df_annotated_chimeric_species_count['Read_count_%'] = ((df_annotated_chimeric_species_count['Read_count']/total_read_count) * 100).round(2)
# df_annotated_chimeric_species_count.to_csv("annotated_chimeric_species_count.csv", index = None)
df_annotated_chimeric_species_count = df_annotated_chimeric_species_count[df_annotated_chimeric_species_count['Read_count'] >50]
Total_chimeric_reads = df_annotated_chimeric_species_count['Read_count'].sum()
Total_chimeric_percentage = df_annotated_chimeric_species_count['Read_count_%'].sum()
sum_row =pd.DataFrame({'Annotated_chimeric_species' : ['Total'], 'Read_count': [Total_chimeric_reads], 'Read_count_%': [Total_chimeric_percentage]})
df_annotated_chimeric_species_count = pd.concat([sum_row,df_annotated_chimeric_species_count], ignore_index = True)
df_annotated_chimeric_species_count.to_csv("annotated_chimeric_species_count_1.csv", index = None)

with open("annotated_chimeric_species_count_1.csv", 'r') as file:
    content = file.read()

output_csv = sample + "_annotated_chimeric_species_count_cassette_whole_2.csv"

content = content.replace(output_csv, output_csv)
with open("", 'w') as file:
    file.write(content)




