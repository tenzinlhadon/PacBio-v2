import csv
import sys 
import pandas as pd 
import matplotlib.pyplot as plt


bed  = sys.argv[1]
print(bed)
read_count = int(sys.argv[2])

sample = bed.split(".bed")[0]

def convert_zero_one_bed_coordinate(input_file):
    output_file = sample + ".csv"
    with open(input_file, 'r') as bed_file, open(output_file, 'w', newline='') as csv_file:
        bed_reader = csv.reader(bed_file, delimiter='\t')
        csv_writer = csv.writer(csv_file)
        for row in bed_reader:
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            rest = row[3:]
            start_1based = start + 1
            modified_row = [chrom, str(start_1based), str(end)] + rest
            csv_writer.writerow(modified_row)

convert_zero_one_bed_coordinate(bed)



def end_points_calculation(input_file,end):
    df = pd.read_csv(input_file)
    if end == "Start": 
        unique_values = df.iloc[:,1].value_counts().reset_index()
    else:
        unique_values = df.iloc[:,2].value_counts().reset_index()
    unique_values.columns = ['Position','Frequency']
    total_count = unique_values['Frequency'].sum()
    unique_values['Percentage'] = (unique_values['Frequency']/read_count) * 100
    unique_values.to_csv(sample + "_" +  end + "_count.csv")


one_based_coordinate = sample + ".csv"
end_points_calculation(one_based_coordinate, "Start")
end_points_calculation(one_based_coordinate , "End")

start_frequencey = sample + "_" +  "Start" + "_count.csv"
end_frequency = sample + "_" +  "End" + "_count.csv"
df_start = pd.read_csv(start_frequencey)
df_end = pd.read_csv(end_frequency)
df_combined = pd.concat([df_start,df_end])
df_combined = df_combined.groupby('Position').agg({'Frequency':'sum', 'Percentage':'sum'}).reset_index()
df_combined = df_combined.sort_values(by = "Frequency", ascending = False)
df_combined.to_csv(sample + "_combined_read_counts.csv", index = False)
plt.scatter(df_combined['Position'], df_combined['Frequency'], s=5)
plt.xlabel('Position')
plt.ylabel('Read count')
plt.title(sample + "_" +  bed )
plt.grid(True)
top_values = df_combined.head(10)
for i, row in top_values.iterrows():
    plt.text(row['Position'], row['Frequency'], f"({row['Position']})", fontsize=10)

# Create a list to hold scatter plot artists
scatter_artists = []


for i, row in top_values.iterrows():
    label = f"{row['Position']}"
    scatter_artist = plt.scatter(row['Position'], row['Frequency'], label=label)
    scatter_artists.append(scatter_artist)
# Add legend with scatter plot artists
legend_labels = [f"({row['Percentage']:.2f}%) - {row['Position']}" for i, row in top_values.iterrows()]
plt.legend(handles=scatter_artists, labels=legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Position')
plt.ylabel('Read count')
plt.title(sample)
# Add grid lines
plt.grid(True)
# Save the plot
plt.savefig(sample + "breakpoint_read_count.png", bbox_inches = "tight")
# Display the plot

