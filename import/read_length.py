from Bio import SeqIO
import csv
import os
import sys

# Define the read length range
min_length = int(sys.argv[1])
max_length = int(sys.argv[2])
fastq_file = sys.argv[3]

# Define the current directory path
current_directory = os.getcwd()

# Find the FASTQ file ending with '.fastq'
# fastq_file = next(file for file in os.listdir(current_directory) if file.endswith('.fastq'))
sample = fastq_file.split('.fastq')[0]



# Read the FASTQ file and calculate the read length statistics
total_reads = 0
target_reads = 0

for record in SeqIO.parse(fastq_file, 'fastq'):
    total_reads += 1
    read_length = len(record.seq)

    if min_length <= read_length <= max_length:
        target_reads += 1

# Calculate the percentage of reads in the target range
if target_reads == 0:
    percentage = 0
else:
    percentage = (target_reads / total_reads) * 100

def create_csv_with_header(csv_file, header_row):
    if not os.path.exists(csv_file):
        with open(csv_file, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(header_row)
            writer.writerow([sample, target_reads, f'{percentage:.2f}%'])
            print("CSV file created with header.")
    else:
        with open(csv_file, 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([sample, target_reads, f'{percentage:.2f}%'])

csv_file_path = '/analysis/sample_files/Intact_cassette_percentage.csv'
header_row_to_add = ['Sample', f'Reads with read lengths: ({min_length}-{max_length})', 'Percentage']
create_csv_with_header(csv_file_path, header_row_to_add)


# Close the CSV file


