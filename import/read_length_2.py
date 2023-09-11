import pandas as pd
import os
def check_and_write_row(csv_file, row_data):
    if not os.path.isfile(csv_file):
        # Create an empty DataFrame with the appropriate column names
        df = pd.DataFrame(columns=['Sample', f'Reads ({min_length}-{max_length})', 'Total Reads', 'Percentage'])
    else:
        df = pd.read_csv(csv_file)
    # Check if the row data already exists in the DataFrame
    row_exists = df.apply(lambda row: all(row == row_data), axis=1).any()
    if not row_exists:
        df = df.append(pd.Series(row_data, index=df.columns), ignore_index=True)
        df.to_csv(csv_file, index=False)
        print("Row added successfully!")
    else:
        print("Row already exists!")
# Usage example:
csv_file_path = 'data.csv'
row_data_to_add = ['Sample', f'Reads ({min_length}-{max_length})', 'Total Reads', 'Percentage']
check_and_write_row(csv_file_path, row_data_to_add)