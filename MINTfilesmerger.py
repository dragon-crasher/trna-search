import pandas as pd
from functools import reduce
import os
import sys
def merge_mint_files(folder_path, set_name):
    """
    Merges multiple MINT files into a single DataFrame and saves it to a CSV file.

    Parameters:
    - folder_path: Path of folder containting the MINT files.
    - output_file: Path to the output CSV file.
    """
    # Initialize an empty list to store DataFrames
    dataframes = []

    # Loop through each file in the list
    for file in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path, sep='\t')
        # Select relevant columns
        filtered_df = df[['MINTbase Unique ID', 'tRF sequence', 'tRF type(s)', 'Unnormalized read counts']]
        # Rename the last column to the sample name
        sample_name = file.split('/')[-1].split('.')[0]
        filtered_df = filtered_df.rename(columns={'Unnormalized read counts': sample_name})
        dataframes.append(filtered_df)

    # Merge all DataFrames on the first three columns
    
    merged_df = reduce(lambda left, right: pd.merge(left, right, on=['MINTbase Unique ID', 'tRF sequence', 'tRF type(s)'], how='outer'), dataframes)

    output_file = os.path.join('/mnt/d/bioinformatics/RNAseq_pipeline/data/', f'{set_name}merged_mint_files.csv')

    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved to: {output_file}")


if len(sys.argv) != 3:
    print("Usage: python MINTfilesmerger.py <folder_path> <set_name>")
    sys.exit(1)
folder_name = sys.argv[1]
set_name = sys.argv[2]

folder_path = os.path.join('/mnt/d/bioinformatics/MINT/outputs/', folder_name)
merge_mint_files(folder_path, set_name)