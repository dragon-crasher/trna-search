import pandas as pd
from functools import reduce
import os
import sys
import re

def merge_mint_files(folder_path, folder_name):
    dataframes = []

    for file in os.listdir(folder_path):
        if file.endswith('.txt') and 'exclusive' in file:
            file_path = os.path.join(folder_path, file)

            with open(file_path, 'r') as f:
                line_count = sum(1 for _ in f)

            if line_count < 8:
                print(f"Skipping file '{file}' because it has less than 7 lines ({line_count} lines).")
                continue

            df = pd.read_csv(file_path, sep='\t', skiprows=6)
            print(f"Columns in {file}: {df.columns.tolist()}")
            
            #filtered_df = df[['MINTbase Unique ID', 'tRF sequence', 'tRF type(s)', 'Unnormalized read counts']]
            filtered_df = df[['License Plate', 'tRF sequence', 'tRF type(s)', 'Unnormalized read counts']]

            sample_ne = os.path.splitext(file)[0]
            sample_name = re.split(r'[._]',sample_ne)
            filtered_df = filtered_df.rename(columns={'Unnormalized read counts': sample_name})
            dataframes.append(filtered_df)

    if not dataframes:
        print("No valid files to merge.")
        sys.exit(1)

    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on=['License Plate', 'tRF sequence', 'tRF type(s)'], how='outer'),
        dataframes
    )

    output_dir = '/raid/anirudh/bioinformatics/RNAseq_pipeline/data/'
    output_file = os.path.join(output_dir, f'{folder_name}merged_mint_files.csv')
    merged_df.to_csv(output_file, index=False)

    print(f"Merged file saved to: {output_file}")  # Important: this line is parsed by bash

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MINTfilesmerger.py <folder_name>")
        sys.exit(1)

    folder_name = sys.argv[1]
    folder_path = os.path.join('/raid/anirudh/bioinformatics/MINT/outputs/', folder_name)

    if not os.path.isdir(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        sys.exit(1)

    merge_mint_files(folder_path, folder_name)
