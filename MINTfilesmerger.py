import pandas as pd
from functools import reduce
import os
import sys

def merge_mint_files(folder_path, set_name):
    """
    Merges multiple MINT files into a single DataFrame and saves it to a CSV file.

    Parameters:
    - folder_path: Path of folder containing the MINT files.
    - set_name: Prefix for the output merged CSV file.
    """
    dataframes = []

    for file in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file)

        # Check number of lines in the file
        with open(file_path, 'r') as f:
            line_count = sum(1 for _ in f)

        if line_count < 7:
            print(f"Skipping file '{file}' because it has less than 7 lines ({line_count} lines).")
            continue

        # Read file skipping first 5 lines
        df = pd.read_csv(file_path, sep='\t', skiprows=6)
        print(f"Columns in {file}: {df.columns.tolist()}")

        # Select relevant columns
        filtered_df = df[['License Plate', 'tRF sequence', 'tRF type(s)', 'Unnormalized read counts']]

        # Use filename without extension as sample name
        sample_name = os.path.splitext(file)[0]
        filtered_df = filtered_df.rename(columns={'Unnormalized read counts': sample_name})
        dataframes.append(filtered_df)

    if not dataframes:
        print("No valid files to merge.")
        return

    # Merge all DataFrames on the first three columns
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on=['License Plate', 'tRF sequence', 'tRF type(s)'], how='outer'),
        dataframes
    )

    output_file = os.path.join('/raid/anirudh/bioinformatics/RNAseq_pipeline/data/', f'{set_name}merged_mint_files.csv')
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python MINTfilesmerger.py <folder_path> <set_name>")
        sys.exit(1)

    folder_name = sys.argv[1]
    set_name = sys.argv[2]

    folder_path = os.path.join('/raid/anirudh/bioinformatics/MINT/outputs/', folder_name)
    merge_mint_files(folder_path, set_name)
