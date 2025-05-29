import pandas as pd
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np

def find_header_line(file_path, header_keyword='License plate', max_lines=20):
    """
    Optimized header line finder using line cache and early termination.
    """
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            if header_keyword in line.lower():
                return i
    return None

def process_file(file_path):
    """
    Process a single file and return filtered DataFrame or None if invalid.
    """
    try:
        header_line = find_header_line(file_path)
        if header_line is None:
            return None

        # Read only required columns if possible
        df = pd.read_csv(
            file_path,
            sep='\t',
            skiprows=header_line,
            usecols=lambda col: col.strip() in {
                'License Plate', 
                'tRF sequence', 
                'tRF type(s)', 
                'Unnormalized read counts'
            },
            dtype={
                'License Plate': 'category',
                'tRF sequence': 'category',
                'tRF type(s)': 'category',
                'Unnormalized read counts': np.uint32
            },
            engine='c'
        )

        # Strip whitespace from column names
        df.columns = df.columns.str.strip()

        # Check required columns
        required_cols = ['License Plate', 'tRF sequence', 'tRF type(s)', 'Unnormalized read counts']
        if not all(col in df.columns for col in required_cols):
            return None

        # Get sample name from filename
        sample_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Rename and filter
        df = df[required_cols].rename(columns={'Unnormalized read counts': sample_name})
        df = df.dropna(subset=['License Plate', 'tRF sequence', 'tRF type(s)'])
        
        return df if not df.empty else None

    except Exception as e:
        print(f"Error processing {file_path}: {str(e)}")
        return None

def merge_mint_files(folder_path, folder_name):
    # Get all relevant files
    files = [
        os.path.join(folder_path, f) 
        for f in os.listdir(folder_path) 
        if f.endswith('.txt') and 'exclusive' in f
    ]
    
    if not files:
        print("No valid files found to process.")
        sys.exit(1)

    # Process files in parallel
    with ThreadPoolExecutor() as executor:
        results = executor.map(process_file, files)
        dataframes = [df for df in results if df is not None and not df.empty]

    if not dataframes:
        print("No valid dataframes to merge.")
        sys.exit(1)

    # Efficient merge using pd.concat after setting index
    indexed_dfs = [
        df.set_index(['License Plate', 'tRF sequence', 'tRF type(s)']) 
        for df in dataframes
    ]
    
    merged_df = pd.concat(indexed_dfs, axis=1).reset_index()

    print(f"Merged DataFrame shape: {merged_df.shape}")

    # Save output
    output_dir = '/raid/anirudh/bioinformatics/RNAseq_pipeline/data/'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{folder_name}_merged_mint_files.csv')
    
    # Use more efficient file format if possible
    merged_df.to_csv(output_file, index=False)
    print(f"Merged file saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python MINTfilesmerger.py <folder_name>")
        sys.exit(1)

    folder_name = sys.argv[1].rstrip('/\\')
    folder_path = os.path.join('/raid/anirudh/bioinformatics/MINT/outputs/', folder_name)

    if not os.path.isdir(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        sys.exit(1)

    merge_mint_files(folder_path, folder_name)
