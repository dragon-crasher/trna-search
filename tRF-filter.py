import pandas as pd
import numpy as np
import os
import shutil

def file_chooser():
    source_dir = '/mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0/outputs/CLL'  
    dest_dir = '/mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0/outputs/CLL/chosen'  

    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    for filename in os.listdir(source_dir):
        if 'exclusive' in filename and filename.endswith('.txt') and 'expression' in filename:
            source_file = os.path.join(source_dir, filename)
            new_filename = filename.split('-')[0] + '.txt'
            dest_file = os.path.join(dest_dir, new_filename)
            shutil.copy2(source_file, dest_file)  

    print("Files copied successfully!")

def read_and_merge_dataframes(folder_path, count_threshold=10, output_csv="merged_data.csv"):
    merged_df = None
    sample_names = []
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            df = pd.read_csv(file_path, sep='\t')
            filename_without_extension = os.path.splitext(filename)[0]

            required_columns = ['tRF sequence', 'Unnormalized read counts']
            if all(col in df.columns for col in required_columns):
                df = df[required_columns].copy()
            else:
                print(f"Warning: Required columns not found in {filename}. Skipping.")
                continue

            df = df[df['Unnormalized read counts'] >= count_threshold]
            df.rename(columns={'Unnormalized read counts': filename_without_extension}, inplace=True)

            if merged_df is None:
                merged_df = df
            else:
                merged_df = pd.merge(merged_df, df, on='tRF sequence', how='outer')

            sample_names.append(filename_without_extension)

        except Exception as e:
            print(f"Error reading {filename}: {e}")

    if merged_df is not None:
        merged_df.fillna(0, inplace=True)
        merged_df.to_csv(output_csv, index=True)  # Export to CSV with index
        
        # Generate colData
        # col_data = pd.DataFrame({
        #     'sample': sample_names,
        #     # Example condition; adjust based on your data
        #     'condition': ['control' if 'control' in name else 'treatment' for name in sample_names]
        # })
        # col_data.to_csv("/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/colData.csv", index=False)
        
        print(f"Merged DataFrame exported to {output_csv}")
        print("ColData exported to colData.csv")
    else:
        print("No DataFrames were merged.")

if __name__ == "__main__":
    folder_path = "/mnt/d/bioinformatics/MINT/MINTmap-release-v1.0/MINTmap-release-v1.0/outputs/CLL/chosen"  
    read_and_merge_dataframes(folder_path, count_threshold= 0, output_csv="/mnt/d/bioinformatics/RNAseq_pipeline/diffexp/merged_tRF_data.csv")
